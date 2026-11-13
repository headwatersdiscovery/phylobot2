"""Object storage (Garage locally, S3 in production), keyed under jobs/<id>/."""

import time
from collections.abc import Callable, Iterator
from functools import cache
from pathlib import Path
from typing import Any

import boto3
from botocore.config import Config
from botocore.exceptions import ClientError
from django.conf import settings

# Files over 8 MB go up in parts (upload_file), and Garage stores a checksum of the parts'
# checksums without marking it as such; botocore then compares it with the whole file's
# and refuses every read of a perfectly good file (FlexibleChecksumError; found
# 2026-10-19 on an 11.7 MB ancestors_posteriors.tsv). Validate responses only where an
# operation requires it, as boto3 did before 1.36.
CONFIG = Config(response_checksum_validation="when_required")


@cache
def client() -> Any:
    """One client per process: making one costs ~6 ms, which a page reading 300 small
    files paid 300 times. boto3 clients are safe to share between threads."""
    return boto3.client("s3", endpoint_url=settings.S3_ENDPOINT_URL, config=CONFIG)


def put_bytes(key: str, data: bytes) -> None:
    client().put_object(Bucket=settings.S3_BUCKET, Key=key, Body=data)


def put_file(key: str, path: Path) -> None:
    client().upload_file(str(path), settings.S3_BUCKET, key)


def get_file(key: str, path: Path) -> None:
    client().download_file(settings.S3_BUCKET, key, str(path))


def delete_key(key: str) -> None:
    client().delete_object(Bucket=settings.S3_BUCKET, Key=key)


def list_keys(prefix: str) -> list[str]:
    keys: list[str] = []
    paginator = client().get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=settings.S3_BUCKET, Prefix=prefix):
        keys += [item["Key"] for item in page.get("Contents", [])]
    return keys


def list_sizes(prefix: str) -> dict[str, int]:
    """{key: bytes} for every object under `prefix`."""
    sizes: dict[str, int] = {}
    paginator = client().get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=settings.S3_BUCKET, Prefix=prefix):
        sizes |= {item["Key"]: item["Size"] for item in page.get("Contents", [])}
    return sizes


def exists(key: str) -> bool:
    try:
        client().head_object(Bucket=settings.S3_BUCKET, Key=key)
    except ClientError as error:
        if error.response.get("Error", {}).get("Code") in ("404", "NoSuchKey", "NotFound"):
            return False
        raise
    return True


def prefix_size(prefix: str) -> int:
    """Total bytes stored under `prefix` (one listing request per 1000 objects)."""
    size = 0
    paginator = client().get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=settings.S3_BUCKET, Prefix=prefix):
        size += sum(item["Size"] for item in page.get("Contents", []))
    return size


def download_prefix(prefix: str, dest: Path) -> list[Path]:
    """Download every object under `prefix` into `dest`, keeping relative paths."""
    paths = []
    s3 = client()
    for key in list_keys(prefix):
        relative = key[len(prefix) :]
        target = (dest / relative).resolve()
        if not target.is_relative_to(dest.resolve()):
            raise ValueError(f"refusing to write outside {dest}: {key}")
        target.parent.mkdir(parents=True, exist_ok=True)
        s3.download_file(settings.S3_BUCKET, key, str(target))
        paths.append(target)
    return paths


def upload_tree(
    directory: Path,
    prefix: str,
    alive: Callable[[], bool] | None = None,
    every: float = 30.0,
) -> int:
    """Upload every file under `directory` to `prefix`, keeping relative paths. `alive`
    is called about every `every` seconds (a job's heartbeat: a large job's upload took
    1.5 h with none); when it returns False the upload stops, raising UploadStopped."""
    count = 0
    last = time.monotonic()
    for path in sorted(directory.rglob("*")):
        if path.is_file():
            put_file(prefix + path.relative_to(directory).as_posix(), path)
            count += 1
            if alive is not None and time.monotonic() - last >= every:
                last = time.monotonic()
                if not alive():
                    raise UploadStopped(f"stopped after {count} files")
    return count


class UploadStopped(Exception):
    """upload_tree's `alive` said to stop."""


def copy_key(source: str, dest: str) -> None:
    """Copy one object (server-side): the bytes never pass through this process, which
    matters for posterior tables that run to megabytes."""
    client().copy_object(
        Bucket=settings.S3_BUCKET,
        Key=dest,
        CopySource={"Bucket": settings.S3_BUCKET, "Key": source},
    )


def copy_prefix(source: str, dest: str) -> int:
    """Copy every object under `source` to `dest` (server-side), keeping relative paths."""
    s3 = client()
    keys = list_keys(source)
    for key in keys:
        s3.copy_object(
            Bucket=settings.S3_BUCKET,
            Key=dest + key[len(source) :],
            CopySource={"Bucket": settings.S3_BUCKET, "Key": key},
        )
    return len(keys)


def stream(key: str, chunk_size: int = 64 * 1024) -> Iterator[bytes]:
    body = client().get_object(Bucket=settings.S3_BUCKET, Key=key)["Body"]
    yield from body.iter_chunks(chunk_size)


def size(key: str) -> int:
    return int(client().head_object(Bucket=settings.S3_BUCKET, Key=key)["ContentLength"])


def read_range(key: str, start: int, length: int) -> bytes:
    """`length` bytes of an object from `start` (one ranged GET)."""
    if length <= 0:
        return b""
    found = client().get_object(
        Bucket=settings.S3_BUCKET, Key=key, Range=f"bytes={start}-{start + length - 1}"
    )
    return bytes(found["Body"].read())


def read_bytes(key: str) -> bytes:
    return b"".join(stream(key))


def delete_prefix(prefix: str) -> int:
    """Delete every object under `prefix`; returns how many."""
    if not prefix.endswith("/") or prefix.count("/") < 2:
        raise ValueError(f"refusing to delete a broad prefix: {prefix!r}")
    keys = list_keys(prefix)
    s3 = client()
    for start in range(0, len(keys), 1000):
        batch = [{"Key": key} for key in keys[start : start + 1000]]
        s3.delete_objects(Bucket=settings.S3_BUCKET, Delete={"Objects": batch, "Quiet": True})
    return len(keys)


def bucket_reachable() -> bool:
    try:
        client().head_bucket(Bucket=settings.S3_BUCKET)
    except Exception:
        return False
    return True
