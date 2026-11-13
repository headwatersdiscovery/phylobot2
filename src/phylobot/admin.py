from typing import Any

from django.contrib import admin
from django.contrib.auth import get_user_model
from django.contrib.auth.admin import UserAdmin

from phylobot.models import (
    Alignment,
    Job,
    JobAttempt,
    Library,
    LibraryExport,
    Tag,
    UserLimits,
)


class UserLimitsInline(admin.StackedInline):  # type: ignore[type-arg]
    model = UserLimits
    can_delete = False
    extra = 1
    max_num = 1
    verbose_name_plural = "Limits (blank = site default)"


User = get_user_model()
admin.site.unregister(User)


@admin.register(User)
class PhyloBotUserAdmin(UserAdmin):  # type: ignore[type-arg]
    inlines = [UserLimitsInline]


class JobAttemptInline(admin.TabularInline):  # type: ignore[type-arg]
    """Each run of the job: the record for billing, never edited."""

    model = JobAttempt
    extra = 0
    can_delete = False
    readonly_fields = [
        "number", "worker", "started_at", "finished_at", "status", "error",
        "completed_cpu_seconds", "failed_cpu_seconds", "tasks_run", "tasks_reused",
        "failed_rate_percent", "billable_cpu_seconds",
    ]  # fmt: skip

    def has_add_permission(self, request: Any, obj: Any = None) -> bool:
        return False


@admin.register(Job)
class JobAdmin(admin.ModelAdmin):  # type: ignore[type-arg]
    inlines = [JobAttemptInline]
    list_display = [
        "id", "name", "owner", "status", "attempt", "created_at", "finished_at",
        "cpu_seconds", "billable_cpu_seconds", "deleted_at",
    ]  # fmt: skip
    list_filter = [
        "status", ("deleted_at", admin.EmptyFieldListFilter),
        ("imported_at", admin.EmptyFieldListFilter), "tags",
    ]  # fmt: skip
    filter_horizontal = ["tags"]
    search_fields = ["id", "name", "owner__username"]
    readonly_fields = [
        "id", "created_at", "started_at", "finished_at", "heartbeat_at", "cpu_seconds",
        "billable_cpu_seconds", "attempt", "deleted_at", "purged_at", "rerun_of", "extends",
        "imported_at", "imported_from",
    ]  # fmt: skip


@admin.register(Tag)
class TagAdmin(admin.ModelAdmin):  # type: ignore[type-arg]
    list_display = ["name", "owner"]
    search_fields = ["name", "owner__username"]


class LibraryExportInline(admin.TabularInline):  # type: ignore[type-arg]
    """Receipts are never edited."""

    model = LibraryExport
    extra = 0
    can_delete = False
    readonly_fields = [
        "number", "created_at", "member_count", "phylobot_version", "spec_version", "spec",
    ]  # fmt: skip

    def has_add_permission(self, request: Any, obj: Any = None) -> bool:
        return False


@admin.register(Library)
class LibraryAdmin(admin.ModelAdmin):  # type: ignore[type-arg]
    list_display = ["id", "name", "code", "owner", "job", "updated_at", "deleted_at"]
    list_filter = [("deleted_at", admin.EmptyFieldListFilter)]
    search_fields = ["id", "name", "code", "owner__username", "job__id"]
    readonly_fields = ["id", "job", "created_at", "updated_at", "deleted_at"]
    inlines = [LibraryExportInline]

    def has_delete_permission(self, request: Any, obj: Any = None) -> bool:
        return False  # never hard-deleted, as jobs


@admin.register(Alignment)
class AlignmentAdmin(admin.ModelAdmin):  # type: ignore[type-arg]
    list_display = ["id", "name", "origin", "owner", "job", "sequences", "columns", "created_at",
                    "deleted_at"]  # fmt: skip
    list_filter = ["origin", "seqtype", ("deleted_at", admin.EmptyFieldListFilter)]
    search_fields = ["id", "name", "owner__username", "job__id"]
    readonly_fields = ["id", "job", "job_key", "path", "derived_from", "edits", "created_at",
                       "updated_at", "deleted_at"]  # fmt: skip

    def has_delete_permission(self, request: Any, obj: Any = None) -> bool:
        return False  # deleting keeps the record (alignments.delete), as for jobs
