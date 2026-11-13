#!/usr/bin/env bash
# Report the CPUs of Docker's Linux VM and suggest PHYLOBOT_DOCKER_CPUSET for .env.
#   scripts/docker-cpus.sh          print the report (make docker-cpus)
#   scripts/docker-cpus.sh --check  quiet when .env is fine; warn when it leaves an
#                                   isolated CPU in use; fail when it names CPUs the VM
#                                   does not have (every task container would refuse to start)
#
# Why: OrbStack boots its VM with isolcpus=<last CPU>. The kernel never load-balances an
# isolated CPU, and emulated (amd64) IQ-TREE ends up on it, so every concurrent IQ-TREE
# task shares one core (each at 20% with five running). Other runtimes may isolate
# nothing, in which case the setting stays empty.
set -euo pipefail

check=false
[ "${1:-}" = "--check" ] && check=true

value=$(grep -E '^PHYLOBOT_DOCKER_CPUSET=' .env 2>/dev/null | tail -1 | cut -d= -f2- | tr -d '"'"'"' ' || true)

# Same base image as docker/*.Dockerfile, run natively: it only reads two kernel files.
docker run --rm -i python:3.14.7-slim python - "$check" "$value" <<'EOF'
import sys

check, value = sys.argv[1] == "true", sys.argv[2]


def cpus(text):
    out = set()
    for part in filter(None, text.strip().split(",")):
        low, _, high = part.partition("-")
        out.update(range(int(low), int(high or low) + 1))
    return out


def show(numbers):
    numbers, runs = sorted(numbers), []
    for n in numbers:
        if runs and n == runs[-1][1] + 1:
            runs[-1][1] = n
        else:
            runs.append([n, n])
    return ",".join(f"{a}" if a == b else f"{a}-{b}" for a, b in runs)


online = cpus(open("/sys/devices/system/cpu/online").read())
isolated = cpus(open("/sys/devices/system/cpu/isolated").read())
suggested = show(online - isolated) if isolated else ""
chosen = cpus(value) if value else online

if check:
    if missing := chosen - online:
        print(f"Error: PHYLOBOT_DOCKER_CPUSET={value} names CPUs Docker's VM does not have "
              f"({show(missing)}; it has {show(online)}). Every pipeline step would fail to "
              f"start. Set it to {suggested or 'nothing (empty)'}: see make docker-cpus.",
              file=sys.stderr)
        sys.exit(1)
    if stuck := chosen & isolated:
        setting = f"PHYLOBOT_DOCKER_CPUSET={value} includes it" if value else (
            "PHYLOBOT_DOCKER_CPUSET is empty")
        print(f"Warning: Docker's VM isolates CPU {show(stuck)} and {setting}: IQ-TREE steps may all "
              f"share that one core. Set PHYLOBOT_DOCKER_CPUSET={suggested} in .env "
              f"(see make docker-cpus).", file=sys.stderr)
    sys.exit(0)

print(f"Docker VM CPUs online:   {show(online)}")
print(f"Isolated (isolcpus):     {show(isolated) or 'none'}")
print(f"PHYLOBOT_DOCKER_CPUSET:  {value or '(empty: all CPUs)'}")
if isolated:
    print(f"Suggested:               PHYLOBOT_DOCKER_CPUSET={suggested}")
else:
    print("Suggested:               leave it empty (no CPU is isolated)")
print("Change .env, then restart the worker when no job runs (make up). Recheck after "
      "changing the CPU count in your Docker runtime's settings.")
EOF
