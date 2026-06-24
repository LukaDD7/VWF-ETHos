#!/usr/bin/env python3
"""Resilient multi-GPU scheduler for resumable 7A6O production MD."""
from __future__ import annotations

import argparse
import csv
import os
import signal
import subprocess
import time
from pathlib import Path


def log(message: str) -> None:
    print(f"[{time.strftime('%F %T')}] {message}", flush=True)


def md_done(out_root: Path, variant: str) -> bool:
    work = out_root / variant / "md_7a6o"
    return (work / "md_prod.gro").is_file() or any(work.glob("md_prod.part*.gro"))


def gpu_free_mib() -> dict[int, int]:
    result = subprocess.run(
        ["nvidia-smi", "--query-gpu=index,memory.free", "--format=csv,noheader,nounits"],
        check=True, capture_output=True, text=True,
    )
    free = {}
    for row in csv.reader(result.stdout.splitlines()):
        free[int(row[0].strip())] = int(row[1].strip())
    return free


def read_variants(path: Path) -> list[str]:
    variants = []
    for line in path.read_text().splitlines():
        clean = line.split("#", 1)[0].strip()
        if clean:
            variants.append(clean.split()[0])
    return variants


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--variants-file", type=Path, required=True)
    parser.add_argument("--gpu-ids", default="1,2,3,4,5,6")
    parser.add_argument("--min-free-mib", type=int, default=12000)
    parser.add_argument("--poll-seconds", type=int, default=30)
    parser.add_argument("--retry-seconds", type=int, default=300)
    parser.add_argument("--ns", type=int, default=50)
    parser.add_argument("--nvt-ps", type=int, default=50)
    parser.add_argument("--npt-ps", type=int, default=200)
    parser.add_argument("--ntomp", type=int, default=8)
    args = parser.parse_args()

    root = args.root.resolve()
    out_root = root / "output/gromacs_md_autoinhib"
    scheduler_pid = out_root / "resilient_md_scheduler.pid"
    if scheduler_pid.is_file():
        try:
            existing_pid = int(scheduler_pid.read_text().strip())
            os.kill(existing_pid, 0)
        except (ValueError, ProcessLookupError, PermissionError):
            pass
        else:
            log(f"FATAL scheduler already running: PID {existing_pid}")
            return 3
    scheduler_pid.write_text(f"{os.getpid()}\n")

    def clear_pid_file(signum=None, _frame=None) -> None:
        try:
            if scheduler_pid.read_text().strip() == str(os.getpid()):
                scheduler_pid.unlink()
        except FileNotFoundError:
            pass
        if signum is not None:
            raise SystemExit(128 + signum)

    signal.signal(signal.SIGTERM, clear_pid_file)
    signal.signal(signal.SIGINT, clear_pid_file)
    mutant_dir = root / "structures/7a6o_mutants"
    runner = root / "scripts/pipeline/run_7a6o_variant_direct.sh"
    gpu_ids = [int(item) for item in args.gpu_ids.split(",") if item.strip()]

    pending = []
    for variant in read_variants(args.variants_file):
        if not (mutant_dir / f"{variant}.pdb").is_file():
            log(f"skip non-runnable variant: {variant}")
        elif md_done(out_root, variant):
            log(f"skip MD done: {variant}")
        elif not (out_root / variant / "relax_pdb/solv_ions_em.gro").is_file():
            log(f"FATAL missing relaxed structure: {variant}")
            return 2
        else:
            pending.append(variant)

    active: dict[int, tuple[subprocess.Popen, str, object]] = {}
    retry_after: dict[str, float] = {}
    attempts: dict[str, int] = {}
    log(f"resilient MD start: pending={len(pending)} GPUs={gpu_ids} min_free={args.min_free_mib} MiB")

    while pending or active:
        for gpu, (process, variant, handle) in list(active.items()):
            rc = process.poll()
            if rc is None:
                continue
            handle.close()
            del active[gpu]
            if rc == 0 and md_done(out_root, variant):
                log(f"complete {variant} on GPU {gpu}")
            else:
                retry_after[variant] = time.time() + args.retry_seconds
                pending.append(variant)
                log(f"retry {variant}: rc={rc}, checkpoint retained, cooldown={args.retry_seconds}s")

        if pending:
            try:
                free = gpu_free_mib()
            except Exception as exc:
                log(f"nvidia-smi unavailable: {exc}; retrying")
                time.sleep(args.poll_seconds)
                continue

            for gpu in gpu_ids:
                if gpu in active or free.get(gpu, 0) < args.min_free_mib:
                    continue
                ready_index = next(
                    (i for i, variant in enumerate(pending)
                     if retry_after.get(variant, 0) <= time.time()),
                    None,
                )
                if ready_index is None:
                    break
                variant = pending.pop(ready_index)
                attempts[variant] = attempts.get(variant, 0) + 1
                log_path = out_root / f"{variant}_md_direct.log"
                handle = log_path.open("a")
                env = os.environ.copy()
                env.update({
                    "NS": str(args.ns), "NVT_PS": str(args.nvt_ps),
                    "NPT_PS": str(args.npt_ps), "NTOMP": str(args.ntomp),
                })
                handle.write(f"\n=== resilient attempt {attempts[variant]} GPU {gpu} {time.strftime('%F %T')} ===\n")
                handle.flush()
                process = subprocess.Popen(
                    ["bash", str(runner), variant, str(gpu)], cwd=root,
                    env=env, stdout=handle, stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
                active[gpu] = (process, variant, handle)
                log(f"launch {variant} attempt={attempts[variant]} GPU={gpu} free={free[gpu]} MiB")

        if pending or active:
            time.sleep(args.poll_seconds)

    log("all runnable variants have completed MD")
    clear_pid_file()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
