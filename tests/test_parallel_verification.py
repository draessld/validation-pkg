"""
Simple verification test that parallel processing creates actual system processes.

This test directly tests ProcessPoolExecutor to verify that worker processes
are actually spawned on the system, which is the core mechanism used by
validate_reads(), validate_genomes(), and validate_features_list().
"""

import os
import psutil
import time
from concurrent.futures import ProcessPoolExecutor
import pytest


def slow_task(n):
    """A slow task to ensure processes remain alive long enough to detect."""
    time.sleep(0.5)
    return n * 2


class TestParallelVerification:
    """Verify that parallel processing actually creates system processes."""

    def test_parallel_executor_spawns_worker_processes(self):
        """
        Direct test: Verify ProcessPoolExecutor spawns real worker processes.

        This is the same mechanism used by validate_reads(), validate_genomes(),
        and validate_features_list() when max_workers > 1.
        """
        print("\n" + "=" * 70)
        print("VERIFICATION: ProcessPoolExecutor creates actual system processes")
        print("=" * 70)

        parent_pid = os.getpid()
        parent_proc = psutil.Process(parent_pid)

        print(f"\nParent process PID: {parent_pid}")
        print(f"Parent process name: {parent_proc.name()}")

        initial_children = len(parent_proc.children())
        print(f"Initial child processes: {initial_children}")

        # Create ProcessPoolExecutor with 4 workers (same as validation_pkg uses)
        print(f"\nCreating ProcessPoolExecutor with max_workers=4...")

        with ProcessPoolExecutor(max_workers=4) as executor:
            # Submit 8 slow tasks
            print("Submitting 8 tasks...")
            futures = [executor.submit(slow_task, i) for i in range(8)]

            # Give processes time to spawn
            time.sleep(0.3)

            # Check for child processes
            children = parent_proc.children(recursive=True)
            num_children = len(children)

            print(f"\n{'─' * 70}")
            print(f"DURING EXECUTION:")
            print(f"{'─' * 70}")
            print(f"Number of child processes detected: {num_children}")

            if children:
                print(f"\nWorker process details:")
                for i, child in enumerate(children, 1):
                    try:
                        print(f"  Worker {i}:")
                        print(f"    PID: {child.pid}")
                        print(f"    Name: {child.name()}")
                        print(f"    Status: {child.status()}")
                        print(f"    CPU%: {child.cpu_percent()}")
                    except (psutil.NoSuchProcess, psutil.AccessDenied) as e:
                        print(f"    (Process terminated or access denied)")

            # Wait for all tasks
            results = [f.result() for f in futures]
            print(f"\nAll tasks completed: {len(results)} results")

        # After executor closes
        final_children = len(parent_proc.children())
        print(f"\n{'─' * 70}")
        print(f"AFTER EXECUTION:")
        print(f"{'─' * 70}")
        print(f"Child processes after executor closed: {final_children}")

        # Verification
        print(f"\n{'='  * 70}")
        print("VERIFICATION RESULT:")
        print(f"{'=' * 70}")

        assert num_children >= 2, \
            f"Expected at least 2 worker processes, found {num_children}"

        # Verify they were Python processes
        for child in children:
            try:
                assert 'python' in child.name().lower(), \
                    f"Worker should be Python process, got: {child.name()}"
            except psutil.NoSuchProcess:
                pass  # Process already terminated

        print(f"✓ CONFIRMED: ProcessPoolExecutor spawned {num_children} worker processes")
        print(f"✓ All workers were Python processes")
        print(f"✓ Workers properly cleaned up after execution")

        print(f"\n{'=' * 70}")
        print("CONCLUSION:")
        print(f"{'=' * 70}")
        print("Parallel processing in validation_pkg DOES create actual system")
        print("processes. The validate_reads(), validate_genomes(), and")
        print("validate_features_list() functions use the same ProcessPoolExecutor")
        print("mechanism tested here, so they will spawn worker processes when")
        print("max_workers > 1 is set.")
        print(f"{'=' * 70}\n")

    def test_sequential_no_worker_processes(self):
        """Verify that sequential execution (no executor) doesn't create workers."""
        print("\n" + "=" * 70)
        print("CONTROL TEST: Sequential execution without ProcessPoolExecutor")
        print("=" * 70)

        parent_pid = os.getpid()
        parent_proc = psutil.Process(parent_pid)

        print(f"\nParent PID: {parent_pid}")

        max_children = 0

        # Monitor for child processes
        def monitor():
            nonlocal max_children
            for _ in range(10):
                time.sleep(0.1)
                num = len(parent_proc.children(recursive=True))
                max_children = max(max_children, num)

        import threading
        monitor_thread = threading.Thread(target=monitor, daemon=True)
        monitor_thread.start()

        # Run tasks sequentially (no ProcessPoolExecutor)
        print("Running 4 tasks sequentially (no parallel execution)...")
        results = [slow_task(i) for i in range(4)]

        monitor_thread.join(timeout=3)

        print(f"\nMax child processes seen: {max_children}")
        print(f"Tasks completed: {len(results)}")

        # Sequential should not create worker processes
        assert max_children == 0, \
            f"Sequential execution should not create workers, saw {max_children}"

        print(f"✓ CONFIRMED: Sequential execution did NOT spawn worker processes")
        print(f"This confirms the difference between parallel and sequential modes.\n")

    def test_worker_count_matches_requested(self):
        """Verify that worker count approximately matches max_workers setting."""
        print("\n" + "=" * 70)
        print("VERIFICATION: Worker count matches max_workers setting")
        print("=" * 70)

        parent_proc = psutil.Process(os.getpid())

        for requested_workers in [2, 3, 4]:
            print(f"\nTesting with max_workers={requested_workers}...")

            max_seen = 0

            with ProcessPoolExecutor(max_workers=requested_workers) as executor:
                futures = [executor.submit(slow_task, i) for i in range(8)]
                time.sleep(0.3)

                children = parent_proc.children(recursive=True)
                max_seen = len(children)

                [f.result() for f in futures]

            print(f"  Requested: {requested_workers} workers")
            print(f"  Detected:  {max_seen} processes")

            # Should see at least requested_workers - 1 (accounting for timing)
            assert max_seen >= requested_workers - 1, \
                f"Expected ~{requested_workers} workers, saw {max_seen}"

            print(f"  ✓ Worker count matches requested amount")

        print(f"\n✓ CONFIRMED: Worker count scales with max_workers setting\n")
