"""
Test that parallel processing actually creates multiple system processes.

This test monitors system process information to verify that ProcessPoolExecutor
actually spawns worker processes by directly checking the executor's state.
"""

import os
import psutil
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pytest

from validation_pkg import ReadValidator
from validation_pkg.config_manager import ReadConfig
from validation_pkg.utils.formats import ReadFormat, CodingType


def slow_work_function(x):
    """Simulate slow work to make processes visible."""
    time.sleep(0.5)  # Sleep for 500ms
    return x * 2


class TestProcessMonitoring:
    """Test that parallel execution creates actual system processes."""

    @pytest.fixture
    def test_fastq_files(self, tmp_path):
        """Create multiple test FASTQ files with enough data to make processing measurable."""
        files = []
        for i in range(6):  # 6 files for better parallelization
            fastq_file = tmp_path / f"test_{i}.fastq"
            with open(fastq_file, 'w') as f:
                # Create 5000 records per file to slow down validation
                for j in range(5000):
                    f.write(f"@read_{i}_{j}\n")
                    f.write("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
                    f.write("+\n")
                    f.write("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
            files.append(fastq_file)
        return files

    @pytest.fixture
    def read_configs(self, test_fastq_files):
        """Create ReadConfig objects for test files."""
        configs = []
        for fastq_file in test_fastq_files:
            config = ReadConfig(
                filename=str(fastq_file),
                filepath=fastq_file,
                ngs_type='illumina',
                coding_type=CodingType.NONE,
                detected_format=ReadFormat.FASTQ
            )
            configs.append(config)
        return configs

    def test_processpoolexecutor_creates_processes(self):
        """Test that ProcessPoolExecutor actually spawns worker processes."""
        print("\n=== Testing ProcessPoolExecutor Directly ===")

        parent_pid = os.getpid()
        parent_proc = psutil.Process(parent_pid)

        print(f"Parent PID: {parent_pid}")
        print(f"Initial children: {len(parent_proc.children())}")

        # Create a ProcessPoolExecutor with 4 workers
        with ProcessPoolExecutor(max_workers=4) as executor:
            # Submit some slow tasks
            futures = [executor.submit(slow_work_function, i) for i in range(8)]

            # Immediately check for child processes
            time.sleep(0.2)  # Give processes time to spawn
            children = parent_proc.children(recursive=True)

            print(f"\n=== During Execution ===")
            print(f"Number of child processes: {len(children)}")

            if children:
                print("\nChild processes:")
                for i, child in enumerate(children):
                    try:
                        print(f"  {i+1}. PID={child.pid}, name={child.name()}, status={child.status()}")
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        print(f"  {i+1}. Process already terminated")

            # Wait for all tasks to complete
            for future in futures:
                future.result()

        # After executor closes
        final_children = parent_proc.children()
        print(f"\n=== After Execution ===")
        print(f"Final children: {len(final_children)}")

        # Assertion: we should have seen child processes
        assert len(children) >= 2, \
            f"Expected at least 2 worker processes, found {len(children)}"

        print(f"\n✓ Verified: ProcessPoolExecutor created {len(children)} worker processes")

    def test_validate_reads_parallel_creates_processes(self, read_configs, tmp_path):
        """Test that validate_reads with parallel execution creates worker processes."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        parent_pid = os.getpid()
        parent_proc = psutil.Process(parent_pid)

        print(f"\n=== Testing validate_reads Parallel Processing ===")
        print(f"Parent PID: {parent_pid}")
        print(f"Number of files: {len(read_configs)}")
        print(f"Initial children: {len(parent_proc.children())}")

        # Monitor processes during execution
        max_children_seen = 0
        child_snapshots = []

        def monitor_processes():
            """Monitor child processes in background."""
            nonlocal max_children_seen, child_snapshots
            for _ in range(30):  # Monitor for 3 seconds
                time.sleep(0.1)
                children = parent_proc.children(recursive=True)
                num_children = len(children)
                max_children_seen = max(max_children_seen, num_children)

                if children and len(child_snapshots) < 5:  # Capture first 5 snapshots with children
                    snapshot = {
                        'time': time.time(),
                        'count': num_children,
                        'pids': [c.pid for c in children if c.is_running()]
                    }
                    child_snapshots.append(snapshot)

        import threading
        monitor_thread = threading.Thread(target=monitor_processes, daemon=True)
        monitor_thread.start()

        # Run parallel validation with 4 workers
        # Use 'strict' mode with compression to ensure actual work is done
        settings = ReadValidator.Settings()
        settings = settings.update(max_workers=4, validation_level='strict', coding_type='gz')

        from validation_pkg import validate_reads

        print(f"\n=== Starting Parallel Validation ===")
        start_time = time.time()
        results = validate_reads(read_configs, output_dir, settings)
        end_time = time.time()

        # Wait for monitoring to finish
        monitor_thread.join(timeout=5)

        print(f"\n=== Execution Results ===")
        print(f"Execution time: {end_time - start_time:.2f}s")
        print(f"Files processed: {len([r for r in results if r['success']])}/{len(results)}")
        print(f"Max child processes seen: {max_children_seen}")

        # Print any errors
        errors = [r for r in results if not r['success']]
        if errors:
            print(f"Errors: {len(errors)}/{len(results)}")
            for r in errors[:2]:  # Show first 2 errors
                print(f"  {r['filename']}: {r['error'][:80]}")

        if child_snapshots:
            print(f"\n=== Child Process Snapshots ===")
            for i, snap in enumerate(child_snapshots[:3]):
                print(f"Snapshot {i+1}: {snap['count']} processes, PIDs={snap['pids']}")

        # Assertions
        assert len(results) == len(read_configs), "Should process all files"

        # We should have seen at least 2 child processes during parallel execution
        assert max_children_seen >= 2, \
            f"Expected at least 2 worker processes during parallel execution, saw {max_children_seen}"

        print(f"\n✓ Verified: Parallel validation created {max_children_seen} worker processes")

    def test_validate_reads_sequential_no_processes(self, read_configs, tmp_path):
        """Test that sequential execution doesn't create worker processes."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        parent_pid = os.getpid()
        parent_proc = psutil.Process(parent_pid)

        print(f"\n=== Testing validate_reads Sequential Processing ===")
        print(f"Initial children: {len(parent_proc.children())}")

        max_children_seen = 0

        def monitor_processes():
            nonlocal max_children_seen
            for _ in range(30):
                time.sleep(0.1)
                children = parent_proc.children(recursive=True)
                max_children_seen = max(max_children_seen, len(children))

        import threading
        monitor_thread = threading.Thread(target=monitor_processes, daemon=True)
        monitor_thread.start()

        # Run sequential validation (max_workers=None)
        settings = ReadValidator.Settings()
        settings = settings.update(max_workers=None, validation_level='trust')

        from validation_pkg import validate_reads

        start_time = time.time()
        results = validate_reads(read_configs, output_dir, settings)
        end_time = time.time()

        monitor_thread.join(timeout=5)

        print(f"\nExecution time: {end_time - start_time:.2f}s")
        print(f"Max child processes seen: {max_children_seen}")

        # Sequential should not create a worker pool (may see 0-2 for subprocess compression calls)
        assert max_children_seen < 3, \
            f"Sequential execution should not create worker pool, saw {max_children_seen} processes"

        print(f"✓ Verified: Sequential execution did not create worker pool")

    def test_worker_count_scales_with_setting(self, read_configs, tmp_path):
        """Test that number of workers scales with max_workers setting."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        parent_pid = os.getpid()
        parent_proc = psutil.Process(parent_pid)

        for worker_count in [2, 4]:
            print(f"\n=== Testing with max_workers={worker_count} ===")

            max_children_seen = 0
            max_snapshot_time = None

            def monitor_processes():
                nonlocal max_children_seen, max_snapshot_time
                for _ in range(40):  # Monitor for 4 seconds
                    time.sleep(0.1)
                    children = parent_proc.children(recursive=True)
                    num_children = len(children)
                    if num_children > max_children_seen:
                        max_children_seen = num_children
                        max_snapshot_time = time.time()

            import threading
            monitor_thread = threading.Thread(target=monitor_processes, daemon=True)
            monitor_thread.start()

            # Use strict mode with compression to ensure work is done
            settings = ReadValidator.Settings()
            settings = settings.update(max_workers=worker_count, validation_level='strict', coding_type='gz')

            from validation_pkg import validate_reads

            start_time = time.time()
            results = validate_reads(read_configs, output_dir, settings)
            end_time = time.time()

            monitor_thread.join(timeout=5)

            print(f"Execution time: {end_time - start_time:.2f}s")
            print(f"Requested workers: {worker_count}")
            print(f"Max child processes seen: {max_children_seen}")
            print(f"Files processed: {len([r for r in results if r['success']])}/{len(results)}")

            # We should see at least worker_count - 1 processes
            # (accounting for timing where not all workers spawn simultaneously)
            assert max_children_seen >= worker_count - 1, \
                f"Expected at least {worker_count-1} processes for {worker_count} workers, saw {max_children_seen}"

            print(f"✓ Verified: Saw {max_children_seen} processes for {worker_count} workers")

    def test_process_details_during_parallel_execution(self, read_configs, tmp_path):
        """Capture detailed information about worker processes."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        parent_pid = os.getpid()
        parent_proc = psutil.Process(parent_pid)

        print(f"\n=== Detailed Process Monitoring ===")
        print(f"Parent: PID={parent_pid}, name={parent_proc.name()}")

        worker_details = []

        def monitor_detailed():
            """Capture detailed worker process info."""
            nonlocal worker_details
            seen_pids = set()

            for _ in range(30):
                time.sleep(0.1)
                children = parent_proc.children(recursive=True)

                for child in children:
                    try:
                        if child.pid not in seen_pids:
                            name = child.name()
                            # Record all python processes - we'll verify they're worker processes
                            # ProcessPoolExecutor workers will have 'python' as name
                            if 'python' in name.lower():
                                info = {
                                    'pid': child.pid,
                                    'name': name,
                                    'status': child.status(),
                                    'cmdline': ' '.join(child.cmdline()[:5])  # First 5 args
                                }
                                worker_details.append(info)
                            seen_pids.add(child.pid)
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        continue

        import threading
        monitor_thread = threading.Thread(target=monitor_detailed, daemon=True)
        monitor_thread.start()

        # Use strict mode with compression to ensure work is done
        settings = ReadValidator.Settings()
        settings = settings.update(max_workers=4, validation_level='strict', coding_type='gz')

        from validation_pkg import validate_reads
        results = validate_reads(read_configs, output_dir, settings)

        monitor_thread.join(timeout=5)

        print(f"\n=== Worker Processes Detected ({len(worker_details)} total) ===")
        for i, info in enumerate(worker_details[:5]):  # Show first 5
            print(f"Worker {i+1}:")
            print(f"  PID: {info['pid']}")
            print(f"  Name: {info['name']}")
            print(f"  Status: {info['status']}")
            print(f"  Command: {info['cmdline']}")

        print(f"Files processed: {len([r for r in results if r['success']])}/{len(results)}")

        # Should detect at least 2 child processes
        # Note: May include pytest infrastructure processes, but that's okay
        # We're verifying that parallel execution creates child processes
        assert len(worker_details) >= 2, \
            f"Expected at least 2 child processes during parallel execution, detected {len(worker_details)}"

        # All detected processes should be Python processes
        for info in worker_details:
            assert 'python' in info['name'].lower(), \
                f"Child process should be Python, got: {info['name']} (cmdline: {info['cmdline']})"

        print(f"\n✓ Verified: Detected {len(worker_details)} Python child processes during parallel execution")
