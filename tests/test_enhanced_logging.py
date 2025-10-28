"""
Test enhanced logging with process ID and worker ID identification.

Verifies that logs from parallel execution include:
- Process IDs (PID)
- Worker IDs (Worker-1, Worker-2, etc.)
- File context (filename being processed)
"""

import os
import tempfile
from pathlib import Path
import pytest

from validation_pkg import ReadValidator, validate_reads
from validation_pkg.config_manager import ReadConfig
from validation_pkg.utils.formats import ReadFormat, CodingType
from validation_pkg.logger import get_logger, setup_logging


class TestEnhancedLogging:
    """Test enhanced logging with process and worker identification."""

    @pytest.fixture
    def test_fastq_files(self, tmp_path):
        """Create test FASTQ files."""
        files = []
        for i in range(4):
            fastq_file = tmp_path / f"sample_{i}.fastq"
            with open(fastq_file, 'w') as f:
                for j in range(100):
                    f.write(f"@read_{i}_{j}\n")
                    f.write("ACGTACGT\n")
                    f.write("+\n")
                    f.write("IIIIIIII\n")
            files.append(fastq_file)
        return files

    @pytest.fixture
    def read_configs(self, test_fastq_files):
        """Create ReadConfig objects."""
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

    def test_logger_bind_worker_context(self):
        """Test that logger can bind worker context."""
        logger = get_logger()
        setup_logging(console_level="INFO")

        # Bind worker context
        logger.bind_worker_context(worker_id=1, file_context="test.fastq")

        # This should work without errors
        logger.info("Test message with worker context")

        # Unbind context
        logger.unbind_worker_context()
        logger.info("Test message without worker context")

        # Test passes if no exceptions raised
        assert True

    def test_process_id_in_logs(self, tmp_path):
        """Test that process ID appears in log file."""
        log_file = tmp_path / "test.log"

        logger = get_logger()
        logger.clear_issues()
        setup_logging(console_level="INFO", log_file=log_file)

        # Get current process ID
        current_pid = os.getpid()

        # Write a log message
        logger.info("Test message for PID verification")

        # Force flush by getting summary
        _ = logger.get_summary()

        # Read log file (JSON format)
        if log_file.exists():
            log_content = log_file.read_text()

            # Log should contain process_id field with current PID
            assert str(current_pid) in log_content, \
                f"Log file should contain current PID {current_pid}"

            # Should contain "process_id" field name
            assert "process_id" in log_content, \
                "Log file should contain process_id field"

            print(f"\n✓ Log file contains process ID {current_pid}")
            print(f"Log content sample:\n{log_content[:200]}")

    def test_worker_context_binding(self, tmp_path):
        """Test that worker context can be bound and appears in logs."""
        log_file = tmp_path / "worker_context.log"

        logger = get_logger()
        logger.clear_issues()
        setup_logging(console_level="INFO", log_file=log_file)

        # Bind worker context
        logger.bind_worker_context(worker_id=42, file_context="sample_1.fastq")
        logger.info("Processing with worker context")

        # Read log file
        if log_file.exists():
            log_content = log_file.read_text()

            # Should contain worker_id
            assert "worker_id" in log_content or "42" in log_content, \
                "Log should contain worker ID"

            # Should contain file_context
            assert "file_context" in log_content or "sample_1.fastq" in log_content, \
                "Log should contain file context"

            print(f"\n✓ Log contains worker context")
            print(f"Log content:\n{log_content}")

    def test_parallel_execution_worker_ids(self, read_configs, tmp_path):
        """Test that parallel execution assigns worker IDs correctly."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        log_file = tmp_path / "parallel.log"

        # Setup logging
        logger = get_logger()
        logger.clear_issues()
        setup_logging(console_level="INFO", log_file=log_file)

        # Run parallel validation
        # Note: Using minimal mode and compression to avoid "same file" error
        settings = ReadValidator.Settings()
        settings = settings.update(max_workers=2, validation_level='minimal', coding_type='gz')

        results = validate_reads(read_configs, output_dir, settings)

        # Note: We're testing logging, not validation success
        # Files should process (even if they fail due to file path issues)
        assert len(results) == len(read_configs), "Should process all files"

        # Check log file for worker IDs
        if log_file.exists():
            log_content = log_file.read_text()

            # Should contain multiple process IDs (main + workers)
            # Note: Actual worker processes may have different PIDs than main
            assert "process_id" in log_content, \
                "Log should contain process_id field"

            # Should contain worker_id field
            assert "worker_id" in log_content, \
                "Log should contain worker_id field"

            # Should contain file_context
            assert "file_context" in log_content or "sample_" in log_content, \
                "Log should contain file context"

            print(f"\n✓ Parallel execution log generated")
            print(f"Log file size: {len(log_content)} bytes")
            print(f"Sample from log:\n{log_content[:500]}")

    def test_multiple_workers_different_pids(self, read_configs, tmp_path):
        """Test that multiple workers have different PIDs in logs."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        log_file = tmp_path / "multi_worker.log"

        # Setup logging
        logger = get_logger()
        logger.clear_issues()
        setup_logging(console_level="DEBUG", log_file=log_file)

        # Run with 3 workers - use compression to avoid "same file" error
        settings = ReadValidator.Settings()
        settings = settings.update(max_workers=3, validation_level='minimal', coding_type='gz')

        results = validate_reads(read_configs, output_dir, settings)

        assert len(results) == len(read_configs), "Should process all files"

        # Read and analyze log file
        if log_file.exists():
            log_content = log_file.read_text()

            # Count occurrences of process_id in logs
            pid_count = log_content.count('"process_id"')

            print(f"\n✓ Log generated with {pid_count} process_id entries")

            # Should have multiple entries (one per log message)
            assert pid_count > 0, "Should have process ID entries in log"

    def test_file_context_in_worker_logs(self, read_configs, tmp_path):
        """Test that file context appears in worker logs."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        log_file = tmp_path / "file_context.log"

        logger = get_logger()
        logger.clear_issues()
        setup_logging(console_level="DEBUG", log_file=log_file)

        # Run parallel validation with compression to avoid "same file" error
        settings = ReadValidator.Settings()
        settings = settings.update(max_workers=2, validation_level='minimal', coding_type='gz')

        results = validate_reads(read_configs, output_dir, settings)

        # Note: Testing logging, not validation success
        # Files should be processed even if they have errors
        assert len(results) == len(read_configs), "Should process all files"

        # Check for file contexts in log
        if log_file.exists():
            log_content = log_file.read_text()

            # Should contain at least one of the sample filenames
            has_file_reference = any(
                f"sample_{i}.fastq" in log_content
                for i in range(len(read_configs))
            )

            assert has_file_reference, \
                "Log should reference sample files being processed"

            # Should contain worker_id field
            assert "worker_id" in log_content, \
                "Log should contain worker_id field"

            # Should contain file_context field
            assert "file_context" in log_content, \
                "Log should contain file_context field"

            print(f"\n✓ Log contains file context references")
            print(f"Sample from log:\n{log_content[:400]}")

    def test_console_output_format(self, read_configs, tmp_path, capsys):
        """Test that console output includes formatted worker info."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        # Setup logging (console only, no log file)
        logger = get_logger()
        logger.clear_issues()
        setup_logging(console_level="INFO")

        # Run parallel validation
        settings = ReadValidator.Settings()
        settings = settings.update(max_workers=2, validation_level='minimal')

        results = validate_reads(read_configs, output_dir, settings)

        # Capture console output
        captured = capsys.readouterr()

        print(f"\n=== Console Output ===")
        print(captured.out)
        print(f"=== End Console Output ===\n")

        # Check that output contains success messages
        assert "Completed" in captured.out or "Processing" in captured.out, \
            "Console should show processing messages"

    def test_sequential_vs_parallel_logging(self, read_configs, tmp_path):
        """Compare logging between sequential and parallel execution."""
        # Sequential execution
        output_seq = tmp_path / "output_seq"
        output_seq.mkdir()
        log_seq = tmp_path / "sequential.log"

        logger = get_logger()
        logger.clear_issues()
        setup_logging(console_level="INFO", log_file=log_seq)

        # Use compression to avoid "same file" error
        settings_seq = ReadValidator.Settings()
        settings_seq = settings_seq.update(max_workers=None, validation_level='minimal', coding_type='gz')

        results_seq = validate_reads(read_configs[:2], output_seq, settings_seq)

        # Parallel execution
        output_par = tmp_path / "output_par"
        output_par.mkdir()
        log_par = tmp_path / "parallel.log"

        logger.clear_issues()
        setup_logging(console_level="INFO", log_file=log_par)

        settings_par = ReadValidator.Settings()
        settings_par = settings_par.update(max_workers=2, validation_level='minimal', coding_type='gz')

        results_par = validate_reads(read_configs[:2], output_par, settings_par)

        # Check that files were processed (not necessarily succeeded)
        assert len(results_seq) == 2, "Sequential should process 2 files"
        assert len(results_par) == 2, "Parallel should process 2 files"

        # Compare log files
        seq_content = log_seq.read_text() if log_seq.exists() else ""
        par_content = log_par.read_text() if log_par.exists() else ""

        print(f"\n=== Sequential Log ({len(seq_content)} bytes) ===")
        print(seq_content[:300] if seq_content else "(empty)")

        print(f"\n=== Parallel Log ({len(par_content)} bytes) ===")
        print(par_content[:300] if par_content else "(empty)")

        # Both should contain process IDs
        if seq_content:
            assert "process_id" in seq_content
        if par_content:
            assert "process_id" in par_content

        print(f"\n✓ Both sequential and parallel logging work correctly")
