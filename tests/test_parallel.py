"""
Tests for parallel file validation functionality.

Tests cover:
- Parallel processing with validate_reads()
- Parallel processing with validate_genomes()
- Parallel processing with validate_features_list()
- Thread-safe logging
- Error handling in parallel mode
- Sequential vs parallel execution
"""

import pytest
from pathlib import Path
import tempfile
import shutil
from validation_pkg import (
    ConfigManager,
    validate_reads,
    validate_genomes,
    validate_features_list,
    ReadValidator,
    GenomeValidator,
    FeatureValidator,
    get_logger,
    setup_logging
)


@pytest.fixture
def temp_output_dir():
    """Create a temporary output directory."""
    temp_dir = tempfile.mkdtemp()
    yield Path(temp_dir)
    # Cleanup after test
    shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture
def sample_fastq_files(tmp_path):
    """Create multiple sample FASTQ files for testing."""
    files = []
    for i in range(4):
        fastq_file = tmp_path / f"sample_{i}.fastq"
        with open(fastq_file, 'w') as f:
            # Write 10 reads per file
            for j in range(10):
                f.write(f"@read_{i}_{j}\n")
                f.write("ATCGATCGATCG\n")
                f.write("+\n")
                f.write("IIIIIIIIIIII\n")
        files.append(fastq_file)
    return files


@pytest.fixture
def sample_fasta_files(tmp_path):
    """Create multiple sample FASTA files for testing."""
    files = []
    for i in range(3):
        fasta_file = tmp_path / f"genome_{i}.fasta"
        with open(fasta_file, 'w') as f:
            f.write(f">sequence_{i}\n")
            f.write("ATCGATCGATCGATCGATCG\n" * 50)  # 1000 bp
        files.append(fasta_file)
    return files


@pytest.fixture
def sample_gff_files(tmp_path):
    """Create multiple sample GFF files for testing."""
    files = []
    for i in range(3):
        gff_file = tmp_path / f"features_{i}.gff"
        with open(gff_file, 'w') as f:
            f.write("##gff-version 3\n")
            for j in range(5):
                start = j * 100 + 1
                end = start + 99
                f.write(f"chr1\t.\tgene\t{start}\t{end}\t.\t+\t.\tID=gene_{i}_{j}\n")
        files.append(gff_file)
    return files


class TestParallelReads:
    """Tests for parallel read validation."""

    def test_parallel_reads_basic(self, sample_fastq_files, temp_output_dir):
        """Test basic parallel processing of read files."""
        # Setup logger
        logger = get_logger()
        logger.clear_issues()

        # Create mock ReadConfig objects
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import ReadFormat, CodingType

        read_configs = []
        for fastq_file in sample_fastq_files:
            read_config = ReadConfig(
                filename=fastq_file.name,
                filepath=fastq_file,
                detected_format=ReadFormat.FASTQ,
                coding_type=CodingType.NONE,
                ngs_type='illumina',
                settings_dict={},
                extra={}
            )
            read_configs.append(read_config)

        # Test with max_workers=2
        settings = ReadValidator.Settings()
        settings = settings.update(max_workers=2)

        results = validate_reads(read_configs, temp_output_dir, settings)

        # Assertions
        assert len(results) == 4
        assert all(r['success'] for r in results)
        assert all((temp_output_dir / 'illumina' / r['filename']).exists() for r in results)

    def test_parallel_reads_with_errors(self, sample_fastq_files, temp_output_dir):
        """Test parallel processing with validation errors."""
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import ReadFormat, CodingType

        # Create one invalid file
        invalid_file = sample_fastq_files[0].parent / "invalid.fastq"
        with open(invalid_file, 'w') as f:
            f.write("Not a valid FASTQ file\n")

        read_configs = []
        # Add valid files
        for fastq_file in sample_fastq_files[1:]:
            read_config = ReadConfig(
                filename=fastq_file.name,
                filepath=fastq_file,
                detected_format=ReadFormat.FASTQ,
                coding_type=CodingType.NONE,
                ngs_type='illumina',
                settings_dict={},
                extra={}
            )
            read_configs.append(read_config)

        # Add invalid file
        invalid_config = ReadConfig(
            filename=invalid_file.name,
            filepath=invalid_file,
            detected_format=ReadFormat.FASTQ,
            coding_type=CodingType.NONE,
            ngs_type='illumina',
            settings_dict={},
            extra={}
        )
        read_configs.append(invalid_config)

        settings = ReadValidator.Settings()
        settings = settings.update(max_workers=2)

        results = validate_reads(read_configs, temp_output_dir, settings)

        # Check that we got results for all files
        assert len(results) == 4

        # Check that 3 succeeded and 1 failed
        successes = [r for r in results if r['success']]
        failures = [r for r in results if not r['success']]
        assert len(successes) == 3
        assert len(failures) == 1
        assert failures[0]['filename'] == invalid_file.name

    def test_sequential_vs_parallel(self, sample_fastq_files, temp_output_dir):
        """Test that sequential and parallel produce same results."""
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import ReadFormat, CodingType

        read_configs = []
        for fastq_file in sample_fastq_files:
            read_config = ReadConfig(
                filename=fastq_file.name,
                filepath=fastq_file,
                detected_format=ReadFormat.FASTQ,
                coding_type=CodingType.NONE,
                ngs_type='illumina',
                settings_dict={},
                extra={}
            )
            read_configs.append(read_config)

        # Sequential execution
        temp_sequential = temp_output_dir / "sequential"
        temp_sequential.mkdir(parents=True, exist_ok=True)

        settings_seq = ReadValidator.Settings()
        results_seq = validate_reads(read_configs, temp_sequential, settings_seq)

        # Parallel execution
        temp_parallel = temp_output_dir / "parallel"
        temp_parallel.mkdir(parents=True, exist_ok=True)

        settings_parallel = ReadValidator.Settings()
        settings_parallel = settings_parallel.update(max_workers=2)
        results_parallel = validate_reads(read_configs, temp_parallel, settings_parallel)

        # Both should succeed
        assert all(r['success'] for r in results_seq)
        assert all(r['success'] for r in results_parallel)

        # Same number of output files
        seq_files = list((temp_sequential / 'illumina').glob('*.fastq'))
        parallel_files = list((temp_parallel / 'illumina').glob('*.fastq'))
        assert len(seq_files) == len(parallel_files) == 4


class TestParallelGenomes:
    """Tests for parallel genome validation."""

    def test_parallel_genomes_basic(self, sample_fasta_files, temp_output_dir):
        """Test basic parallel processing of genome files."""
        from validation_pkg.config_manager import GenomeConfig
        from validation_pkg.utils.formats import GenomeFormat, CodingType

        genome_configs = []
        for fasta_file in sample_fasta_files:
            genome_config = GenomeConfig(
                filename=fasta_file.name,
                filepath=fasta_file,
                detected_format=GenomeFormat.FASTA,
                coding_type=CodingType.NONE,
                settings_dict={},
                extra={}
            )
            genome_configs.append(genome_config)

        settings = GenomeValidator.Settings()
        settings = settings.update(max_workers=2)

        results = validate_genomes(genome_configs, temp_output_dir, settings)

        # Assertions
        assert len(results) == 3
        assert all(r['success'] for r in results)

    def test_parallel_genomes_single_file(self, sample_fasta_files, temp_output_dir):
        """Test that single file falls back to sequential."""
        from validation_pkg.config_manager import GenomeConfig
        from validation_pkg.utils.formats import GenomeFormat, CodingType

        # Only one file
        genome_config = GenomeConfig(
            filename=sample_fasta_files[0].name,
            filepath=sample_fasta_files[0],
            detected_format=GenomeFormat.FASTA,
            coding_type=CodingType.NONE,
            settings_dict={},
            extra={}
        )

        settings = GenomeValidator.Settings()
        settings = settings.update(max_workers=4)  # Request parallel but only 1 file

        results = validate_genomes([genome_config], temp_output_dir, settings)

        assert len(results) == 1
        assert results[0]['success']


class TestParallelFeatures:
    """Tests for parallel feature validation."""

    def test_parallel_features_basic(self, sample_gff_files, temp_output_dir):
        """Test basic parallel processing of feature files."""
        from validation_pkg.config_manager import FeatureConfig
        from validation_pkg.utils.formats import FeatureFormat, CodingType

        feature_configs = []
        for gff_file in sample_gff_files:
            feature_config = FeatureConfig(
                filename=gff_file.name,
                filepath=gff_file,
                detected_format=FeatureFormat.GFF,
                coding_type=CodingType.NONE,
                settings_dict={},
                extra={}
            )
            feature_configs.append(feature_config)

        settings = FeatureValidator.Settings()
        settings = settings.update(max_workers=2)

        results = validate_features_list(feature_configs, temp_output_dir, settings)

        # Assertions
        assert len(results) == 3
        assert all(r['success'] for r in results)


class TestThreadSafeLogging:
    """Tests for thread-safe logging in parallel mode."""

    def test_logging_thread_safety(self, sample_fastq_files, temp_output_dir):
        """Test that logging works correctly in parallel mode."""
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import ReadFormat, CodingType

        # Setup logger
        logger = get_logger()
        logger.clear_issues()

        read_configs = []
        for fastq_file in sample_fastq_files:
            read_config = ReadConfig(
                filename=fastq_file.name,
                filepath=fastq_file,
                detected_format=ReadFormat.FASTQ,
                coding_type=CodingType.NONE,
                ngs_type='illumina',
                settings_dict={},
                extra={}
            )
            read_configs.append(read_config)

        settings = ReadValidator.Settings()
        settings = settings.update(max_workers=4)

        results = validate_reads(read_configs, temp_output_dir, settings)

        # All should succeed
        assert all(r['success'] for r in results)

        # Logger should not crash or lose data
        summary = logger.get_summary()
        assert isinstance(summary, dict)
        assert 'total_issues' in summary


class TestMaxWorkersValidation:
    """Tests for max_workers parameter validation."""

    def test_max_workers_none(self, sample_fastq_files, temp_output_dir):
        """Test that None max_workers uses sequential processing."""
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import ReadFormat, CodingType

        read_configs = []
        for fastq_file in sample_fastq_files:
            read_config = ReadConfig(
                filename=fastq_file.name,
                filepath=fastq_file,
                detected_format=ReadFormat.FASTQ,
                coding_type=CodingType.NONE,
                ngs_type='illumina',
                settings_dict={},
                extra={}
            )
            read_configs.append(read_config)

        settings = ReadValidator.Settings()
        # max_workers is None by default
        assert settings.max_workers is None

        results = validate_reads(read_configs, temp_output_dir, settings)
        assert all(r['success'] for r in results)

    def test_max_workers_one(self, sample_fastq_files, temp_output_dir):
        """Test that max_workers=1 uses sequential processing."""
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import ReadFormat, CodingType

        read_configs = []
        for fastq_file in sample_fastq_files[:2]:  # Just 2 files
            read_config = ReadConfig(
                filename=fastq_file.name,
                filepath=fastq_file,
                detected_format=ReadFormat.FASTQ,
                coding_type=CodingType.NONE,
                ngs_type='illumina',
                settings_dict={},
                extra={}
            )
            read_configs.append(read_config)

        settings = ReadValidator.Settings()
        settings = settings.update(max_workers=1)

        results = validate_reads(read_configs, temp_output_dir, settings)
        assert all(r['success'] for r in results)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
