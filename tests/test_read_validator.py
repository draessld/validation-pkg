"""
Comprehensive tests for ReadValidator.

Tests cover:
- File format detection and parsing (FASTQ, BAM)
- Compression handling (gzip, bzip2, uncompressed)
- Validation rules (empty IDs, duplicate IDs, invalid characters)
- BAM file handling (pass-through, keep_bam, conversion)
- Statistics collection (read counts, length, GC content, quality scores)
- Output generation with various compression formats
- Error handling and edge cases
"""

import pytest
import tempfile
import gzip
import bz2
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from validation_pkg.config_manager import ReadConfig
from validation_pkg.validators.read_validator import ReadValidator
from validation_pkg.exceptions import (
    ReadValidationError,
    FastqFormatError,
    BamFormatError,
    CompressionError,
)
from validation_pkg.utils.formats import ReadFormat, CodingType


class TestReadValidatorInitialization:
    """Test ReadValidator initialization."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        """Create output directory."""
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def simple_fastq(self, temp_dir):
        """Create a simple FASTQ file."""
        fastq_file = temp_dir / "reads.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIII\n")
        return fastq_file

    def test_init_with_defaults(self, simple_fastq, output_dir):
        """Test initialization with default settings."""
        read_config = ReadConfig(
            filename="reads.fastq",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir)

        assert validator.input_path == simple_fastq
        assert validator.output_dir == output_dir
        assert validator.settings is not None
        assert validator.sequences == []
        assert isinstance(validator.statistics, dict)

    def test_init_with_custom_settings(self, simple_fastq, output_dir):
        """Test initialization with custom settings."""
        settings = ReadValidator.Settings(
            check_invalid_chars=True,
            allow_duplicate_ids=False,
            keep_bam=False
        )

        read_config = ReadConfig(
            filename="reads.fastq",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)

        assert validator.settings.check_invalid_chars is True
        assert validator.settings.allow_duplicate_ids is False
        assert validator.settings.keep_bam is False


class TestReadValidatorParsing:
    """Test file parsing functionality."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def simple_fastq(self, temp_dir):
        """Create a simple FASTQ file with two reads."""
        fastq_file = temp_dir / "reads.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")
            f.write("@read2\n")
            f.write("GCTAGCTAGCTAGCTAGCTA\n")
            f.write("+\n")
            f.write("HHHHHHHHHHHHHHHHHHHH\n")
        return fastq_file

    @pytest.fixture
    def empty_fastq(self, temp_dir):
        """Create an empty FASTQ file."""
        fastq_file = temp_dir / "empty.fastq"
        with open(fastq_file, "w") as f:
            f.write("")
        return fastq_file

    @pytest.fixture
    def invalid_fastq(self, temp_dir):
        """Create an invalid FASTQ file (missing quality scores)."""
        fastq_file = temp_dir / "invalid.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGATCG\n")
            # Missing + and quality lines
        return fastq_file

    def test_parse_simple_fastq(self, simple_fastq, output_dir):
        """Test parsing a simple FASTQ file."""
        read_config = ReadConfig(
            filename="reads.fastq",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir)
        validator.validate()

        assert len(validator.sequences) == 2
        assert validator.sequences[0].id == "read1"
        assert validator.sequences[1].id == "read2"
        assert str(validator.sequences[0].seq) == "ATCGATCGATCGATCGATCG"

    def test_parse_empty_file_raises_error(self, empty_fastq, output_dir):
        """Test that empty FASTQ file raises ReadValidationError."""
        read_config = ReadConfig(
            filename="empty.fastq",
            filepath=empty_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir)

        with pytest.raises(ReadValidationError, match="No sequences found"):
            validator.validate()

    def test_parse_invalid_fastq_raises_error(self, invalid_fastq, output_dir):
        """Test that invalid FASTQ raises FastqFormatError."""
        read_config = ReadConfig(
            filename="invalid.fastq",
            filepath=invalid_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir)

        with pytest.raises((FastqFormatError, ReadValidationError)):
            validator.validate()


class TestReadValidatorCompression:
    """Test compression handling."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def compressed_fastq_gz(self, temp_dir):
        """Create a gzip compressed FASTQ file."""
        fastq_file = temp_dir / "reads.fastq.gz"
        with gzip.open(fastq_file, "wt") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")
        return fastq_file

    @pytest.fixture
    def compressed_fastq_bz2(self, temp_dir):
        """Create a bzip2 compressed FASTQ file."""
        fastq_file = temp_dir / "reads.fastq.bz2"
        with bz2.open(fastq_file, "wt") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")
        return fastq_file

    def test_parse_gzip_compressed(self, compressed_fastq_gz, output_dir):
        """Test parsing gzip compressed FASTQ."""
        read_config = ReadConfig(
            filename="reads.fastq.gz",
            filepath=compressed_fastq_gz,
            ngs_type="illumina",
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir)
        validator.validate()

        assert len(validator.sequences) == 1
        assert str(validator.sequences[0].seq) == "ATCGATCGATCGATCGATCG"

    def test_parse_bzip2_compressed(self, compressed_fastq_bz2, output_dir):
        """Test parsing bzip2 compressed FASTQ."""
        read_config = ReadConfig(
            filename="reads.fastq.bz2",
            filepath=compressed_fastq_bz2,
            ngs_type="illumina",
            coding_type=CodingType.BZIP2,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir)
        validator.validate()

        assert len(validator.sequences) == 1
        assert str(validator.sequences[0].seq) == "ATCGATCGATCGATCGATCG"


class TestReadValidatorValidation:
    """Test validation rules."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def fastq_with_duplicates(self, temp_dir):
        """Create FASTQ with duplicate read IDs."""
        fastq_file = temp_dir / "duplicates.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")
            f.write("@read1\n")
            f.write("GCTAGCTAGCTAGCTAGCTA\n")
            f.write("+\n")
            f.write("HHHHHHHHHHHHHHHHHHHH\n")
        return fastq_file

    @pytest.fixture
    def fastq_with_invalid_chars(self, temp_dir):
        """Create FASTQ with invalid characters."""
        fastq_file = temp_dir / "invalid_chars.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGXYZATCG\n")  # X, Y, Z are invalid
            f.write("+\n")
            f.write("IIIIIIIIIII\n")
        return fastq_file

    def test_duplicate_ids_allowed_by_default(self, fastq_with_duplicates, output_dir):
        """Test that duplicate IDs are allowed by default."""
        read_config = ReadConfig(
            filename="duplicates.fastq",
            filepath=fastq_with_duplicates,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        # Default settings allow duplicates
        validator = ReadValidator(read_config, output_dir)
        validator.validate()  # Should not raise

        assert len(validator.sequences) == 2

    def test_duplicate_ids_rejected_when_disabled(self, fastq_with_duplicates, output_dir):
        """Test that duplicate IDs are rejected when allow_duplicate_ids=False."""
        settings = ReadValidator.Settings(allow_duplicate_ids=False)

        read_config = ReadConfig(
            filename="duplicates.fastq",
            filepath=fastq_with_duplicates,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)

        with pytest.raises(ReadValidationError, match="Duplicate sequence IDs"):
            validator.validate()

    def test_invalid_chars_detected(self, fastq_with_invalid_chars, output_dir):
        """Test that invalid characters are detected when check_invalid_chars=True."""
        settings = ReadValidator.Settings(check_invalid_chars=True)

        read_config = ReadConfig(
            filename="invalid_chars.fastq",
            filepath=fastq_with_invalid_chars,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)

        with pytest.raises(ReadValidationError, match="invalid chars"):
            validator.validate()

    def test_invalid_chars_ignored_by_default(self, fastq_with_invalid_chars, output_dir):
        """Test that invalid characters are ignored by default."""
        read_config = ReadConfig(
            filename="invalid_chars.fastq",
            filepath=fastq_with_invalid_chars,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        # Default settings don't check invalid chars
        validator = ReadValidator(read_config, output_dir)
        validator.validate()  # Should not raise


class TestReadValidatorBAMHandling:
    """Test BAM file handling."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def fake_bam(self, temp_dir):
        """Create a fake BAM file (for testing pass-through)."""
        bam_file = temp_dir / "reads.bam"
        # Create a dummy file - actual BAM parsing tested separately
        with open(bam_file, "wb") as f:
            f.write(b"BAM\x01")  # BAM magic number
        return bam_file

    def test_bam_pass_through(self, fake_bam, output_dir):
        """Test that BAM files are passed through without processing."""
        read_config = ReadConfig(
            filename="reads.bam",
            filepath=fake_bam,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.BAM
        )

        validator = ReadValidator(read_config, output_dir)

        # Should raise ReadValidationError about BAM not being fully implemented
        with pytest.raises(ReadValidationError, match="BAM file processing not fully implemented"):
            validator.validate()

        # But the file should be copied to output
        copied_files = list(output_dir.glob("*.bam"))
        assert len(copied_files) == 1


class TestReadValidatorStatistics:
    """Test statistics collection."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def fastq_for_stats(self, temp_dir):
        """Create FASTQ with known statistics."""
        fastq_file = temp_dir / "stats.fastq"
        with open(fastq_file, "w") as f:
            # Read 1: 8bp, 50% GC, quality 40 (I = 73 - 33 = 40)
            f.write("@read1\n")
            f.write("ATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIII\n")
            # Read 2: 12bp, 100% GC, quality 39 (H = 72 - 33 = 39)
            f.write("@read2\n")
            f.write("GCGCGCGCGCGC\n")
            f.write("+\n")
            f.write("HHHHHHHHHHHH\n")
        return fastq_file

    def test_statistics_collection(self, fastq_for_stats, output_dir):
        """Test that statistics are collected correctly."""
        read_config = ReadConfig(
            filename="stats.fastq",
            filepath=fastq_for_stats,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir)
        validator.validate()

        stats = validator.get_statistics()

        assert stats['num_sequences'] == 2
        assert stats['total_length'] == 20
        assert stats['min_length'] == 8
        assert stats['max_length'] == 12
        assert stats['avg_length'] == 10.0
        # GC content: (4 + 12) / 20 = 80%
        assert 79 < stats['gc_content'] < 81

        # Quality statistics should be present for FASTQ
        if 'mean_quality' in stats:
            assert 38 < stats['mean_quality'] < 41


class TestReadValidatorOutput:
    """Test output generation."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def simple_fastq(self, temp_dir):
        """Create a simple FASTQ file."""
        fastq_file = temp_dir / "reads.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")
        return fastq_file

    def test_output_uncompressed(self, simple_fastq, output_dir):
        """Test generating uncompressed output."""
        settings = ReadValidator.Settings(coding_type=None)

        read_config = ReadConfig(
            filename="reads.fastq",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)
        validator.validate()

        # Check output file exists
        output_files = list(output_dir.glob("*.fastq"))
        assert len(output_files) == 1
        assert output_files[0].suffix == ".fastq"

    def test_output_gzip_compressed(self, simple_fastq, output_dir):
        """Test generating gzip compressed output."""
        settings = ReadValidator.Settings(coding_type='gz')

        read_config = ReadConfig(
            filename="reads.fastq",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)
        validator.validate()

        # Check gzip output file exists
        output_files = list(output_dir.glob("*.fastq.gz"))
        assert len(output_files) == 1

        # Verify it's actually gzipped and contains correct data
        with gzip.open(output_files[0], 'rt') as f:
            records = list(SeqIO.parse(f, 'fastq'))
            assert len(records) == 1

    def test_output_bzip2_compressed(self, simple_fastq, output_dir):
        """Test generating bzip2 compressed output."""
        settings = ReadValidator.Settings(coding_type='bz2')

        read_config = ReadConfig(
            filename="reads.fastq",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)
        validator.validate()

        # Check bzip2 output file exists
        output_files = list(output_dir.glob("*.fastq.bz2"))
        assert len(output_files) == 1

        # Verify it's actually bzip2'd and contains correct data
        with bz2.open(output_files[0], 'rt') as f:
            records = list(SeqIO.parse(f, 'fastq'))
            assert len(records) == 1

    def test_output_with_suffix(self, simple_fastq, output_dir):
        """Test output filename with custom suffix."""
        settings = ReadValidator.Settings(output_filename_suffix="filtered")

        read_config = ReadConfig(
            filename="reads.fastq",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)
        validator.validate()

        output_files = list(output_dir.glob("*_filtered.fastq"))
        assert len(output_files) == 1

    def test_output_with_subdirectory(self, simple_fastq, output_dir):
        """Test output to custom subdirectory."""
        settings = ReadValidator.Settings(output_subdir_name="reads")

        read_config = ReadConfig(
            filename="reads.fastq",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)
        validator.validate()

        subdir = output_dir / "reads"
        assert subdir.exists()
        output_files = list(subdir.glob("*.fastq"))
        assert len(output_files) == 1


class TestReadValidatorNGSTypes:
    """Test different NGS types."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def simple_fastq(self, temp_dir):
        """Create a simple FASTQ file."""
        fastq_file = temp_dir / "reads.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")
        return fastq_file

    def test_illumina_reads(self, simple_fastq, output_dir):
        """Test processing Illumina reads."""
        read_config = ReadConfig(
            filename="reads.fastq",
            filepath=simple_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir)
        validator.validate()

        assert len(validator.sequences) == 1

    def test_ont_reads(self, simple_fastq, output_dir):
        """Test processing Oxford Nanopore reads."""
        read_config = ReadConfig(
            filename="reads.fastq",
            filepath=simple_fastq,
            ngs_type="ont",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir)
        validator.validate()

        assert len(validator.sequences) == 1

    def test_pacbio_reads(self, simple_fastq, output_dir):
        """Test processing PacBio reads."""
        read_config = ReadConfig(
            filename="reads.fastq",
            filepath=simple_fastq,
            ngs_type="pacbio",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir)
        validator.validate()

        assert len(validator.sequences) == 1


class TestReadValidatorEdgeCases:
    """Test edge cases and error handling."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        out_dir = temp_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def fastq_with_varying_lengths(self, temp_dir):
        """Create FASTQ with reads of varying lengths."""
        fastq_file = temp_dir / "varying.fastq"
        with open(fastq_file, "w") as f:
            f.write("@short_read\n")
            f.write("ATCG\n")
            f.write("+\n")
            f.write("IIII\n")
            f.write("@medium_read\n")
            f.write("ATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIII\n")
            f.write("@long_read\n")
            f.write("ATCGATCG" * 50 + "\n")  # 400bp
            f.write("+\n")
            f.write("I" * 400 + "\n")
        return fastq_file

    @pytest.fixture
    def fastq_with_low_quality(self, temp_dir):
        """Create FASTQ with low quality scores."""
        fastq_file = temp_dir / "low_quality.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("!!!!!!!!!!!!!!!!!!!\n")  # Quality score 0
        return fastq_file

    def test_varying_read_lengths(self, fastq_with_varying_lengths, output_dir):
        """Test handling reads of varying lengths."""
        read_config = ReadConfig(
            filename="varying.fastq",
            filepath=fastq_with_varying_lengths,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir)
        validator.validate()

        assert len(validator.sequences) == 3
        stats = validator.get_statistics()
        assert stats['min_length'] == 4
        assert stats['max_length'] == 400

    def test_low_quality_reads(self, fastq_with_low_quality, output_dir):
        """Test handling reads with low quality scores."""
        read_config = ReadConfig(
            filename="low_quality.fastq",
            filepath=fastq_with_low_quality,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir)
        validator.validate()

        assert len(validator.sequences) == 1
        # Low quality is allowed, just reported in stats
        stats = validator.get_statistics()
        if 'min_quality' in stats:
            assert stats['min_quality'] == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
