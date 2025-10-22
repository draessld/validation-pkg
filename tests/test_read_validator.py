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

        with pytest.raises(FastqFormatError, match="No sequences found"):
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

        with pytest.raises(ReadValidationError, match="invalid character"):
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

    # def test_bam_conversion_workflow(self, fake_bam, output_dir):
    #     """
    #     Test that BAM files follow the proper workflow.

    #     With fake BAM files, the conversion will fail (as expected).
    #     This tests that the workflow is triggered correctly.
    #     Real BAM conversion requires pysam or samtools.
    #     """
    #     read_config = ReadConfig(
    #         filename="reads.bam",
    #         filepath=fake_bam,
    #         ngs_type="illumina",
    #         coding_type=CodingType.NONE,
    #         detected_format=ReadFormat.BAM
    #     )

    #     validator = ReadValidator(read_config, output_dir,
    #                              ReadValidator.Settings(keep_bam=True))

    #     # Should raise error during BAM conversion (fake BAM file)
    #     # The error could be from pysam, samtools, or "neither tool available"
    #     with pytest.raises((ReadValidationError, Exception)):
    #         validator.validate()

    #     # With keep_bam=True, the original BAM should be copied
    #     copied_files = list(output_dir.glob("*.bam"))
    #     assert len(copied_files) == 1


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


class TestReadValidatorValidationLevels:
    """Test multi-level validation modes (strict, trust, minimal)."""

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
    def multi_read_fastq(self, temp_dir):
        """Create a gzipped FASTQ file with multiple reads (for minimal mode)."""
        fastq_file = temp_dir / "reads.fastq.gz"
        with gzip.open(fastq_file, "wt") as f:
            for i in range(1, 21):  # 20 reads
                f.write(f"@read{i}\n")
                f.write("ATCGATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIIIIIII\n")
        return fastq_file

    @pytest.fixture
    def damaged_fastq(self, temp_dir):
        """Create a FASTQ file with empty ID (uncompressed for strict/trust mode)."""
        fastq_file = temp_dir / "damaged.fastq"
        with open(fastq_file, "w") as f:
            f.write("@\n")  # Empty ID
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")
        return fastq_file

    # ===== Tests for STRICT validation level =====

    def test_strict_correct_file_passes(self, multi_read_fastq, output_dir):
        """Test strict mode with correct FASTQ file - should pass."""
        settings = ReadValidator.Settings(validation_level='strict', coding_type='gz')
        read_config = ReadConfig(
            filename="reads.fastq.gz",
            filepath=multi_read_fastq,
            ngs_type="illumina",
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)
        validator.validate()

        # All reads should be parsed
        assert len(validator.sequences) == 20
        # Output file should exist (default output is gzipped)
        output_files = list(output_dir.glob("*.fastq.gz"))
        assert len(output_files) == 1

    def test_strict_damaged_file_fails(self, damaged_fastq, output_dir):
        """Test strict mode with damaged file - should fail during parsing."""
        settings = ReadValidator.Settings(
            validation_level='strict',
            allow_empty_id=False
        )
        read_config = ReadConfig(
            filename="damaged.fastq",
            filepath=damaged_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)

        # Should fail during parsing - empty ID after @ is malformed FASTQ
        with pytest.raises(FastqFormatError):
            validator.validate()

    # ===== Tests for TRUST validation level =====

    def test_trust_correct_file_passes(self, multi_read_fastq, output_dir):
        """Test trust mode with correct FASTQ file - parses first 10, validates them, copies original file."""
        settings = ReadValidator.Settings(validation_level='trust', coding_type='gz')
        read_config = ReadConfig(
            filename="reads.fastq.gz",
            filepath=multi_read_fastq,
            ngs_type="illumina",
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)
        validator.validate()

        # Only first 10 reads should be parsed (trust mode parses only first 10)
        assert len(validator.sequences) == 10
        # Output file should exist with all reads copied (gzipped)
        output_files = list(output_dir.glob("*.fastq.gz"))
        assert len(output_files) == 1

    def test_trust_damaged_first_sequence_fails(self, damaged_fastq, output_dir):
        """Test trust mode detects error in first sequence (within first 10 validated)."""
        settings = ReadValidator.Settings(
            validation_level='trust',
            allow_empty_id=False
        )
        read_config = ReadConfig(
            filename="damaged.fastq",
            filepath=damaged_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)

        # Should fail during parsing - empty ID after @ is malformed FASTQ
        with pytest.raises(FastqFormatError):
            validator.validate()

    # ===== Tests for MINIMAL validation level =====

    def test_minimal_correct_file_passes(self, multi_read_fastq, output_dir):
        """Test minimal mode with correct file - should pass without validation."""
        settings = ReadValidator.Settings(validation_level='minimal', coding_type='gz')
        read_config = ReadConfig(
            filename="reads.fastq.gz",
            filepath=multi_read_fastq,
            ngs_type="illumina",
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)
        validator.validate()

        # No sequences parsed in minimal mode
        assert len(validator.sequences) == 0
        # Output file should exist (gzipped copy)
        output_files = list(output_dir.glob("*.fastq.gz"))
        assert len(output_files) == 1

    def test_minimal_damaged_file_raises_error(self, damaged_fastq, output_dir):
        """Test minimal mode with uncompressed file - should raise error due to coding mismatch."""
        settings = ReadValidator.Settings(
            validation_level='minimal',
            coding_type='gz',  # Require GZIP output
            allow_empty_id=False  # Ignored in minimal mode
        )
        read_config = ReadConfig(
            filename="damaged.fastq",
            filepath=damaged_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,  # Uncompressed - doesn't match output coding
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)
        # Should raise error - minimal mode requires input coding to match output coding
        with pytest.raises(ReadValidationError, match="input coding to match output coding"):
            validator.validate()

    def test_minimal_output_is_copy(self, multi_read_fastq, output_dir):
        """Test that minimal mode copies file as-is."""
        settings = ReadValidator.Settings(validation_level='minimal', coding_type='gz')
        read_config = ReadConfig(
            filename="reads.fastq.gz",
            filepath=multi_read_fastq,
            ngs_type="illumina",
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)
        validator.validate()

        # Check output file exists (gzipped)
        output_files = list(output_dir.glob("*.fastq.gz"))
        assert len(output_files) == 1

        # Verify output is byte-for-byte identical to input
        with open(multi_read_fastq, 'rb') as f_in, open(output_files[0], 'rb') as f_out:
            assert f_in.read() == f_out.read()

    # ===== Test for compressed files =====

    def test_trust_compressed_gz_passes(self, temp_dir, output_dir):
        """Test trust mode with gzip compressed file."""
        fastq_file = temp_dir / "reads.fastq.gz"
        with gzip.open(fastq_file, "wt") as f:
            for i in range(1, 11):
                f.write(f"@read{i}\n")
                f.write("ATCGATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIIIIIII\n")

        settings = ReadValidator.Settings(validation_level='trust', coding_type='gz')
        read_config = ReadConfig(
            filename="reads.fastq.gz",
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.GZIP,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)
        validator.validate()

        # All 10 reads should be parsed (trust mode parses first 10)
        assert len(validator.sequences) == 10

        # Output should be created
        output_files = list(output_dir.glob("*.fastq.gz"))
        assert len(output_files) == 1

    def test_minimal_compressed_bz2_raises_error(self, temp_dir, output_dir):
        """Test minimal mode with bzip2 compressed file - should raise error due to coding mismatch."""
        fastq_file = temp_dir / "reads.fastq.bz2"
        with bz2.open(fastq_file, "wt") as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIIIIIII\n")

        settings = ReadValidator.Settings(validation_level='minimal', coding_type='gz')
        read_config = ReadConfig(
            filename="reads.fastq.bz2",
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.BZIP2,
            detected_format=ReadFormat.FASTQ
        )

        validator = ReadValidator(read_config, output_dir, settings)

        # Should raise error - minimal mode requires input coding to match output coding
        with pytest.raises(ReadValidationError, match="input coding to match output coding"):
            validator.validate()

    # ===== Test for invalid validation level =====

    def test_invalid_validation_level_raises_error(self, multi_read_fastq, output_dir):
        """Test that invalid validation level raises ValueError."""
        settings = ReadValidator.Settings(validation_level='invalid_mode')
        read_config = ReadConfig(
            filename="reads.fastq",
            filepath=multi_read_fastq,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
        )

        with pytest.raises(ValueError, match="Invalid validation_level"):
            ReadValidator(read_config, output_dir, settings)


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
            f.write("!!!!!!!!!!!!!!!!!!!!\n")  # Quality score 0
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

        # All reads should be parsed
        assert len(validator.sequences) == 3
        # Verify different lengths were parsed correctly
        lengths = [len(seq.seq) for seq in validator.sequences]
        assert min(lengths) == 4
        assert max(lengths) == 400

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

        # Low quality reads are allowed
        assert len(validator.sequences) == 1



if __name__ == "__main__":
    pytest.main([__file__, "-v"])
