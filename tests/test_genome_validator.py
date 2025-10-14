"""
Comprehensive tests for GenomeValidator.

Tests cover:
- File format detection and parsing (FASTA, GenBank)
- Compression handling (gzip, bzip2, uncompressed)
- Validation rules (empty sequences, duplicate IDs, invalid characters)
- Editing specifications (min length, sequence prefix, plasmid split)
- Statistics collection
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

from validation_pkg.config_manager import GenomeConfig
from validation_pkg.validators.genome_validator import GenomeValidator
from validation_pkg.exceptions import (
    GenomeValidationError,
    FastaFormatError,
    GenBankFormatError,
    CompressionError,
)
from validation_pkg.utils.formats import GenomeFormat, CodingType


class TestGenomeValidatorInitialization:
    """Test GenomeValidator initialization."""

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
    def simple_fasta(self, temp_dir):
        """Create a simple FASTA file."""
        fasta_file = temp_dir / "genome.fasta"
        with open(fasta_file, "w") as f:
            f.write(">seq1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
        return fasta_file

    def test_init_with_defaults(self, simple_fasta, output_dir):
        """Test initialization with default settings."""
        genome_config = GenomeConfig(
            filename="genome.fasta",
            filepath=simple_fasta,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        # Use min_sequence_length=0 to avoid filtering test sequences
        settings = GenomeValidator.Settings(min_sequence_length=0)
        validator = GenomeValidator(genome_config, output_dir, settings)

        assert validator.input_path == simple_fasta
        assert validator.output_dir == output_dir
        assert validator.settings is not None
        assert validator.sequences == []
        assert isinstance(validator.statistics, dict)

    def test_init_with_custom_settings(self, simple_fasta, output_dir):
        """Test initialization with custom settings."""
        settings = GenomeValidator.Settings(
            min_sequence_length=500,
            sequence_prefix="chr",
            allow_duplicate_ids=False
        )

        genome_config = GenomeConfig(
            filename="genome.fasta",
            filepath=simple_fasta,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, settings)

        assert validator.settings.min_sequence_length == 500
        assert validator.settings.sequence_prefix == "chr"
        assert validator.settings.allow_duplicate_ids is False


class TestGenomeValidatorParsing:
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
    def default_settings(self):
        """Default settings that don't filter sequences."""
        return GenomeValidator.Settings(min_sequence_length=0)

    @pytest.fixture
    def simple_fasta(self, temp_dir):
        """Create a simple FASTA file with two sequences."""
        fasta_file = temp_dir / "genome.fasta"
        with open(fasta_file, "w") as f:
            f.write(">seq1 description1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write(">seq2 description2\n")
            f.write("GCTAGCTAGCTAGCTAGCTA\n")
        return fasta_file

    @pytest.fixture
    def simple_genbank(self, temp_dir):
        """Create a simple GenBank file."""
        gb_file = temp_dir / "genome.gbk"
        record = SeqRecord(
            Seq("ATCGATCGATCGATCGATCG"),
            id="TEST001",
            name="test_sequence",
            description="Test sequence",
        )
        # GenBank requires molecule_type annotation
        record.annotations["molecule_type"] = "DNA"
        SeqIO.write([record], gb_file, "genbank")
        return gb_file

    @pytest.fixture
    def empty_fasta(self, temp_dir):
        """Create an empty FASTA file."""
        fasta_file = temp_dir / "empty.fasta"
        with open(fasta_file, "w") as f:
            f.write("")
        return fasta_file

    @pytest.fixture
    def invalid_fasta(self, temp_dir):
        """Create an invalid FASTA file (no headers)."""
        fasta_file = temp_dir / "invalid.fasta"
        with open(fasta_file, "w") as f:
            f.write("This is not a valid FASTA file\n")
            f.write("No headers at all\n")
        return fasta_file

    def test_parse_simple_fasta(self, simple_fasta, output_dir, default_settings):
        """Test parsing a simple FASTA file."""
        genome_config = GenomeConfig(
            filename="genome.fasta",
            filepath=simple_fasta,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, default_settings)
        validator.validate()

        assert len(validator.sequences) == 2
        assert validator.sequences[0].id == "seq1"
        assert validator.sequences[1].id == "seq2"
        assert str(validator.sequences[0].seq) == "ATCGATCGATCGATCGATCG"

    def test_parse_genbank(self, simple_genbank, output_dir, default_settings):
        """Test parsing a GenBank file."""
        genome_config = GenomeConfig(
            filename="genome.gbk",
            filepath=simple_genbank,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.GENBANK
        )

        validator = GenomeValidator(genome_config, output_dir, default_settings)
        validator.validate()

        assert len(validator.sequences) == 1
        assert validator.sequences[0].id == "TEST001"

    def test_parse_empty_file_raises_error(self, empty_fasta, output_dir, default_settings):
        """Test that empty FASTA file raises GenomeValidationError."""
        genome_config = GenomeConfig(
            filename="empty.fasta",
            filepath=empty_fasta,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, default_settings)

        with pytest.raises(GenomeValidationError, match="No sequences found"):
            validator.validate()

    def test_parse_invalid_fasta_raises_error(self, invalid_fasta, output_dir, default_settings):
        """Test that invalid FASTA raises FastaFormatError."""
        genome_config = GenomeConfig(
            filename="invalid.fasta",
            filepath=invalid_fasta,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, default_settings)

        with pytest.raises((FastaFormatError, GenomeValidationError)):
            validator.validate()


class TestGenomeValidatorCompression:
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
    def default_settings(self):
        """Default settings that don't filter sequences."""
        return GenomeValidator.Settings(min_sequence_length=0)

    @pytest.fixture
    def compressed_fasta_gz(self, temp_dir):
        """Create a gzip compressed FASTA file."""
        fasta_file = temp_dir / "genome.fasta.gz"
        with gzip.open(fasta_file, "wt") as f:
            f.write(">seq1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
        return fasta_file

    @pytest.fixture
    def compressed_fasta_bz2(self, temp_dir):
        """Create a bzip2 compressed FASTA file."""
        fasta_file = temp_dir / "genome.fasta.bz2"
        with bz2.open(fasta_file, "wt") as f:
            f.write(">seq1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
        return fasta_file

    def test_parse_gzip_compressed(self, compressed_fasta_gz, output_dir, default_settings):
        """Test parsing gzip compressed FASTA."""
        genome_config = GenomeConfig(
            filename="genome.fasta.gz",
            filepath=compressed_fasta_gz,
            coding_type=CodingType.GZIP,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, default_settings)
        validator.validate()

        assert len(validator.sequences) == 1
        assert str(validator.sequences[0].seq) == "ATCGATCGATCGATCGATCG"

    def test_parse_bzip2_compressed(self, compressed_fasta_bz2, output_dir, default_settings):
        """Test parsing bzip2 compressed FASTA."""
        genome_config = GenomeConfig(
            filename="genome.fasta.bz2",
            filepath=compressed_fasta_bz2,
            coding_type=CodingType.BZIP2,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, default_settings)
        validator.validate()

        assert len(validator.sequences) == 1
        assert str(validator.sequences[0].seq) == "ATCGATCGATCGATCGATCG"


class TestGenomeValidatorValidation:
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
    def default_settings(self):
        """Default settings that don't filter sequences."""
        return GenomeValidator.Settings(min_sequence_length=0)

    @pytest.fixture
    def fasta_with_duplicates(self, temp_dir):
        """Create FASTA with duplicate sequence IDs."""
        fasta_file = temp_dir / "duplicates.fasta"
        with open(fasta_file, "w") as f:
            f.write(">seq1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write(">seq1\n")
            f.write("GCTAGCTAGCTAGCTAGCTA\n")
        return fasta_file

    @pytest.fixture
    def fasta_with_invalid_chars(self, temp_dir):
        """Create FASTA with invalid characters."""
        fasta_file = temp_dir / "invalid_chars.fasta"
        with open(fasta_file, "w") as f:
            f.write(">seq1\n")
            f.write("ATCGXYZATCG\n")  # X, Y, Z are invalid
        return fasta_file

    @pytest.fixture
    def fasta_with_empty_sequence(self, temp_dir):
        """Create FASTA with an empty sequence."""
        fasta_file = temp_dir / "empty_seq.fasta"
        with open(fasta_file, "w") as f:
            f.write(">seq1\n")
            f.write("\n")
            f.write(">seq2\n")
            f.write("ATCGATCG\n")
        return fasta_file

    def test_duplicate_ids_allowed_by_default(self, fasta_with_duplicates, output_dir, default_settings):
        """Test that duplicate IDs are allowed by default."""
        genome_config = GenomeConfig(
            filename="duplicates.fasta",
            filepath=fasta_with_duplicates,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        # Default settings allow duplicates
        validator = GenomeValidator(genome_config, output_dir, default_settings)
        validator.validate()  # Should not raise

        assert len(validator.sequences) == 2

    def test_duplicate_ids_rejected_when_disabled(self, fasta_with_duplicates, output_dir):
        """Test that duplicate IDs are rejected when allow_duplicate_ids=False."""
        settings = GenomeValidator.Settings(allow_duplicate_ids=False)

        genome_config = GenomeConfig(
            filename="duplicates.fasta",
            filepath=fasta_with_duplicates,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, settings)

        with pytest.raises(GenomeValidationError, match="Duplicate sequence IDs"):
            validator.validate()

    def test_invalid_chars_detected(self, fasta_with_invalid_chars, output_dir):
        """Test that invalid characters are detected when check_invalid_chars=True."""
        settings = GenomeValidator.Settings(check_invalid_chars=True)

        genome_config = GenomeConfig(
            filename="invalid_chars.fasta",
            filepath=fasta_with_invalid_chars,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, settings)

        with pytest.raises(GenomeValidationError, match="invalid chars"):
            validator.validate()

    def test_invalid_chars_ignored_by_default(self, fasta_with_invalid_chars, output_dir, default_settings):
        """Test that invalid characters are ignored by default."""
        genome_config = GenomeConfig(
            filename="invalid_chars.fasta",
            filepath=fasta_with_invalid_chars,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        # Default settings don't check invalid chars
        validator = GenomeValidator(genome_config, output_dir, default_settings)
        validator.validate()  # Should not raise

    def test_empty_sequence_rejected(self, fasta_with_empty_sequence, output_dir, default_settings):
        """Test that empty sequences are rejected by default."""
        genome_config = GenomeConfig(
            filename="empty_seq.fasta",
            filepath=fasta_with_empty_sequence,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, default_settings)

        with pytest.raises(GenomeValidationError, match="zero length"):
            validator.validate()


class TestGenomeValidatorEditing:
    """Test editing specifications."""

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
    def fasta_with_mixed_lengths(self, temp_dir):
        """Create FASTA with sequences of various lengths."""
        fasta_file = temp_dir / "mixed_lengths.fasta"
        with open(fasta_file, "w") as f:
            f.write(">short_seq\n")
            f.write("ATCG\n")  # 4bp
            f.write(">medium_seq\n")
            f.write("ATCGATCG" * 10 + "\n")  # 80bp
            f.write(">long_seq\n")
            f.write("ATCGATCG" * 50 + "\n")  # 400bp
        return fasta_file

    def test_min_sequence_length_filter(self, fasta_with_mixed_lengths, output_dir):
        """Test filtering sequences by minimum length."""
        settings = GenomeValidator.Settings(min_sequence_length=100)

        genome_config = GenomeConfig(
            filename="mixed_lengths.fasta",
            filepath=fasta_with_mixed_lengths,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, settings)
        validator.validate()

        # Only long_seq (400bp) should remain
        assert len(validator.sequences) == 1
        assert validator.sequences[0].id == "long_seq"

    def test_sequence_prefix_added(self, fasta_with_mixed_lengths, output_dir):
        """Test adding prefix to sequence IDs."""
        settings = GenomeValidator.Settings(
            sequence_prefix="chr",
            min_sequence_length=0
        )

        genome_config = GenomeConfig(
            filename="mixed_lengths.fasta",
            filepath=fasta_with_mixed_lengths,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, settings)
        validator.validate()

        # All sequences should have prefix
        assert all(seq.id.startswith("chr_") for seq in validator.sequences)
        assert "chr_short_seq" in [seq.id for seq in validator.sequences]
        assert "chr_medium_seq" in [seq.id for seq in validator.sequences]
        assert "chr_long_seq" in [seq.id for seq in validator.sequences]


class TestGenomeValidatorStatistics:
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
    def default_settings(self):
        """Default settings that don't filter sequences."""
        return GenomeValidator.Settings(min_sequence_length=0)

    @pytest.fixture
    def fasta_for_stats(self, temp_dir):
        """Create FASTA with known statistics."""
        fasta_file = temp_dir / "stats.fasta"
        with open(fasta_file, "w") as f:
            f.write(">seq1\n")
            f.write("ATCGATCG\n")  # 8bp, 50% GC
            f.write(">seq2\n")
            f.write("GCGCGCGCGCGC\n")  # 12bp, 100% GC
        return fasta_file

    def test_statistics_collection(self, fasta_for_stats, output_dir, default_settings):
        """Test that statistics are collected correctly."""
        genome_config = GenomeConfig(
            filename="stats.fasta",
            filepath=fasta_for_stats,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, default_settings)
        validator.validate()

        stats = validator.get_statistics()

        assert stats['num_sequences'] == 2
        assert stats['total_length'] == 20
        assert stats['min_length'] == 8
        assert stats['max_length'] == 12
        assert stats['avg_length'] == 10.0
        # GC content: (4 + 12) / 20 = 80%
        assert 79 < stats['gc_content'] < 81


class TestGenomeValidatorOutput:
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
    def default_settings(self):
        """Default settings that don't filter sequences."""
        return GenomeValidator.Settings(min_sequence_length=0)

    @pytest.fixture
    def simple_genbank(self, temp_dir):
        """Create a simple GenBank file for conversion test."""
        gb_file = temp_dir / "genome.gbk"
        record = SeqRecord(
            Seq("ATCGATCGATCGATCGATCG"),
            id="TEST001",
            name="test_sequence",
            description="Test sequence",
        )
        # GenBank requires molecule_type annotation
        record.annotations["molecule_type"] = "DNA"
        SeqIO.write([record], gb_file, "genbank")
        return gb_file

    @pytest.fixture
    def simple_fasta(self, temp_dir):
        """Create a simple FASTA file."""
        fasta_file = temp_dir / "genome.fasta"
        with open(fasta_file, "w") as f:
            f.write(">seq1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
        return fasta_file

    def test_output_uncompressed(self, simple_fasta, output_dir, default_settings):
        """Test generating uncompressed output."""
        settings = GenomeValidator.Settings(coding_type=None, min_sequence_length=0)

        genome_config = GenomeConfig(
            filename="genome.fasta",
            filepath=simple_fasta,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, settings)
        validator.validate()

        # Check output file exists
        output_files = list(output_dir.glob("*.fasta"))
        assert len(output_files) == 1
        assert output_files[0].suffix == ".fasta"

    def test_output_gzip_compressed(self, simple_fasta, output_dir, default_settings):
        """Test generating gzip compressed output."""
        settings = GenomeValidator.Settings(coding_type='gz', min_sequence_length=0)

        genome_config = GenomeConfig(
            filename="genome.fasta",
            filepath=simple_fasta,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, settings)
        validator.validate()

        # Check gzip output file exists
        output_files = list(output_dir.glob("*.fasta.gz"))
        assert len(output_files) == 1

        # Verify it's actually gzipped and contains correct data
        with gzip.open(output_files[0], 'rt') as f:
            records = list(SeqIO.parse(f, 'fasta'))
            assert len(records) == 1

    def test_output_bzip2_compressed(self, simple_fasta, output_dir, default_settings):
        """Test generating bzip2 compressed output."""
        settings = GenomeValidator.Settings(coding_type='bz2', min_sequence_length=0)

        genome_config = GenomeConfig(
            filename="genome.fasta",
            filepath=simple_fasta,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, settings)
        validator.validate()

        # Check bzip2 output file exists
        output_files = list(output_dir.glob("*.fasta.bz2"))
        assert len(output_files) == 1

        # Verify it's actually bzip2'd and contains correct data
        with bz2.open(output_files[0], 'rt') as f:
            records = list(SeqIO.parse(f, 'fasta'))
            assert len(records) == 1

    def test_output_with_suffix(self, simple_fasta, output_dir, default_settings):
        """Test output filename with custom suffix."""
        settings = GenomeValidator.Settings(output_filename_suffix="processed", min_sequence_length=0)

        genome_config = GenomeConfig(
            filename="genome.fasta",
            filepath=simple_fasta,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, settings)
        validator.validate()

        output_files = list(output_dir.glob("*_processed.fasta"))
        assert len(output_files) == 1

    def test_output_with_subdirectory(self, simple_fasta, output_dir, default_settings):
        """Test output to custom subdirectory."""
        settings = GenomeValidator.Settings(output_subdir_name="genomes", min_sequence_length=0)

        genome_config = GenomeConfig(
            filename="genome.fasta",
            filepath=simple_fasta,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
        )

        validator = GenomeValidator(genome_config, output_dir, settings)
        validator.validate()

        subdir = output_dir / "genomes"
        assert subdir.exists()
        output_files = list(subdir.glob("*.fasta"))
        assert len(output_files) == 1

    def test_genbank_to_fasta_conversion(self, simple_genbank, output_dir, default_settings):
        """Test converting GenBank to FASTA output."""
        genome_config = GenomeConfig(
            filename="genome.gbk",
            filepath=simple_genbank,
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.GENBANK
        )

        validator = GenomeValidator(genome_config, output_dir, default_settings)
        validator.validate()

        # Output should be FASTA
        output_files = list(output_dir.glob("*.fasta"))
        assert len(output_files) == 1

        # Verify content
        with open(output_files[0], 'r') as f:
            records = list(SeqIO.parse(f, 'fasta'))
            assert len(records) == 1
            assert records[0].id == "TEST001"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
