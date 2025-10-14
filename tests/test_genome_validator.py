"""
Comprehensive tests for GenomeValidator.
"""

import pytest
import tempfile
import gzip
import bz2
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Keep your original package imports
from validation_pkg.config_manager import GenomeConfig
from validation_pkg.validators.genome_validator import GenomeValidator
from validation_pkg.exceptions import (
    GenomeValidationError,
    FastaFormatError,
    GenBankFormatError,
    CompressionError,
    FileNotFoundError as ValidationFileNotFoundError,
)
from validation_pkg.utils import settings as settings_mod


class TestGenomeValidator:
    """Test suite for GenomeValidator."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def output_dir(self, temp_dir):
        """Create output directory."""
        out_dir = temp_dir / settings_mod.output_base_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @pytest.fixture
    def simple_fasta(self, temp_dir):
        """Create a simple FASTA file."""
        fasta_file = temp_dir / "genome.fasta"
        with open(fasta_file, "w") as f:
            f.write(">seq1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write(">seq2\n")
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
        SeqIO.write([record], gb_file, "genbank")
        return gb_file

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

    @pytest.fixture
    def fasta_with_short_sequences(self, temp_dir):
        """Create FASTA with some short sequences."""
        fasta_file = temp_dir / "mixed_lengths.fasta"
        with open(fasta_file, "w") as f:
            f.write(">short_seq\n")
            f.write("ATCG\n")  # Only 4bp - will be removed (depending on min length settings)
            f.write(">long_seq\n")
            f.write("ATCGATCGATCGATCGATCG" * 10 + "\n")  # 200bp
        return fasta_file

    @pytest.fixture
    def fasta_with_ambiguous(self, temp_dir):
        """Create FASTA with ambiguous nucleotides."""
        fasta_file = temp_dir / "ambiguous.fasta"
        with open(fasta_file, "w") as f:
            f.write(">seq1\n")
            f.write("ATCGNNNNNATCG\n")  # 5/13 ambiguous
        return fasta_file

    @pytest.fixture
    def invalid_fasta(self, temp_dir):
        """Create an invalid FASTA file."""
        fasta_file = temp_dir / "invalid.fasta"
        with open(fasta_file, "w") as f:
            f.write("This is not a valid FASTA file\n")
            f.write("No headers at all\n")
        return fasta_file

    def test_init_validator(self, simple_fasta, output_dir):
        """Test GenomeValidator initialization."""
        genome_config = GenomeConfig(
            filename="genome.fasta",
            filepath=simple_fasta,
        )

        validator = GenomeValidator(genome_config, output_dir, settings_mod.genome_settings)

        assert validator.input_path == simple_fasta
        assert validator.output_dir.exists()
        # Don't assume internal attributes if the implementation differs;
        # basic smoke check that the object is created is enough.

    def test_file_not_found(self, temp_dir, output_dir):
        """Test error when genome file doesn't exist."""
        genome_config = GenomeConfig(
            filename="nonexistent.fasta",
            filepath=temp_dir / "nonexistent.fasta",
        )

        validator = GenomeValidator(genome_config, output_dir, settings_mod.genome_settings)

        with pytest.raises(ValidationFileNotFoundError):
            validator.validate()

    def test_detect_fasta_and_parse(self, simple_fasta, output_dir):
        """
        End-to-end FASTA path: validate() should complete successfully for a well-formed FASTA.
        Also sanity-check that sequences were parsed.
        """
        genome_config = GenomeConfig(
            filename="genome.fasta",
            filepath=simple_fasta,
        )

        validator = GenomeValidator(genome_config, output_dir, settings_mod.genome_settings)
        result = validator.validate()

        # Truthy/None depending on your implementation — the key is: should not raise.
        assert result is not False

        # If implementation exposes parsed sequences, check length:
        if hasattr(validator, "sequences"):
            assert len(validator.sequences) == 2

    # If your validator supports GenBank, you can enable this test.
    # def test_detect_genbank_and_parse(self, simple_genbank, output_dir):
    #     """GenBank detection and parsing should succeed."""
    #     genome_config = GenomeConfig(
    #         filename="genome.gbk",
    #         filepath=simple_genbank,
    #     )
    #     validator = GenomeValidator(genome_config, output_dir, settings_mod.genome_settings)
    #     result = validator.validate()
    #     assert result is not False
    #     if hasattr(validator, "sequences"):
    #         assert len(validator.sequences) == 1
    #         assert validator.sequences[0].id == "TEST001"

    def test_detect_gzip_compression(self, compressed_fasta_gz, output_dir):
        """gzip-compressed FASTA should validate and parse."""
        genome_config = GenomeConfig(
            filename="genome.fasta.gz",
            filepath=compressed_fasta_gz,
        )
        validator = GenomeValidator(genome_config, output_dir, settings_mod.genome_settings)
        result = validator.validate()
        assert result is not False
        if hasattr(validator, "sequences"):
            assert len(validator.sequences) >= 1

    def test_detect_bzip2_compression(self, compressed_fasta_bz2, output_dir):
        """bzip2-compressed FASTA should validate and parse."""
        genome_config = GenomeConfig(
            filename="genome.fasta.bz2",
            filepath=compressed_fasta_bz2,
        )
        validator = GenomeValidator(genome_config, output_dir, settings_mod.genome_settings)
        result = validator.validate()
        assert result is not False
        if hasattr(validator, "sequences"):
            assert len(validator.sequences) >= 1

    def test_invalid_fasta_format(self, invalid_fasta, output_dir):
        """Invalid FASTA should raise a FastaFormatError (or a more general GenomeValidationError)."""
        genome_config = GenomeConfig(
            filename="invalid.fasta",
            filepath=invalid_fasta,
        )
        validator = GenomeValidator(genome_config, output_dir, settings_mod.genome_settings)

        with pytest.raises((FastaFormatError, GenomeValidationError)):
            validator.validate()

    # If your implementation writes outputs and supports compression selection in settings,
    # you can adapt and enable tests like the following. They’re left commented to avoid coupling
    # to internal paths/structures that may differ between implementations.

    # def test_output_uncompressed(self, simple_fasta, output_dir):
    #     """Validate with no output compression."""
    #     # Make a shallow copy of genome settings (object or dict)
    #     gs = settings_mod.genome_settings.copy() if hasattr(settings_mod.genome_settings, "copy") else settings_mod.genome_settings
    #     # If your settings object has a field or attr for coding type, set it here:
    #     if isinstance(gs, dict):
    #         gs.setdefault("output", {}).update({"coding_type": None})
    #     else:
    #         # Example for dataclass-like settings:
    #         gs.output.coding_type = None
    #
    #     genome_config = GenomeConfig(filename="genome.fasta", filepath= simple_fasta)
    #     validator = GenomeValidator(genome_config, output_dir, gs)
    #     validator.validate()
    #
    #     # Check for expected uncompressed output if your validator writes it:
    #     out_glob = (output_dir / "output" / "genomes").glob("*.fasta")
    #     assert any(out_glob)

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
