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

from validation_pkg.coordinator import GenomeConfig
from validation_pkg.validators.genome_validator import GenomeValidator
from validation_pkg.exceptions import (
    GenomeValidationError,
    FastaFormatError,
    GenBankFormatError,
    CompressionError,
    FileNotFoundError as ValidationFileNotFoundError
)
from validation_pkg.utils import settings


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
        out_dir = temp_dir / settings.output_base_dir
        out_dir.mkdir()
        return out_dir
    
    @pytest.fixture
    def simple_fasta(self, temp_dir):
        """Create a simple FASTA file."""
        fasta_file = temp_dir / "genome.fasta"
        with open(fasta_file, 'w') as f:
            f.write(">seq1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
            f.write(">seq2\n")
            f.write("GCTAGCTAGCTAGCTAGCTA\n")
        return fasta_file
    
    @pytest.fixture
    def simple_genbank(self, temp_dir):
        """Create a simple GenBank file."""
        gb_file = temp_dir / "genome.gbk"
        # Create a simple GenBank record using BioPython
        record = SeqRecord(
            Seq("ATCGATCGATCGATCGATCG"),
            id="TEST001",
            name="test_sequence",
            description="Test sequence"
        )
        SeqIO.write([record], gb_file, "genbank")
        return gb_file
    
    @pytest.fixture
    def compressed_fasta_gz(self, temp_dir):
        """Create a gzip compressed FASTA file."""
        fasta_file = temp_dir / "genome.fasta.gz"
        with gzip.open(fasta_file, 'wt') as f:
            f.write(">seq1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
        return fasta_file
    
    @pytest.fixture
    def compressed_fasta_bz2(self, temp_dir):
        """Create a bzip2 compressed FASTA file."""
        fasta_file = temp_dir / "genome.fasta.bz2"
        with bz2.open(fasta_file, 'wt') as f:
            f.write(">seq1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
        return fasta_file
    
    @pytest.fixture
    def fasta_with_short_sequences(self, temp_dir):
        """Create FASTA with some short sequences."""
        fasta_file = temp_dir / "mixed_lengths.fasta"
        with open(fasta_file, 'w') as f:
            f.write(">short_seq\n")
            f.write("ATCG\n")  # Only 4bp - will be removed
            f.write(">long_seq\n")
            f.write("ATCGATCGATCGATCGATCG" * 10 + "\n")  # 200bp
        return fasta_file
    
    @pytest.fixture
    def fasta_with_ambiguous(self, temp_dir):
        """Create FASTA with ambiguous nucleotides."""
        fasta_file = temp_dir / "ambiguous.fasta"
        with open(fasta_file, 'w') as f:
            f.write(">seq1\n")
            f.write("ATCGNNNNNATCG\n")  # 5/13 = ~38% ambiguous
        return fasta_file
    
    @pytest.fixture
    def invalid_fasta(self, temp_dir):
        """Create an invalid FASTA file."""
        fasta_file = temp_dir / "invalid.fasta"
        with open(fasta_file, 'w') as f:
            f.write("This is not a valid FASTA file\n")
            f.write("No headers at all\n")
        return fasta_file
    
    def test_init_validator(self, simple_fasta, output_dir):
        """Test GenomeValidator initialization."""
        genome_config = GenomeConfig(
            filename="genome.fasta",
            filepath=simple_fasta
        )
        
        validator = GenomeValidator(genome_config,output_dir,settings.genome_settings)
        
        assert validator.input_path == simple_fasta
        assert validator.output_dir.exists()
        assert validator.detected_format is None
        assert len(validator.sequences) == 0
    
    def test_file_not_found(self, temp_dir, output_dir):
        """Test error when genome file doesn't exist."""
        genome_config = GenomeConfig(
            filename="nonexistent.fasta",
            filepath=temp_dir / "nonexistent.fasta"
        )
        
        validator = GenomeValidator(genome_config,output_dir,settings.genome_settings)
        
        with pytest.raises(ValidationFileNotFoundError):
            validator.validate()
    
    def test_detect_fasta_format(self, simple_fasta, output_dir):
        """Test FASTA format detection."""
        genome_config = GenomeConfig(
            filename="genome.fasta",
            filepath=simple_fasta
        )
        
        validator = GenomeValidator(genome_config, output_dir,settings.genome_settings)
        validator._check_file_exists()
        validator._detect_compression()
        validator._detect_format()
        
        assert validator.detected_format == "fasta"
        assert validator.is_compressed == False
    
    # def test_detect_genbank_format(self, simple_genbank, output_dir):
    #     """Test GenBank format detection."""
    #     genome_config = GenomeConfig(
    #         filename="genome.gbk",
    #         filepath=simple_genbank
    #     )
        
    #     validator = GenomeValidator(genome_config, output_dir,globals.genome_settings)
    #     validator._check_file_exists()
    #     validator._detect_compression()
    #     validator._detect_format()
        
    #     assert validator.detected_format == "genbank"
    
    def test_detect_gzip_compression(self, compressed_fasta_gz, output_dir):
        """Test gzip compression detection."""
        genome_config = GenomeConfig(
            filename="genome.fasta.gz",
            filepath=compressed_fasta_gz
        )
        
        validator = GenomeValidator(genome_config, output_dir,settings.genome_settings)
        validator._check_file_exists()
        validator._detect_compression()
        
        assert validator.is_compressed == True
        assert validator.compression_type == "gz"
    
    def test_detect_bzip2_compression(self, compressed_fasta_bz2, output_dir):
        """Test bzip2 compression detection."""
        genome_config = GenomeConfig(
            filename="genome.fasta.bz2",
            filepath=compressed_fasta_bz2
        )
        
        validator = GenomeValidator(genome_config, output_dir,settings.genome_settings)
        validator._check_file_exists()
        validator._detect_compression()
        
        assert validator.is_compressed == True
        assert validator.compression_type == "bz2"
    
    def test_parse_fasta_file(self, simple_fasta, output_dir):
        """Test parsing FASTA file."""
        genome_config = GenomeConfig(
            filename="genome.fasta",
            filepath=simple_fasta
        )
        
        validator = GenomeValidator(genome_config, output_dir,settings.genome_settings)
        validator._check_file_exists()
        validator._detect_compression()
        validator._detect_format()
        validator._parse_file()
        
        assert len(validator.sequences) == 2
        assert validator.sequences[0].id == "seq1"
        assert validator.sequences[1].id == "seq2"
    
    # def test_parse_genbank_file(self, simple_genbank, output_dir):
    #     """Test parsing GenBank file."""
    #     genome_config = GenomeConfig(
    #         filename="genome.gbk",
    #         filepath=simple_genbank
    #     )
        
    #     validator = GenomeValidator(genome_config, output_dir)
    #     validator._check_file_exists()
    #     validator._detect_compression()
    #     validator._detect_format()
    #     validator._parse_file()
        
    #     assert len(validator.sequences) == 1
    #     assert validator.sequences[0].id == "TEST001"
    
    def test_parse_compressed_fasta(self, compressed_fasta_gz, output_dir):
        """Test parsing compressed FASTA file."""
        genome_config = GenomeConfig(
            filename="genome.fasta.gz",
            filepath=compressed_fasta_gz
        )
        
        validator = GenomeValidator(genome_config, output_dir,settings.genome_settings)
        validator._check_file_exists()
        validator._detect_compression()
        validator._detect_format()
        validator._parse_file()
        
        assert len(validator.sequences) == 1
    
    def test_invalid_fasta_format(self, invalid_fasta, output_dir):
        """Test error with invalid FASTA format."""
        genome_config = GenomeConfig(
            filename="invalid.fasta",
            filepath=invalid_fasta
        )
        
        validator = GenomeValidator(genome_config, output_dir,settings.genome_settings)
        validator._check_file_exists()
        validator._detect_compression()
        validator._detect_format()
        
        with pytest.raises(FastaFormatError):
            validator._parse_file()
    
    # def test_genbank_to_fasta_conversion(self, simple_genbank, output_dir):
    #     """Test GenBank to FASTA conversion."""
    #     genome_config = GenomeConfig(
    #         filename="genome.gbk",
    #         filepath=simple_genbank
    #     )
        
    #     validator = GenomeValidator(genome_config, output_dir)
    #     validator.validate()
        
    #     # Check output is FASTA
    #     output_files = list((output_dir / "output" / "genomes").glob("*.fasta"))
    #     assert len(output_files) == 1
        
    #     # Verify it's valid FASTA
    #     output_file = output_files[0]
    #     with open(output_file) as f:
    #         records = list(SeqIO.parse(f, "fasta"))
    #     assert len(records) == 1
    
    def test_output_compression_none(self, simple_fasta, output_dir):
        """Test output without compression."""
        genome_config = GenomeConfig(
            filename="genome.fasta",
            filepath=simple_fasta
        )
        
        settings = settings.genome_settings
        settings["output"]['coding_type'] = None
        
        validator = GenomeValidator(genome_config, output_dir, settings)
        validator.validate()
        
        # Check uncompressed output exists
        output_files = list((output_dir / "output" / "genomes").glob("genome.fasta"))
        assert len(output_files) == 1
    
    # def test_output_compression_gz(self, simple_fasta, output_dir):
    #     """Test output with gzip compression."""
    #     genome_config = GenomeConfig(
    #         filename="genome.fasta",
    #         filepath=simple_fasta
    #     )
        
    #     settings = ValidatorSettings()
    #     settings.genome_settings["output"]['coding_type'] = 'gz'
        
    #     validator = GenomeValidator(genome_config, output_dir, settings)
    #     validator.validate()
        
    #     # Check compressed output exists
    #     output_files = list((output_dir / "output" / "genomes").glob("genome.fasta.gz"))
    #     assert len(output_files) == 1
        
    #     # Verify it can be decompressed and read
    #     output_file = output_files[0]
    #     with gzip.open(output_file, 'rt') as f:
    #         records = list(SeqIO.parse(f, "fasta"))
    #     assert len(records) > 0
    
    # def test_output_compression_bzip(self, simple_fasta, output_dir):
    #     """Test output with bzip2 compression."""
    #     genome_config = GenomeConfig(
    #         filename="genome.fasta",
    #         filepath=simple_fasta
    #     )
        
    #     settings = ValidatorSettings()
    #     settings.genome_settings["output"]['coding_type'] = 'bzip'
        
    #     validator = GenomeValidator(genome_config, output_dir, settings)
    #     validator.validate()
        
    #     # Check compressed output exists
    #     output_files = list((output_dir / "output" / "genomes").glob("genome.fasta.bz2"))
    #     assert len(output_files) == 1
    
    # def test_duplicate_sequence_ids_warning(self, temp_dir, output_dir):
    #     """Test warning for duplicate sequence IDs."""
    #     fasta_file = temp_dir / "duplicates.fasta"
    #     with open(fasta_file, 'w') as f:
    #         f.write(">seq1\n")
    #         f.write("ATCGATCG\n")
    #         f.write(">seq1\n")  # Duplicate ID
    #         f.write("GCTAGCTA\n")
        
    #     genome_config = GenomeConfig(
    #         filename="duplicates.fasta",
    #         filepath=fasta_file
    #     )
        
    #     validator = GenomeValidator(genome_config, output_dir)
    #     validator.validate()
        
    #     # Should complete but with warnings
    #     # (Check logs would show warnings about duplicates)
    #     assert len(validator.sequences) > 0
    
    # def test_empty_sequence_error(self, temp_dir, output_dir):
    #     """Test error with empty sequence."""
    #     fasta_file = temp_dir / "empty_seq.fasta"
    #     with open(fasta_file, 'w') as f:
    #         f.write(">seq1\n")
    #         f.write("\n")  # Empty sequence
        
    #     genome_config = GenomeConfig(
    #         filename="empty_seq.fasta",
    #         filepath=fasta_file
    #     )
        
    #     validator = GenomeValidator(genome_config, output_dir)
    #     validator._check_file_exists()
    #     validator._detect_compression()
    #     validator._detect_format()
    #     validator._parse_file()
        
    #     with pytest.raises(GenomeValidationError, match="zero length"):
    #         validator._validate_sequences()
    
    # def test_multiple_file_extensions(self, temp_dir, output_dir):
    #     """Test handling of different FASTA extensions."""
    #     for ext in ['.fa', '.fasta', '.fna']:
    #         fasta_file = temp_dir / f"genome{ext}"
    #         with open(fasta_file, 'w') as f:
    #             f.write(">seq1\n")
    #             f.write("ATCGATCG\n")
            
    #         genome_config = GenomeConfig(
    #             filename=f"genome{ext}",
    #             filepath=fasta_file
    #         )
            
    #         validator = GenomeValidator(genome_config, output_dir)
    #         validator._detect_format()
            
    #         assert validator.detected_format == "fasta"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])