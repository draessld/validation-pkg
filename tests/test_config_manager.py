"""
Comprehensive tests for Config and ConfigManager classes.
"""

import pytest
import json
import tempfile
from pathlib import Path
from validation_pkg.config_manager import ConfigManager, Config, GenomeConfig, ReadConfig, FeatureConfig
from validation_pkg.exceptions import ConfigurationError, FileNotFoundError as ValidationFileNotFoundError
from validation_pkg.utils.formats import ReadFormat,FeatureFormat,GenomeFormat,CodingType


class TestConfigManager:
    """Test suite for ConfigManager configuration loading and validation."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    @pytest.fixture
    def valid_config_minimal(self, temp_dir):
        """Create minimal valid configuration with actual files."""
        # Create dummy files
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads_R1.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [
                {"filename": "reads_R1.fastq", "ngs_type": "illumina"}
            ]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))
        return config_file
    
    @pytest.fixture
    def valid_config_full(self, temp_dir):
        """Create full valid configuration with all optional fields."""
        # Create dummy files
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "ref_plasmid.gbk").write_text("LOCUS\n")
        (temp_dir / "mod_plasmid.gbk").write_text("LOCUS\n")
        (temp_dir / "reads_R1.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        (temp_dir / "reads_R2.fastq").write_text("@read2\nATCG\n+\nIIII\n")
        (temp_dir / "ref_features.gff").write_text("##gff-version 3\n")
        (temp_dir / "mod_features.bed").write_text("chr1\t100\t200\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "ref_plasmid_filename": {"filename": "ref_plasmid.gbk"},
            "mod_plasmid_filename": {"filename": "mod_plasmid.gbk"},
            "reads": [
                {"filename": "reads_R1.fastq", "ngs_type": "illumina"},
                {"filename": "reads_R2.fastq", "ngs_type": "illumina"}
            ],
            "ref_feature_filename": {"filename": "ref_features.gff"},
            "mod_feature_filename": {"filename": "mod_features.bed"},
            "options": {"threads": 8, "verbose": True}
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))
        return config_file
    
    def test_load_minimal_valid_config(self, valid_config_minimal):
        """Test loading minimal valid configuration."""
        config = ConfigManager.load(str(valid_config_minimal))
        
        assert config.ref_genome is not None
        assert config.mod_genome is not None
        assert len(config.reads) == 1
        assert config.reads[0].ngs_type == "illumina"
        assert config.ref_plasmid is None
        assert config.mod_plasmid is None
        
        # Check absolute paths are set
        assert config.ref_genome.filepath.is_absolute()
        assert config.ref_genome.filepath.exists()
        assert config.mod_genome.filepath.is_absolute()
        assert config.reads[0].filepath.is_absolute()
    
    def test_load_full_valid_config(self, valid_config_full):
        """Test loading full configuration with all fields."""
        config = ConfigManager.load(str(valid_config_full))
        
        assert config.ref_plasmid is not None
        assert config.mod_plasmid is not None
        assert len(config.reads) == 2
        assert config.ref_feature is not None
        assert config.mod_feature is not None
        assert config.options["threads"] == 8
        
        # Check all absolute paths are set
        assert config.ref_genome.filepath.is_absolute()
        assert config.mod_genome.filepath.is_absolute()
        assert config.ref_plasmid.filepath.is_absolute()
        assert config.mod_plasmid.filepath.is_absolute()
        assert config.ref_feature.filepath.is_absolute()
        assert config.mod_feature.filepath.is_absolute()
        assert all(read.filepath.is_absolute() for read in config.reads)
        
    
    def test_missing_config_file(self, temp_dir):
        """Test error when config file doesn't exist."""
        with pytest.raises(ValidationFileNotFoundError, match="Configuration file not found"):
            ConfigManager.load(str(temp_dir / "nonexistent.json"))
    
    def test_missing_required_field_ref_genome(self, temp_dir):
        """Test error when ref_genome_filename is missing."""
        config = {
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ConfigurationError, ValueError), match="Missing required field: ref_genome_filename"):
            ConfigManager.load(str(config_file))
    
    def test_missing_required_field_mod_genome(self, temp_dir):
        """Test error when mod_genome_filename is missing."""
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ConfigurationError, ValueError), match="Missing required field: mod_genome_filename"):
            ConfigManager.load(str(config_file))
    
    def test_missing_required_field_reads(self, temp_dir):
        """Test error when reads field is missing."""
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"}
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ConfigurationError, ValueError), match="Missing required field: reads"):
            ConfigManager.load(str(config_file))
    
    def test_empty_reads_list(self, temp_dir):
        """Test error when reads list is empty."""
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": []
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ConfigurationError, ValueError), match="'reads' must be a non-empty list"):
            ConfigManager.load(str(config_file))
    
    def test_invalid_ngs_type(self, temp_dir):
        """Test error with invalid ngs_type value."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "454"}]  # Invalid type
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ConfigurationError, ValueError), match="Invalid ngs_type"):
            ConfigManager.load(str(config_file))
    
    def test_missing_genome_file(self, temp_dir):
        """Test error when referenced genome file doesn't exist."""
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},  # File doesn't exist
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ValidationFileNotFoundError, FileNotFoundError), match="ref_genome|ref.fasta"):
            ConfigManager.load(str(config_file))
    
    def test_missing_read_file(self, temp_dir):
        """Test error when referenced read file doesn't exist."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "missing_reads.fastq", "ngs_type": "illumina"}]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ValidationFileNotFoundError, FileNotFoundError), match="reads\\[0\\]|missing_reads.fastq"):
            ConfigManager.load(str(config_file))
    
    def test_directory_reads(self, temp_dir):
        """Test reads specified as directory."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        
        # Create reads directory
        reads_dir = temp_dir / "reads"
        reads_dir.mkdir()
        (reads_dir / "R1.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        (reads_dir / "R2.fastq").write_text("@read2\nATCG\n+\nIIII\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"directory": "reads", "ngs_type": "illumina"}]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        loaded_config = ConfigManager.load(str(config_file))
        assert len(loaded_config.reads) == 2
        assert loaded_config.reads[0].filepath == (reads_dir / "R2.fastq")
        assert loaded_config.reads[0].filename == "R2.fastq"
        assert loaded_config.reads[0].detected_format == ReadFormat.FASTQ
        assert loaded_config.reads[0].ngs_type == "illumina"
        assert loaded_config.reads[1].filepath == (reads_dir / "R1.fastq")
        assert loaded_config.reads[1].filename == "R1.fastq"
        assert loaded_config.reads[1].detected_format == ReadFormat.FASTQ
        assert loaded_config.reads[1].ngs_type == "illumina"
    
    def test_missing_directory_reads(self, temp_dir):
        """Test error when reads directory doesn't exist."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"directory": "nonexistent_dir", "ngs_type": "illumina"}]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ValidationFileNotFoundError, FileNotFoundError), match="Directory not found"):
            ConfigManager.load(str(config_file))
    
    def test_neither_filename_nor_directory_in_reads(self, temp_dir):
        """Test error when reads has neither filename nor directory."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"ngs_type": "illumina"}]  # Missing filename/directory
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ConfigurationError, ValueError), match="filename' or 'directory' field is required"):
            ConfigManager.load(str(config_file))
    
    def test_malformed_json(self, temp_dir):
        """Test error with malformed JSON."""
        config_file = temp_dir / "config.json"
        config_file.write_text("{ invalid json }")
        
        with pytest.raises((ConfigurationError, json.JSONDecodeError)):
            ConfigManager.load(str(config_file))
    
    def test_genome_config_missing_filename_field(self, temp_dir):
        """Test error when genome config dict is missing filename."""
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        
        config = {
            "ref_genome_filename": {"some_key": "value"},  # Missing filename
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ConfigurationError, ValueError), match="must contain 'filename' field"):
            ConfigManager.load(str(config_file))
    
    def test_backwards_compatibility_string_filenames(self, temp_dir):
        """Test that plain strings are accepted for backwards compatibility."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        
        config = {
            "ref_genome_filename": "ref.fasta",  # Plain string
            "mod_genome_filename": "mod.fasta",  # Plain string
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        loaded_config = ConfigManager.load(str(config_file))
        assert loaded_config.ref_genome.filename == "ref.fasta"
        assert loaded_config.ref_genome.filepath.is_absolute()
    
    def test_path_resolution(self, valid_config_minimal):
        """Test that paths are resolved relative to config directory."""
        config = ConfigManager.load(str(valid_config_minimal))
        
        # Check that filepath is in the same directory as config
        config_dir = valid_config_minimal.parent
        assert config.ref_genome.filepath.parent == config_dir
        assert config.config_dir == config_dir
    
    def test_multiple_reads_different_types(self, temp_dir):
        """Test configuration with multiple read types."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "illumina.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        (temp_dir / "ont.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        (temp_dir / "pacbio.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [
                {"filename": "illumina.fastq", "ngs_type": "illumina"},
                {"filename": "ont.fastq", "ngs_type": "ont"},
                {"filename": "pacbio.fastq", "ngs_type": "pacbio"}
            ]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        loaded_config = ConfigManager.load(str(config_file))
        assert len(loaded_config.reads) == 3
        assert loaded_config.reads[0].ngs_type == "illumina"
        assert loaded_config.reads[1].ngs_type == "ont"
        assert loaded_config.reads[2].ngs_type == "pacbio"
        
        # Check all paths are absolute
        assert all(read.filepath.is_absolute() for read in loaded_config.reads)
    
    def test_extra_keys_stored(self, temp_dir):
        """Test that extra keys in config are stored in extra dict."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        
        config = {
            "ref_genome_filename": {
                "filename": "ref.fasta",
                "custom_key": "custom_value",
                "quality_threshold": 30
            },
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [
                {
                    "filename": "reads.fastq",
                    "ngs_type": "illumina",
                    "custom_read_key": "value"
                }
            ]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        loaded_config = ConfigManager.load(str(config_file))
        
        # Check extra keys are stored
        assert "custom_key" in loaded_config.ref_genome.extra
        assert loaded_config.ref_genome.extra["custom_key"] == "custom_value"
        assert loaded_config.ref_genome.extra["quality_threshold"] == 30
        
        assert "custom_read_key" in loaded_config.reads[0].extra
        assert loaded_config.reads[0].extra["custom_read_key"] == "value"

    def test_genome_file_genbank_format(self,temp_dir):
        """Genome files in GenBank format (.gbk) are accepted and paths resolved."""
        (temp_dir / "ref.gbk").write_text("LOCUS       X  10 bp DNA\nORIGIN\natcg\n//\n")
        (temp_dir / "mod.gbk").write_text("LOCUS       Y  10 bp DNA\nORIGIN\ngcta\n//\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.gbk"},
            "mod_genome_filename": {"filename": "mod.gbk"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        cfg_file = temp_dir / "config.json"
        cfg_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(cfg_file))
        assert loaded.ref_genome.detected_format == GenomeFormat.GENBANK
        assert loaded.mod_genome.detected_format == GenomeFormat.GENBANK
        assert loaded.ref_genome.filepath.is_absolute()
        assert loaded.mod_genome.filepath.is_absolute()

    def test_genome_file_fasta_format(self,temp_dir):
        """Genome files in FASTA format (.fasta) are accepted and paths resolved."""
        (temp_dir / "ref.fasta").write_text(">seq\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        cfg_file = temp_dir / "config.json"
        cfg_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(cfg_file))
        assert loaded.ref_genome.detected_format == GenomeFormat.FASTA
        assert loaded.mod_genome.detected_format == GenomeFormat.FASTA
        assert loaded.ref_genome.filepath.is_absolute()
        assert loaded.mod_genome.filepath.is_absolute()

    def test_genome_file_gz_coding(self,temp_dir):
        """
        A genome entry that requests gz output coding should store that preference.
        (The input file can be any extension; this tests config parsing, not IO.)
        """
        (temp_dir / "ref.fasta.gz").write_text(">seq\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta.gz"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        cfg_file = temp_dir / "config.json"
        cfg_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(cfg_file))
        assert loaded.ref_genome.coding_type == CodingType.GZIP
        assert loaded.mod_genome.coding_type == CodingType.NONE

    def test_genome_file_no_extension(self,temp_dir):
        """Genome files without an extension are still accepted; paths resolve relative to config."""
        (temp_dir / "ref").write_text(">seq\nATCG\n")  # no suffix
        (temp_dir / "mod").write_text(">seq\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref"},
            "mod_genome_filename": {"filename": "mod"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        cfg_file = temp_dir / "config.json"
        cfg_file.write_text(json.dumps(config))

        with pytest.raises((ConfigurationError)):
            ConfigManager.load(str(cfg_file))

    def test_genome_file_fasta_gz_coding(self,temp_dir):
        """FASTA genomes plus explicit gz output coding on both ref and mod."""
        (temp_dir / "ref.fasta.gzip").write_text(">seq\nATCG\n")
        (temp_dir / "mod.fasta.gzip").write_text(">seq\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta.gzip"},
            "mod_genome_filename": {"filename": "mod.fasta.gzip"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        cfg_file = temp_dir / "config.json"
        cfg_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(cfg_file))
        assert loaded.ref_genome.detected_format == GenomeFormat.FASTA
        assert loaded.ref_genome.coding_type == CodingType.GZIP
        assert loaded.mod_genome.coding_type == CodingType.GZIP

    def test_read_gz_coding(self,temp_dir):
        """Reads entry with gz output coding is preserved in the model."""
        # Minimal input genomes (required)
        (temp_dir / "ref.fasta").write_text(">s\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">s\nATCG\n")
        (temp_dir / "reads.fastq.gz").write_text("@r1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [
                {"filename": "reads.fastq.gz", "ngs_type": "illumina"}
            ],
        }
        cfg_file = temp_dir / "config.json"
        cfg_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(cfg_file))
        assert len(loaded.reads) == 1
        assert loaded.reads[0].filepath.is_absolute()
        assert loaded.reads[0].ngs_type == "illumina"
        assert loaded.reads[0].coding_type == CodingType.GZIP


    def test_feature_gz_coding(self,temp_dir):
        """
        Feature file with gz output coding is preserved.
        Assumes features are given as a list under the 'features' key.
        """
        (temp_dir / "ref.fasta").write_text(">s\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">s\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")
        (temp_dir / "features.gff.gz").write_text("##gff-version 3\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
            "ref_feature_filename": {"filename": "features.gff.gz"}
        }
        cfg_file = temp_dir / "config.json"
        cfg_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(cfg_file))
        # If your schema uses a different key than 'features', rename here.
        assert loaded.ref_feature.filepath.is_absolute()
        assert loaded.ref_feature.detected_format == FeatureFormat.GFF
        assert loaded.ref_feature.coding_type == CodingType.GZIP

class TestConfigManagerUtilities:
    """Test utility methods in ConfigManager."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_detect_compression_type_none(self):
        """Test detecting no compression."""
        filepath = Path("genome.fasta")
        coding_type = ConfigManager._detect_compression_type(filepath)
        assert coding_type == CodingType.NONE

    def test_detect_compression_type_gzip(self):
        """Test detecting gzip compression."""
        filepath = Path("genome.fasta.gz")
        coding_type = ConfigManager._detect_compression_type(filepath)
        assert coding_type == CodingType.GZIP

    def test_detect_compression_type_gzip_alt(self):
        """Test detecting gzip with .gzip extension."""
        filepath = Path("genome.fasta.gzip")
        coding_type = ConfigManager._detect_compression_type(filepath)
        assert coding_type == CodingType.GZIP

    def test_detect_compression_type_bzip2(self):
        """Test detecting bzip2 compression."""
        filepath = Path("genome.fasta.bz2")
        coding_type = ConfigManager._detect_compression_type(filepath)
        assert coding_type == CodingType.BZIP2

    def test_detect_compression_type_bzip2_alt(self):
        """Test detecting bzip2 with .bzip2 extension."""
        filepath = Path("genome.fasta.bzip2")
        coding_type = ConfigManager._detect_compression_type(filepath)
        assert coding_type == CodingType.BZIP2

    def test_detect_file_format_fasta(self):
        """Test detecting FASTA format."""
        filepath = Path("genome.fasta")
        file_format = ConfigManager._detect_file_format(filepath, GenomeFormat)
        assert file_format == GenomeFormat.FASTA

    def test_detect_file_format_fasta_compressed(self):
        """Test detecting FASTA format from compressed file."""
        filepath = Path("genome.fasta.gz")
        file_format = ConfigManager._detect_file_format(filepath, GenomeFormat)
        assert file_format == GenomeFormat.FASTA

    def test_detect_file_format_genbank(self):
        """Test detecting GenBank format."""
        filepath = Path("genome.gbk")
        file_format = ConfigManager._detect_file_format(filepath, GenomeFormat)
        assert file_format == GenomeFormat.GENBANK

    def test_detect_file_format_fastq(self):
        """Test detecting FASTQ format."""
        filepath = Path("reads.fastq.gz")
        file_format = ConfigManager._detect_file_format(filepath, ReadFormat)
        assert file_format == ReadFormat.FASTQ

    def test_detect_file_format_gff(self):
        """Test detecting GFF format."""
        filepath = Path("features.gff")
        file_format = ConfigManager._detect_file_format(filepath, FeatureFormat)
        assert file_format == FeatureFormat.GFF

    def test_detect_file_format_no_extension(self):
        """Test error when no extension found."""
        filepath = Path("genome")
        with pytest.raises(ValueError, match="no extension found"):
            ConfigManager._detect_file_format(filepath, GenomeFormat)

    def test_parse_config_file_value_dict(self):
        """Test parsing config value as dict."""
        value = {"filename": "genome.fasta", "extra_key": "value"}
        filename, extra = ConfigManager._parse_config_file_value(value, "ref_genome")
        assert filename == "genome.fasta"
        assert extra == {"extra_key": "value"}

    def test_parse_config_file_value_string(self):
        """Test parsing config value as string."""
        value = "genome.fasta"
        filename, extra = ConfigManager._parse_config_file_value(value, "ref_genome")
        assert filename == "genome.fasta"
        assert extra == {}

    def test_parse_config_file_value_missing_filename(self):
        """Test error when dict missing filename field."""
        value = {"other_key": "value"}
        with pytest.raises(ValueError, match="must contain 'filename' field"):
            ConfigManager._parse_config_file_value(value, "ref_genome")

    def test_parse_config_file_value_invalid_type(self):
        """Test error with invalid value type."""
        value = 123  # Not a dict or string
        with pytest.raises(ValueError, match="must be a dict or string"):
            ConfigManager._parse_config_file_value(value, "ref_genome")


class TestConfigManagerOutputDirectory:
    """Test output directory setup functionality."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_output_directory_created(self, temp_dir):
        """Test that output directory is created during config load."""
        # Create minimal config files
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        loaded_config = ConfigManager.load(str(config_file))

        # Check output_dir was set
        assert loaded_config.output_dir is not None
        assert loaded_config.output_dir.exists()
        assert loaded_config.output_dir.is_dir()

        # Check it's at the expected location (config_dir.parent / "valid")
        expected_output_dir = config_file.parent.parent / "valid"
        assert loaded_config.output_dir == expected_output_dir

    def test_output_directory_path_structure(self, temp_dir):
        """Test output directory is at config_dir.parent / 'valid'."""
        # Create config in a subdirectory
        config_subdir = temp_dir / "config"
        config_subdir.mkdir()

        (config_subdir / "ref.fasta").write_text(">seq1\nATCG\n")
        (config_subdir / "mod.fasta").write_text(">seq1\nATCG\n")
        (config_subdir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = config_subdir / "config.json"
        config_file.write_text(json.dumps(config))

        loaded_config = ConfigManager.load(str(config_file))

        # Output dir should be at temp_dir / "valid" (parent of config_subdir)
        expected_output_dir = temp_dir / "valid"
        assert loaded_config.output_dir == expected_output_dir
        assert expected_output_dir.exists()

    def test_genome_configs_have_output_dir(self, temp_dir):
        """Test that genome configs have output_dir set."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        loaded_config = ConfigManager.load(str(config_file))

        # Check genome configs have output_dir
        expected_output_dir = temp_dir.parent / "valid"
        assert loaded_config.ref_genome.output_dir == expected_output_dir
        assert loaded_config.mod_genome.output_dir == expected_output_dir

    def test_read_configs_have_output_dir(self, temp_dir):
        """Test that read configs have output_dir set."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        loaded_config = ConfigManager.load(str(config_file))

        # Check read configs have output_dir
        expected_output_dir = temp_dir.parent / "valid"
        assert all(read.output_dir == expected_output_dir for read in loaded_config.reads)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

