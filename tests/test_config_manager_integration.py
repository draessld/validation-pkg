"""
Integration tests for ConfigManager using real test fixtures.

These tests use realistic configuration files in the fixtures/ directory to test
the complete configuration loading workflow end-to-end.
"""

import pytest
import json
from pathlib import Path

from validation_pkg.config_manager import ConfigManager
from validation_pkg.exceptions import ConfigurationError, FileNotFoundError as ValidationFileNotFoundError
from validation_pkg.utils.formats import GenomeFormat, ReadFormat, FeatureFormat, CodingType


class TestConfigManagerIntegration:
    """Integration tests using realistic test fixtures."""

    @pytest.fixture
    def fixtures_dir(self):
        """Get path to fixtures directory."""
        return Path(__file__).parent / "fixtures"

    def test_config_01_minimal_valid(self, fixtures_dir):
        """Test Case 1: Valid minimal configuration."""
        test_case = fixtures_dir / "config_test_01_minimal_valid" / "data"
        config_file = test_case / "config.json"

        config = ConfigManager.load(str(config_file))

        # Verify required fields
        assert config.ref_genome is not None
        assert config.mod_genome is not None
        assert len(config.reads) == 1

        # Verify paths are absolute
        assert config.ref_genome.filepath.is_absolute()
        assert config.mod_genome.filepath.is_absolute()
        assert config.reads[0].filepath.is_absolute()

        # Verify files exist
        assert config.ref_genome.filepath.exists()
        assert config.mod_genome.filepath.exists()
        assert config.reads[0].filepath.exists()

        # Verify formats detected
        assert config.ref_genome.detected_format == GenomeFormat.FASTA
        assert config.mod_genome.detected_format == GenomeFormat.FASTA
        assert config.reads[0].detected_format == ReadFormat.FASTQ

        # Verify compression types
        assert config.ref_genome.coding_type == CodingType.NONE
        assert config.mod_genome.coding_type == CodingType.NONE
        assert config.reads[0].coding_type == CodingType.NONE

        # Verify read type
        assert config.reads[0].ngs_type == "illumina"

        # Verify optional fields are None
        assert config.ref_plasmid is None
        assert config.mod_plasmid is None
        assert config.ref_feature is None
        assert config.mod_feature is None

    def test_config_02_full_config(self, fixtures_dir):
        """Test Case 2: Full configuration with all optional fields."""
        test_case = fixtures_dir / "config_test_02_full_config" / "data"
        config_file = test_case / "config.json"

        config = ConfigManager.load(str(config_file))

        # Verify all fields are populated
        assert config.ref_genome is not None
        assert config.mod_genome is not None
        assert config.ref_plasmid is not None
        assert config.mod_plasmid is not None
        assert config.ref_feature is not None
        assert config.mod_feature is not None
        assert len(config.reads) == 2

        # Verify formats
        assert config.ref_genome.detected_format == GenomeFormat.FASTA
        assert config.mod_genome.detected_format == GenomeFormat.FASTA
        assert config.ref_plasmid.detected_format == GenomeFormat.GENBANK
        assert config.mod_plasmid.detected_format == GenomeFormat.GENBANK
        assert config.ref_feature.detected_format == FeatureFormat.GFF
        assert config.mod_feature.detected_format == FeatureFormat.GFF

        # Verify options
        assert config.options["threads"] == 8
        assert config.options["verbose"] is True

        # Verify all files exist
        assert config.ref_genome.filepath.exists()
        assert config.mod_genome.filepath.exists()
        assert config.ref_plasmid.filepath.exists()
        assert config.mod_plasmid.filepath.exists()
        assert config.ref_feature.filepath.exists()
        assert config.mod_feature.filepath.exists()
        assert all(read.filepath.exists() for read in config.reads)

    def test_config_03_directory_reads(self, fixtures_dir):
        """Test Case 3: Directory-based reads."""
        test_case = fixtures_dir / "config_test_03_directory_reads" / "data"
        config_file = test_case / "config.json"

        config = ConfigManager.load(str(config_file))

        # Verify reads were loaded from directory
        assert len(config.reads) == 3

        # Verify all reads have correct NGS type
        assert all(read.ngs_type == "ont" for read in config.reads)

        # Verify all read files exist
        assert all(read.filepath.exists() for read in config.reads)

        # Verify all reads are from the directory
        reads_dir = test_case / "reads_dir"
        for read in config.reads:
            assert read.filepath.parent == reads_dir

    def test_config_04_compressed_files(self, fixtures_dir):
        """Test Case 4: Compressed files."""
        test_case = fixtures_dir / "config_test_04_compressed" / "data"
        config_file = test_case / "config.json"

        config = ConfigManager.load(str(config_file))

        # Verify compression types detected
        assert config.ref_genome.coding_type == CodingType.GZIP
        assert config.mod_genome.coding_type == CodingType.BZIP2
        assert config.reads[0].coding_type == CodingType.GZIP

        # Verify file extensions
        assert config.ref_genome.filepath.suffix == ".gz"
        assert config.mod_genome.filepath.suffix == ".bz2"
        assert config.reads[0].filepath.suffix == ".gz"

        # Verify files exist
        assert config.ref_genome.filepath.exists()
        assert config.mod_genome.filepath.exists()
        assert config.reads[0].filepath.exists()

    def test_config_05_multiple_ngs_types(self, fixtures_dir):
        """Test Case 5: Multiple NGS types."""
        test_case = fixtures_dir / "config_test_05_multiple_ngs_types" / "data"
        config_file = test_case / "config.json"

        config = ConfigManager.load(str(config_file))

        # Verify we have 3 read files
        assert len(config.reads) == 3

        # Verify each NGS type is present
        ngs_types = [read.ngs_type for read in config.reads]
        assert "illumina" in ngs_types
        assert "ont" in ngs_types
        assert "pacbio" in ngs_types

        # Verify all files exist
        assert all(read.filepath.exists() for read in config.reads)

    def test_config_06_missing_config_file(self, fixtures_dir):
        """Test Case 6: Missing config file (should fail)."""
        test_case = fixtures_dir / "config_test_06_missing_config" / "data"
        config_file = test_case / "config.json"

        # Should raise FileNotFoundError
        with pytest.raises(ValidationFileNotFoundError, match="Configuration file not found"):
            ConfigManager.load(str(config_file))

    def test_config_07_missing_required_field(self, fixtures_dir):
        """Test Case 7: Missing required field (should fail)."""
        test_case = fixtures_dir / "config_test_07_missing_required_field" / "data"
        config_file = test_case / "config.json"

        # Should raise ConfigurationError about missing 'reads'
        with pytest.raises((ConfigurationError, ValueError), match="reads"):
            ConfigManager.load(str(config_file))

    def test_config_08_missing_file(self, fixtures_dir):
        """Test Case 8: Missing genome file (should fail)."""
        test_case = fixtures_dir / "config_test_08_missing_file" / "data"
        config_file = test_case / "config.json"

        # Should raise FileNotFoundError about missing ref_genome.fasta
        with pytest.raises(ValidationFileNotFoundError, match="ref_genome|ref_genome.fasta"):
            ConfigManager.load(str(config_file))

    def test_config_09_invalid_json(self, fixtures_dir):
        """Test Case 9: Invalid JSON (should fail)."""
        test_case = fixtures_dir / "config_test_09_invalid_json" / "data"
        config_file = test_case / "config.json"

        # Should raise JSONDecodeError or ConfigurationError
        with pytest.raises((ConfigurationError, json.JSONDecodeError)):
            ConfigManager.load(str(config_file))

    def test_config_10_invalid_ngs_type(self, fixtures_dir):
        """Test Case 10: Invalid NGS type (should fail)."""
        test_case = fixtures_dir / "config_test_10_invalid_ngs_type" / "data"
        config_file = test_case / "config.json"

        # Should raise ValueError about invalid ngs_type
        with pytest.raises((ConfigurationError, ValueError), match="Invalid ngs_type|454"):
            ConfigManager.load(str(config_file))

    def test_config_output_directory_created(self, fixtures_dir):
        """Test that output directory is created during config load."""
        test_case = fixtures_dir / "config_test_01_minimal_valid" / "data"
        config_file = test_case / "config.json"

        config = ConfigManager.load(str(config_file))

        # Verify output_dir was created
        assert config.output_dir is not None
        assert config.output_dir.exists()
        assert config.output_dir.is_dir()

        # Verify it's at expected location (config_dir.parent / "valid")
        expected_output_dir = config_file.parent.parent / "valid"
        assert config.output_dir == expected_output_dir

    def test_config_paths_relative_to_config_dir(self, fixtures_dir):
        """Test that paths are resolved relative to config directory."""
        test_case = fixtures_dir / "config_test_01_minimal_valid" / "data"
        config_file = test_case / "config.json"

        config = ConfigManager.load(str(config_file))

        # All file paths should be in the same directory as config
        config_dir = config_file.parent
        assert config.ref_genome.filepath.parent == config_dir
        assert config.mod_genome.filepath.parent == config_dir
        assert config.reads[0].filepath.parent == config_dir

    def test_config_formats_detected_correctly(self, fixtures_dir):
        """Test that file formats are detected correctly from extensions."""
        test_case = fixtures_dir / "config_test_02_full_config" / "data"
        config_file = test_case / "config.json"

        config = ConfigManager.load(str(config_file))

        # Genome formats
        assert config.ref_genome.detected_format == GenomeFormat.FASTA
        assert config.ref_plasmid.detected_format == GenomeFormat.GENBANK

        # Read format
        assert all(read.detected_format == ReadFormat.FASTQ for read in config.reads)

        # Feature format
        assert config.ref_feature.detected_format == FeatureFormat.GFF
        assert config.mod_feature.detected_format == FeatureFormat.GFF

    def test_config_filenames_preserved(self, fixtures_dir):
        """Test that original filenames are preserved in config."""
        test_case = fixtures_dir / "config_test_01_minimal_valid" / "data"
        config_file = test_case / "config.json"

        config = ConfigManager.load(str(config_file))

        # Verify filenames match what's in the JSON
        assert config.ref_genome.filename == "ref_genome.fasta"
        assert config.mod_genome.filename == "mod_genome.fasta"
        assert config.reads[0].filename == "reads.fastq"

    def test_config_extra_fields_ignored(self, fixtures_dir):
        """Test that extra fields in config are stored but don't cause errors."""
        # This test would require a fixture with extra fields
        # For now, we verify the behavior with the existing fixtures
        test_case = fixtures_dir / "config_test_02_full_config" / "data"
        config_file = test_case / "config.json"

        # Should load successfully even with extra fields like options
        config = ConfigManager.load(str(config_file))
        assert config is not None

    def test_config_read_configs_have_output_dir(self, fixtures_dir):
        """Test that read configs have output_dir set."""
        test_case = fixtures_dir / "config_test_01_minimal_valid" / "data"
        config_file = test_case / "config.json"

        config = ConfigManager.load(str(config_file))

        # All reads should have output_dir set
        expected_output_dir = config_file.parent.parent / "valid"
        assert all(read.output_dir == expected_output_dir for read in config.reads)

    def test_config_genome_configs_have_output_dir(self, fixtures_dir):
        """Test that genome configs have output_dir set."""
        test_case = fixtures_dir / "config_test_01_minimal_valid" / "data"
        config_file = test_case / "config.json"

        config = ConfigManager.load(str(config_file))

        # Genomes should have output_dir set
        expected_output_dir = config_file.parent.parent / "valid"
        assert config.ref_genome.output_dir == expected_output_dir
        assert config.mod_genome.output_dir == expected_output_dir


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
