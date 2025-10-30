"""
Tests for 4-layer settings priority system with global options.

This module tests that the 4-layer settings precedence system works correctly:
Layer 1 (lowest): Defaults
Layer 2: Global options (validation_level and threads only)
Layer 3: File-level settings (all Settings fields)
Layer 4 (highest): User settings (complete replacement)
"""

import pytest
from pathlib import Path
from validation_pkg.validators.genome_validator import GenomeValidator
from validation_pkg.validators.read_validator import ReadValidator
from validation_pkg.validators.feature_validator import FeatureValidator
from validation_pkg.config_manager import GenomeConfig, ReadConfig, FeatureConfig, ConfigManager
from validation_pkg.utils.formats import GenomeFormat, ReadFormat, FeatureFormat, CodingType
from validation_pkg.exceptions import ConfigurationError


class TestGlobalOptionsValidation:
    """Tests for global options validation in ConfigManager."""

    def _create_minimal_config_files(self, tmp_path):
        """Helper to create minimal required files for config"""
        ref_file = tmp_path / "ref.fasta"
        ref_file.write_text(">seq1\nACGT\n")
        mod_file = tmp_path / "mod.fasta"
        mod_file.write_text(">seq1\nACGT\n")
        read_file = tmp_path / "read.fastq"
        read_file.write_text("@read1\nACGT\n+\nIIII\n")
        return str(ref_file), str(mod_file), str(read_file)

    def test_valid_global_options_threads_only(self, tmp_path):
        """Only threads in global options should be valid"""
        ref_file, mod_file, read_file = self._create_minimal_config_files(tmp_path)

        config_data = {
            "ref_genome_filename": {"filename": ref_file},
            "mod_genome_filename": {"filename": mod_file},
            "reads": [{"filename": read_file, "ngs_type": "illumina"}],
            "options": {"threads": 4}
        }
        config_file = tmp_path / "config.json"
        import json
        config_file.write_text(json.dumps(config_data))

        config = ConfigManager.load(str(config_file))
        assert config.options == {"threads": 4}

    def test_valid_global_options_validation_level_only(self, tmp_path):
        """Only validation_level in global options should be valid"""
        ref_file, mod_file, read_file = self._create_minimal_config_files(tmp_path)

        config_data = {
            "ref_genome_filename": {"filename": ref_file},
            "mod_genome_filename": {"filename": mod_file},
            "reads": [{"filename": read_file, "ngs_type": "illumina"}],
            "options": {"validation_level": "trust"}
        }
        config_file = tmp_path / "config.json"
        import json
        config_file.write_text(json.dumps(config_data))

        config = ConfigManager.load(str(config_file))
        assert config.options == {"validation_level": "trust"}

    def test_valid_global_options_both(self, tmp_path):
        """Both threads and validation_level together should be valid"""
        ref_file, mod_file, read_file = self._create_minimal_config_files(tmp_path)

        config_data = {
            "ref_genome_filename": {"filename": ref_file},
            "mod_genome_filename": {"filename": mod_file},
            "reads": [{"filename": read_file, "ngs_type": "illumina"}],
            "options": {"threads": 4, "validation_level": "minimal"}
        }
        config_file = tmp_path / "config.json"
        import json
        config_file.write_text(json.dumps(config_data))

        config = ConfigManager.load(str(config_file))
        assert config.options == {"threads": 4, "validation_level": "minimal"}

    def test_invalid_global_option_field(self, tmp_path):
        """Invalid field names in global options should raise ConfigurationError"""
        ref_file, mod_file, read_file = self._create_minimal_config_files(tmp_path)

        config_data = {
            "ref_genome_filename": {"filename": ref_file},
            "mod_genome_filename": {"filename": mod_file},
            "reads": [{"filename": read_file, "ngs_type": "illumina"}],
            "options": {"threads": 4, "plasmid_split": True}  # plasmid_split not allowed
        }
        config_file = tmp_path / "config.json"
        import json
        config_file.write_text(json.dumps(config_data))

        with pytest.raises(ConfigurationError, match="Invalid global options.*plasmid_split"):
            ConfigManager.load(str(config_file))

    def test_invalid_validation_level_value(self, tmp_path):
        """Invalid validation_level value should raise ConfigurationError"""
        ref_file, mod_file, read_file = self._create_minimal_config_files(tmp_path)

        config_data = {
            "ref_genome_filename": {"filename": ref_file},
            "mod_genome_filename": {"filename": mod_file},
            "reads": [{"filename": read_file, "ngs_type": "illumina"}],
            "options": {"validation_level": "invalid_level"}
        }
        config_file = tmp_path / "config.json"
        import json
        config_file.write_text(json.dumps(config_data))

        with pytest.raises(ConfigurationError, match="Invalid validation_level"):
            ConfigManager.load(str(config_file))

    def test_invalid_threads_value_negative(self, tmp_path):
        """Negative threads value should raise ConfigurationError"""
        ref_file, mod_file, read_file = self._create_minimal_config_files(tmp_path)

        config_data = {
            "ref_genome_filename": {"filename": ref_file},
            "mod_genome_filename": {"filename": mod_file},
            "reads": [{"filename": read_file, "ngs_type": "illumina"}],
            "options": {"threads": -1}
        }
        config_file = tmp_path / "config.json"
        import json
        config_file.write_text(json.dumps(config_data))

        with pytest.raises(ConfigurationError, match="must be a positive integer"):
            ConfigManager.load(str(config_file))

    def test_invalid_custom_field_raises_error(self, tmp_path):
        """Custom fields in global options should raise ConfigurationError"""
        ref_file, mod_file, read_file = self._create_minimal_config_files(tmp_path)

        config_data = {
            "ref_genome_filename": {"filename": ref_file},
            "mod_genome_filename": {"filename": mod_file},
            "reads": [{"filename": read_file, "ngs_type": "illumina"}],
            "options": {
                "threads": 4,
                "custom_field": "not_allowed"  # Should raise error
            }
        }
        config_file = tmp_path / "config.json"
        import json
        config_file.write_text(json.dumps(config_data))

        # Should RAISE because custom_field is not in ALLOWED_GLOBAL_OPTIONS
        with pytest.raises(ConfigurationError, match="Invalid global options.*custom_field"):
            ConfigManager.load(str(config_file))


class TestGenomeValidator4LayerPrecedence:
    """Tests for 4-layer settings precedence in GenomeValidator."""

    def test_layer1_defaults_only(self, tmp_path):
        """Layer 1: When nothing is specified, use defaults"""
        config = GenomeConfig(
            filename="genome.fasta",
            filepath=Path(__file__).parent / "fixtures" / "sample.fasta",
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA,
            settings_dict={},
            global_options={}
        )

        validator = GenomeValidator(config, tmp_path, settings=None)

        # Should have default values
        assert validator.settings.validation_level == 'strict'  # Default
        assert validator.settings.threads is None  # Default
        assert validator.settings.min_sequence_length == 100  # Default

    def test_layer2_global_options_override_defaults(self, tmp_path):
        """Layer 2: Global options override defaults"""
        config = GenomeConfig(
            filename="genome.fasta",
            filepath=Path(__file__).parent / "fixtures" / "sample.fasta",
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA,
            settings_dict={},
            global_options={"threads": 6, "validation_level": "trust"}
        )

        validator = GenomeValidator(config, tmp_path, settings=None)

        # Global options should override defaults
        assert validator.settings.validation_level == 'trust'  # From global
        assert validator.settings.threads == 6  # From global
        assert validator.settings.min_sequence_length == 100  # Still default

    def test_layer3_file_settings_override_global(self, tmp_path):
        """Layer 3: File-level settings override global options"""
        config = GenomeConfig(
            filename="genome.fasta",
            filepath=Path(__file__).parent / "fixtures" / "sample.fasta",
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA,
            settings_dict={
                "validation_level": "minimal",  # Override global
                "min_sequence_length": 500  # New setting
            },
            global_options={"threads": 6, "validation_level": "trust"}
        )

        validator = GenomeValidator(config, tmp_path, settings=None)

        # File-level should override global
        assert validator.settings.validation_level == 'minimal'  # File overrides global
        assert validator.settings.threads == 6  # From global (not overridden)
        assert validator.settings.min_sequence_length == 500  # From file-level

    def test_layer4_user_settings_override_all(self, tmp_path):
        """Layer 4: User settings have highest priority"""
        config = GenomeConfig(
            filename="genome.fasta",
            filepath=Path(__file__).parent / "fixtures" / "sample.fasta",
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA,
            settings_dict={
                "validation_level": "minimal",
                "min_sequence_length": 500
            },
            global_options={"threads": 6, "validation_level": "trust"}
        )

        # User provides settings - should replace everything
        user_settings = GenomeValidator.Settings(
            validation_level='strict',
            threads=2
        )
        validator = GenomeValidator(config, tmp_path, user_settings)

        # User settings have highest priority
        assert validator.settings.validation_level == 'strict'  # User value
        assert validator.settings.threads == 2  # User value
        assert validator.settings.min_sequence_length == 100  # User default (NOT from file/global)

    def test_file_level_overrides_global_warning(self, tmp_path):
        """File-level overriding global should log warning"""
        from validation_pkg.logger import get_logger
        import logging

        # Get logger and set it to capture warnings
        logger = get_logger()

        config = GenomeConfig(
            filename="genome.fasta",
            filepath=Path(__file__).parent / "fixtures" / "sample.fasta",
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA,
            settings_dict={"validation_level": "minimal"},  # Override global
            global_options={"validation_level": "trust"}
        )

        # Create validator - this should trigger the warning
        validator = GenomeValidator(config, tmp_path, settings=None)

        # Check that validation_level was correctly applied (file-level wins)
        assert validator.settings.validation_level == "minimal"

    def test_user_settings_override_warning(self, tmp_path):
        """User settings overriding config should log warning with source"""
        from validation_pkg.logger import get_logger

        logger = get_logger()

        config = GenomeConfig(
            filename="genome.fasta",
            filepath=Path(__file__).parent / "fixtures" / "sample.fasta",
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA,
            settings_dict={"validation_level": "minimal"},
            global_options={"threads": 6}
        )

        user_settings = GenomeValidator.Settings(validation_level='strict')
        validator = GenomeValidator(config, tmp_path, user_settings)

        # Check that user settings won (user settings REPLACE, not merge)
        assert validator.settings.validation_level == 'strict'
        # threads is None because user settings completely replace (they didn't specify threads)
        assert validator.settings.threads is None


class TestReadValidator4LayerPrecedence:
    """Tests for 4-layer settings precedence in ReadValidator."""

    def test_layer2_global_options_override_defaults(self, tmp_path):
        """Layer 2: Global options override defaults"""
        config = ReadConfig(
            filename="read.fastq",
            filepath=Path(__file__).parent / "fixtures" / "sample_illumina_1.fastq",
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            settings_dict={},
            global_options={"threads": 4, "validation_level": "minimal"}
        )

        validator = ReadValidator(config, tmp_path, settings=None)

        assert validator.settings.validation_level == 'minimal'  # From global
        assert validator.settings.threads == 4  # From global

    def test_layer3_file_settings_override_global(self, tmp_path):
        """Layer 3: File-level settings override global options"""
        config = ReadConfig(
            filename="read.fastq",
            filepath=Path(__file__).parent / "fixtures" / "sample_illumina_1.fastq",
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            settings_dict={
                "validation_level": "strict",  # Override global
                "check_invalid_chars": False
            },
            global_options={"validation_level": "minimal"}
        )

        validator = ReadValidator(config, tmp_path, settings=None)

        assert validator.settings.validation_level == 'strict'  # File overrides global
        assert validator.settings.check_invalid_chars is False  # From file-level

    def test_layer4_user_settings_override_all(self, tmp_path):
        """Layer 4: User settings have highest priority"""
        config = ReadConfig(
            filename="read.fastq",
            filepath=Path(__file__).parent / "fixtures" / "sample_illumina_1.fastq",
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            settings_dict={"validation_level": "strict"},
            global_options={"threads": 4}
        )

        user_settings = ReadValidator.Settings(validation_level='minimal', threads=2)
        validator = ReadValidator(config, tmp_path, user_settings)

        assert validator.settings.validation_level == 'minimal'  # User value
        assert validator.settings.threads == 2  # User value


class TestFeatureValidator4LayerPrecedence:
    """Tests for 4-layer settings precedence in FeatureValidator."""

    def test_layer2_global_options_override_defaults(self, tmp_path):
        """Layer 2: Global options override defaults"""
        config = FeatureConfig(
            filename="features.gff",
            filepath=Path(__file__).parent / "fixtures" / "sample.gff",
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            settings_dict={},
            global_options={"threads": 8, "validation_level": "trust"}
        )

        validator = FeatureValidator(config, tmp_path, settings=None)

        assert validator.settings.validation_level == 'trust'  # From global
        assert validator.settings.threads == 8  # From global

    def test_layer3_file_settings_override_global(self, tmp_path):
        """Layer 3: File-level settings override global options"""
        config = FeatureConfig(
            filename="features.gff",
            filepath=Path(__file__).parent / "fixtures" / "sample.gff",
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            settings_dict={
                "validation_level": "minimal",  # Override global
                "sort_by_position": True
            },
            global_options={"validation_level": "trust", "threads": 8}
        )

        validator = FeatureValidator(config, tmp_path, settings=None)

        assert validator.settings.validation_level == 'minimal'  # File overrides global
        assert validator.settings.threads == 8  # From global (not overridden)
        assert validator.settings.sort_by_position is True  # From file-level

    def test_layer4_user_settings_override_all(self, tmp_path):
        """Layer 4: User settings have highest priority"""
        config = FeatureConfig(
            filename="features.gff",
            filepath=Path(__file__).parent / "fixtures" / "sample.gff",
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            settings_dict={"sort_by_position": True},
            global_options={"threads": 8}
        )

        user_settings = FeatureValidator.Settings(sort_by_position=False, threads=2)
        validator = FeatureValidator(config, tmp_path, user_settings)

        assert validator.settings.sort_by_position is False  # User value
        assert validator.settings.threads == 2  # User value


class TestGlobalOptionsPropagation:
    """Tests that global_options are properly propagated to all config objects."""

    def test_global_options_propagated_to_genome_config(self, tmp_path):
        """Global options should be copied to GenomeConfig.global_options"""
        ref_file = tmp_path / "ref.fasta"
        ref_file.write_text(">seq1\nACGT\n")
        mod_file = tmp_path / "mod.fasta"
        mod_file.write_text(">seq1\nACGT\n")
        read_file = tmp_path / "read.fastq"
        read_file.write_text("@read1\nACGT\n+\nIIII\n")

        config_data = {
            "ref_genome_filename": {"filename": str(ref_file)},
            "mod_genome_filename": {"filename": str(mod_file)},
            "reads": [{"filename": str(read_file), "ngs_type": "illumina"}],
            "options": {"threads": 4, "validation_level": "trust"}
        }
        config_file = tmp_path / "config.json"
        import json
        config_file.write_text(json.dumps(config_data))

        config = ConfigManager.load(str(config_file))

        assert config.ref_genome is not None
        assert config.ref_genome.global_options == {"threads": 4, "validation_level": "trust"}

    def test_global_options_propagated_to_read_configs(self, tmp_path):
        """Global options should be copied to all ReadConfig.global_options"""
        ref_file = tmp_path / "ref.fasta"
        ref_file.write_text(">seq1\nACGT\n")
        mod_file = tmp_path / "mod.fasta"
        mod_file.write_text(">seq1\nACGT\n")
        read1_file = tmp_path / "read1.fastq"
        read1_file.write_text("@read1\nACGT\n+\nIIII\n")
        read2_file = tmp_path / "read2.fastq"
        read2_file.write_text("@read2\nACGT\n+\nIIII\n")

        config_data = {
            "ref_genome_filename": {"filename": str(ref_file)},
            "mod_genome_filename": {"filename": str(mod_file)},
            "reads": [
                {"filename": str(read1_file), "ngs_type": "illumina"},
                {"filename": str(read2_file), "ngs_type": "ont"}
            ],
            "options": {"threads": 6, "validation_level": "minimal"}
        }
        config_file = tmp_path / "config.json"
        import json
        config_file.write_text(json.dumps(config_data))

        config = ConfigManager.load(str(config_file))

        assert len(config.reads) == 2
        for read_config in config.reads:
            assert read_config.global_options == {"threads": 6, "validation_level": "minimal"}

    def test_global_options_propagated_to_feature_config(self, tmp_path):
        """Global options should be copied to FeatureConfig.global_options"""
        ref_file = tmp_path / "ref.fasta"
        ref_file.write_text(">seq1\nACGT\n")
        mod_file = tmp_path / "mod.fasta"
        mod_file.write_text(">seq1\nACGT\n")
        read_file = tmp_path / "read.fastq"
        read_file.write_text("@read1\nACGT\n+\nIIII\n")
        feature_file = tmp_path / "features.gff"
        feature_file.write_text("##gff-version 3\nseq1\t.\tgene\t1\t4\t.\t+\t.\tID=gene1\n")

        config_data = {
            "ref_genome_filename": {"filename": str(ref_file)},
            "mod_genome_filename": {"filename": str(mod_file)},
            "reads": [{"filename": str(read_file), "ngs_type": "illumina"}],
            "ref_feature_filename": {"filename": str(feature_file)},
            "options": {"threads": 2, "validation_level": "strict"}
        }
        config_file = tmp_path / "config.json"
        import json
        config_file.write_text(json.dumps(config_data))

        config = ConfigManager.load(str(config_file))

        assert config.ref_feature is not None
        assert config.ref_feature.global_options == {"threads": 2, "validation_level": "strict"}


class TestBackwardsCompatibility:
    """Ensure backwards compatibility with code not using global options."""

    def test_config_without_global_options_still_works(self, tmp_path):
        """Configs without global_options field should use empty dict"""
        config = GenomeConfig(
            filename="genome.fasta",
            filepath=Path(__file__).parent / "fixtures" / "sample.fasta",
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
            # No global_options provided
        )

        # Should default to empty dict
        assert config.global_options == {}

        # Validator should work fine
        validator = GenomeValidator(config, tmp_path, settings=None)
        assert validator.settings.validation_level == 'strict'  # Default

    def test_config_with_none_global_options(self, tmp_path):
        """Configs with global_options=None should convert to empty dict"""
        config = GenomeConfig(
            filename="genome.fasta",
            filepath=Path(__file__).parent / "fixtures" / "sample.fasta",
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA,
            global_options=None
        )

        # Should convert None to empty dict in __post_init__
        assert config.global_options == {}
