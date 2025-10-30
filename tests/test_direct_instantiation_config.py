"""
Tests for direct validator instantiation with config settings.

This module tests that validators properly apply settings from config.json
when instantiated directly (without using the functional API).
"""

import pytest
from pathlib import Path
from validation_pkg.validators.genome_validator import GenomeValidator
from validation_pkg.validators.read_validator import ReadValidator
from validation_pkg.validators.feature_validator import FeatureValidator
from validation_pkg.config_manager import GenomeConfig, ReadConfig, FeatureConfig
from validation_pkg.utils.formats import GenomeFormat, ReadFormat, FeatureFormat, CodingType


class TestGenomeValidatorDirectInstantiation:
    """Tests for GenomeValidator direct instantiation with config settings."""

    def test_applies_config_settings_when_no_user_settings(self, tmp_path):
        """Config settings should be applied when settings=None"""
        # Create config with settings_dict
        config = GenomeConfig(
            filename="genome.fasta",
            filepath=Path(__file__).parent / "fixtures" / "sample.fasta",
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA,
            settings_dict={
                'validation_level': 'trust',
                'threads': 4,
                'min_sequence_length': 500
            }
        )

        # Instantiate without settings
        validator = GenomeValidator(config, tmp_path, settings=None)

        # Check settings from config were applied
        assert validator.settings.validation_level == 'trust'
        assert validator.settings.threads == 4
        assert validator.settings.min_sequence_length == 500

    def test_user_settings_replace_config_settings(self, tmp_path):
        """User settings completely replace config settings (no merge)"""
        # Create config with settings_dict
        config = GenomeConfig(
            filename="genome.fasta",
            filepath=Path(__file__).parent / "fixtures" / "sample.fasta",
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA,
            settings_dict={
                'validation_level': 'trust',
                'threads': 4,
                'min_sequence_length': 500
            }
        )

        # Instantiate with user settings
        # When user provides Settings, it REPLACES config settings (not merged)
        user_settings = GenomeValidator.Settings(validation_level='strict')
        validator = GenomeValidator(config, tmp_path, user_settings)

        # User settings completely replace config settings
        assert validator.settings.validation_level == 'strict'  # User value
        assert validator.settings.threads is None  # User default (config NOT merged)
        assert validator.settings.min_sequence_length == 100  # User default (config NOT merged)

    def test_works_without_config_settings(self, tmp_path):
        """Should still work when config has no settings_dict"""
        # Create config without settings_dict
        config = GenomeConfig(
            filename="genome.fasta",
            filepath=Path(__file__).parent / "fixtures" / "sample.fasta",
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA,
            settings_dict={}  # Empty
        )

        # Instantiate without settings
        validator = GenomeValidator(config, tmp_path)

        # Check defaults are used
        assert validator.settings.validation_level == 'strict'  # Default
        assert validator.settings.threads is None  # Default

    def test_use_settings_without_passing_to_validator(self, tmp_path):
        """Best practice: Don't pass settings if you want config settings"""
        # Create config with multiple settings
        config = GenomeConfig(
            filename="genome.fasta",
            filepath=Path(__file__).parent / "fixtures" / "sample.fasta",
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA,
            settings_dict={
                'validation_level': 'trust',
                'threads': 4,
                'min_sequence_length': 500,
                'plasmid_split': False
            }
        )

        # Best practice: Don't pass settings parameter to use config settings
        validator = GenomeValidator(config, tmp_path, settings=None)

        # All config settings are applied
        assert validator.settings.validation_level == 'trust'
        assert validator.settings.threads == 4
        assert validator.settings.min_sequence_length == 500
        assert validator.settings.plasmid_split is False


class TestReadValidatorDirectInstantiation:
    """Tests for ReadValidator direct instantiation with config settings."""

    def test_applies_config_settings_when_no_user_settings(self, tmp_path):
        """Config settings should be applied when settings=None"""
        # Create config with settings_dict
        config = ReadConfig(
            filename="read.fastq",
            filepath=Path(__file__).parent / "fixtures" / "sample_illumina_1.fastq",
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            settings_dict={
                'validation_level': 'minimal',
                'threads': 2,
                'check_invalid_chars': False
            }
        )

        # Instantiate without settings
        validator = ReadValidator(config, tmp_path, settings=None)

        # Check settings from config were applied
        assert validator.settings.validation_level == 'minimal'
        assert validator.settings.threads == 2
        assert validator.settings.check_invalid_chars is False

    def test_user_settings_replace_config_settings(self, tmp_path):
        """User settings completely replace config settings (no merge)"""
        # Create config with settings_dict
        config = ReadConfig(
            filename="read.fastq",
            filepath=Path(__file__).parent / "fixtures" / "sample_illumina_1.fastq",
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            settings_dict={
                'validation_level': 'minimal',
                'threads': 2
            }
        )

        # Instantiate with user settings
        user_settings = ReadValidator.Settings(validation_level='strict')
        validator = ReadValidator(config, tmp_path, user_settings)

        # User settings completely replace config settings
        assert validator.settings.validation_level == 'strict'  # User value
        assert validator.settings.threads is None  # User default (config NOT merged)


class TestFeatureValidatorDirectInstantiation:
    """Tests for FeatureValidator direct instantiation with config settings."""

    def test_applies_config_settings_when_no_user_settings(self, tmp_path):
        """Config settings should be applied when settings=None"""
        # Create config with settings_dict
        config = FeatureConfig(
            filename="features.gff",
            filepath=Path(__file__).parent / "fixtures" / "sample.gff",
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            settings_dict={
                'validation_level': 'trust',
                'threads': 6,
                'sort_by_position': True
            }
        )

        # Instantiate without settings
        validator = FeatureValidator(config, tmp_path, settings=None)

        # Check settings from config were applied
        assert validator.settings.validation_level == 'trust'
        assert validator.settings.threads == 6
        assert validator.settings.sort_by_position is True

    def test_user_settings_replace_config_settings(self, tmp_path):
        """User settings completely replace config settings (no merge)"""
        # Create config with settings_dict
        config = FeatureConfig(
            filename="features.gff",
            filepath=Path(__file__).parent / "fixtures" / "sample.gff",
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF,
            settings_dict={
                'validation_level': 'trust',
                'sort_by_position': True
            }
        )

        # Instantiate with user settings
        user_settings = FeatureValidator.Settings(sort_by_position=False)
        validator = FeatureValidator(config, tmp_path, user_settings)

        # User settings completely replace config settings
        assert validator.settings.validation_level == 'strict'  # User default (config NOT merged)
        assert validator.settings.sort_by_position is False  # User value


class TestBackwardsCompatibility:
    """Ensure backwards compatibility with existing code."""

    def test_genome_validator_without_config_dict(self, tmp_path):
        """Old code without settings_dict should still work"""
        config = GenomeConfig(
            filename="genome.fasta",
            filepath=Path(__file__).parent / "fixtures" / "sample.fasta",
            coding_type=CodingType.NONE,
            detected_format=GenomeFormat.FASTA
            # No settings_dict provided (None by default)
        )

        # Should work with defaults
        validator = GenomeValidator(config, tmp_path)
        assert validator.settings.validation_level == 'strict'

    def test_read_validator_without_config_dict(self, tmp_path):
        """Old code without settings_dict should still work"""
        config = ReadConfig(
            filename="read.fastq",
            filepath=Path(__file__).parent / "fixtures" / "sample_illumina_1.fastq",
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ
            # No settings_dict provided (None by default)
        )

        # Should work with defaults
        validator = ReadValidator(config, tmp_path)
        assert validator.settings.validation_level == 'strict'

    def test_feature_validator_without_config_dict(self, tmp_path):
        """Old code without settings_dict should still work"""
        config = FeatureConfig(
            filename="features.gff",
            filepath=Path(__file__).parent / "fixtures" / "sample.gff",
            coding_type=CodingType.NONE,
            detected_format=FeatureFormat.GFF
            # No settings_dict provided (None by default)
        )

        # Should work with defaults
        validator = FeatureValidator(config, tmp_path)
        assert validator.settings.validation_level == 'strict'
