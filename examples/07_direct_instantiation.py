#!/usr/bin/env python3
"""
Example 7: Direct Validator Instantiation

This example demonstrates how to create validators directly without
using the functional API. This gives you maximum control over the
validation process.

This is useful when:
- You need to build configs programmatically
- You want to integrate with existing pipelines
- You need to validate files not specified in a config file
- You want fine-grained control over validator instantiation

Key Concepts:
1. Config Settings: Validators automatically apply settings from config.json
2. User Settings: When you pass a Settings object, it REPLACES config settings
3. Best Practice: Don't pass settings=... if you want to use config settings
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from validation_pkg.validators.genome_validator import GenomeValidator
from validation_pkg.validators.read_validator import ReadValidator
from validation_pkg.validators.feature_validator import FeatureValidator
from validation_pkg.config_manager import GenomeConfig, ReadConfig, FeatureConfig
from validation_pkg.utils.formats import GenomeFormat, ReadFormat, FeatureFormat, CodingType

def main():
    # Base paths
    data_dir = Path(__file__).parent / "data"
    output_dir = Path(__file__).parent / "output" / "direct_instantiation"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Direct Validator Instantiation Example")
    print("=" * 70)
    print()

    # Example 1: Genome Validator
    print("1. Creating and using GenomeValidator directly")
    print("-" * 70)

    # Create genome config manually
    genome_config = GenomeConfig(
        filename="genome.fasta",
        filepath=data_dir / "genome.fasta",
        coding_type=CodingType.NONE,  # Uncompressed
        detected_format=GenomeFormat.FASTA,
        extra={}
    )

    # Create settings
    genome_settings = GenomeValidator.Settings()
    genome_settings = genome_settings.update(
        validation_level='strict',
        min_sequence_length=100,
        replace_id_with='seq',
        coding_type='gz'
    )

    # Create validator and validate
    genome_validator = GenomeValidator(genome_config, output_dir, genome_settings)
    genome_validator.validate()

    print(f"Validated: {genome_config.filename}")
    print(f"Settings: {genome_settings}")
    print()

    # Example 2: Read Validator
    print("2. Creating and using ReadValidator directly")
    print("-" * 70)

    # Create read config manually
    read_config = ReadConfig(
        filename="read.fastq",
        filepath=data_dir / "read.fastq",
        ngs_type="illumina",
        coding_type=CodingType.NONE,  # Uncompressed
        detected_format=ReadFormat.FASTQ,
        extra={}
    )

    # Create settings
    read_settings = ReadValidator.Settings()
    read_settings = read_settings.update(
        validation_level='trust',
        outdir_by_ngs_type=True,
        coding_type='gz'
    )

    # Create validator and validate
    read_validator = ReadValidator(read_config, output_dir, read_settings)
    read_validator.validate()

    print(f"Validated: {read_config.filename}")
    print(f"Settings: {read_settings}")
    print()

    # Example 3: Feature Validator
    print("3. Creating and using FeatureValidator directly")
    print("-" * 70)

    # Create feature config manually
    feature_config = FeatureConfig(
        filename="feature.gff",
        filepath=data_dir / "feature.gff",
        coding_type=CodingType.NONE,  # Uncompressed
        detected_format=FeatureFormat.GFF,
        extra={}
    )

    # Create settings
    feature_settings = FeatureValidator.Settings()
    feature_settings = feature_settings.update(
        validation_level='strict',
        sort_by_position=True,
        replace_id_with='chr1',
        coding_type='gz'
    )

    # Create validator and validate
    feature_validator = FeatureValidator(feature_config, output_dir, feature_settings)
    feature_validator.validate()

    print(f"Validated: {feature_config.filename}")
    print(f"Settings: {feature_settings}")
    print()

    # Example 4: Batch processing with compression formats
    print("4. Batch processing files with different compressions")
    print("-" * 70)

    # Define files to process
    genome_files = [
        ("genome.fasta", CodingType.NONE, GenomeFormat.FASTA),
        ("genome.fa.gz", CodingType.GZIP, GenomeFormat.FASTA),
    ]

    for filename, coding, format_type in genome_files:
        config = GenomeConfig(
            filename=filename,
            filepath=data_dir / filename,
            coding_type=coding,
            detected_format=format_type,
            extra={}
        )

        settings = GenomeValidator.Settings()
        settings = settings.update(
            validation_level='trust',
            coding_type='gz'  # Convert all to gzip
        )

        validator = GenomeValidator(config, output_dir, settings)
        validator.validate()

        print(f"  Processed: {filename} -> converted to .gz")

    print()
    print("=" * 70)
    print("5. Using Config Settings from config.json (NEW)")
    print("=" * 70)

    # Load a config file
    from validation_pkg.config_manager import ConfigManager
    config_file = Path(__file__).parent / "data" / "config.json"

    if config_file.exists():
        config = ConfigManager.load(str(config_file))

        # Example 5a: Use config settings automatically (RECOMMENDED)
        print("\nExample 5a: Config settings applied automatically:")
        print("-" * 70)
        if config.ref_genome:
            # Don't pass settings parameter - config settings will be used!
            validator = GenomeValidator(config.ref_genome, output_dir)
            print(f"  Settings applied from config:")
            print(f"    - validation_level: {validator.settings.validation_level}")
            print(f"    - threads: {validator.settings.threads}")
            print("  Validator will use all settings from config.json")

        # Example 5b: User settings REPLACE config settings (be careful!)
        print("\nExample 5b: User settings replace config settings:")
        print("-" * 70)
        if config.ref_genome:
            # When you pass settings, config settings are IGNORED
            user_settings = GenomeValidator.Settings(coding_type='bz2')
            validator = GenomeValidator(config.ref_genome, output_dir, user_settings)
            print(f"  User settings override config:")
            print(f"    - coding_type: {validator.settings.coding_type} (from user)")
            print(f"    - validation_level: {validator.settings.validation_level} (default, NOT from config)")
            print("  Config settings are NOT merged - user settings replace them entirely")

    print()
    print("=" * 70)
    print("Key Points:")
    print("=" * 70)
    print("1. You need to create Config objects manually (GenomeConfig, etc.)")
    print("2. You must specify: filename, filepath, coding_type, detected_format")
    print("3. CodingType: NONE, GZIP, or BZIP2")
    print("4. Format detection must be done manually or use ConfigManager")
    print("5. After creating the validator, call validate() to run")
    print()
    print("NEW: Config Settings Support")
    print("-" * 70)
    print("• Validators automatically apply settings from config.json")
    print("• Don't pass settings=... to use config settings")
    print("• Passing Settings object REPLACES config settings (no merge)")
    print()
    print(f"Output directory: {output_dir}")

if __name__ == "__main__":
    main()
