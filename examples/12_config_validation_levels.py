#!/usr/bin/env python3
"""
Example 12: Using Validation Levels in Config Files

This example demonstrates how to specify validation_level and other validator
settings directly in config.json. This is useful for production pipelines where
you want to configure validation behavior without changing code.

Key Features:
- Specify validation_level (strict/trust/minimal) per file in config.json
- Add any validator setting (plasmid_split, check_invalid_chars, etc.)
- Settings are automatically applied when using functional API
- Code-based settings override config settings if needed
"""

from validation_pkg import ConfigManager, validate_genome, validate_read, validate_feature
from validation_pkg import GenomeValidator, ReadValidator
from pathlib import Path

def example_basic_config_settings():
    """Example: Basic usage with validation_level in config."""
    print("=" * 70)
    print("Example 1: Basic Config-Level Validation Levels")
    print("=" * 70)

    # This config has validation_level specified for different files:
    # - ref_genome: trust mode
    # - reads: mix of strict, trust, and minimal
    # - features: strict mode with sorting

    config = ConfigManager.load("examples/data/config_with_validation_levels.json")

    # Check what settings were extracted
    print(f"\nRef genome settings from config:")
    print(f"  {config.ref_genome.settings_dict}")

    print(f"\nRead settings from config:")
    for idx, read in enumerate(config.reads):
        print(f"  Read {idx + 1}: {read.settings_dict}")

    # Validate - settings are applied automatically!
    print("\n✓ Validating with config-level settings...")

    if config.ref_genome:
        validate_genome(config.ref_genome, config.output_dir)
        print(f"  - Genome validated with validation_level: {config.ref_genome.settings_dict.get('validation_level', 'strict')}")

    if config.reads:
        for idx, read in enumerate(config.reads):
            validate_read(read, config.output_dir)
            level = read.settings_dict.get('validation_level', 'strict')
            print(f"  - Read {idx + 1} validated with validation_level: {level}")

    if config.ref_feature:
        validate_feature(config.ref_feature, config.output_dir)
        print(f"  - Feature validated with validation_level: {config.ref_feature.settings_dict.get('validation_level', 'strict')}")


def example_override_config_with_code():
    """Example: Override config settings with code-based settings."""
    print("\n" + "=" * 70)
    print("Example 2: Override Config Settings with Code")
    print("=" * 70)

    config = ConfigManager.load("examples/data/config_with_validation_levels.json")

    print(f"\nConfig says: validation_level = {config.ref_genome.settings_dict.get('validation_level')}")

    # Override with code settings (code takes precedence)
    code_settings = GenomeValidator.Settings(
        validation_level='strict',  # Override config's 'trust'
        min_sequence_length=1000    # Override config's 100
    )

    print(f"Code says: validation_level = strict, min_sequence_length = 1000")
    print("\n✓ Code settings take precedence!")

    # This will use 'strict' mode, not 'trust'
    validate_genome(config.ref_genome, config.output_dir, settings=code_settings)


def example_partial_override():
    """Example: Partially override config settings."""
    print("\n" + "=" * 70)
    print("Example 3: Partial Override - Mix Config and Code Settings")
    print("=" * 70)

    config = ConfigManager.load("examples/data/config_with_validation_levels.json")

    print(f"\nConfig provides:")
    print(f"  - validation_level: {config.ref_genome.settings_dict.get('validation_level')}")
    print(f"  - min_sequence_length: {config.ref_genome.settings_dict.get('min_sequence_length')}")

    # Start with defaults, then update just what you want
    code_settings = GenomeValidator.Settings()
    code_settings = code_settings.update(
        plasmid_split=True  # Add this setting in code
        # validation_level and min_sequence_length come from config
    )

    print(f"\nCode adds: plasmid_split = True")
    print("✓ Result: Uses config's validation_level + min_sequence_length, plus code's plasmid_split")

    # Config settings are applied, then code settings on top
    validate_genome(config.ref_genome, config.output_dir, settings=code_settings)


def example_directory_reads_inherit_settings():
    """Example: Directory-based reads inherit validation_level."""
    print("\n" + "=" * 70)
    print("Example 4: Directory Reads Inherit Settings")
    print("=" * 70)

    config = ConfigManager.load("examples/data/config_with_directories.json")

    print(f"\nConfig specifies directory: 'illumina_reads/' with:")
    print(f"  - ngs_type: illumina")
    print(f"  - validation_level: trust")

    print(f"\nAll files in directory inherit these settings:")
    for idx, read in enumerate(config.reads):
        if 'illumina' in str(read.filepath.parent):
            print(f"  - {read.filename}: {read.settings_dict}")


def example_multiple_settings_per_file():
    """Example: Specify multiple settings per file."""
    print("\n" + "=" * 70)
    print("Example 5: Multiple Settings Per File")
    print("=" * 70)

    config = ConfigManager.load("examples/data/config_with_validation_levels.json")

    print("\nGenome file with multiple settings:")
    print(f"  {config.ref_genome.settings_dict}")
    print("\n  All these settings are applied automatically!")

    # Find a read with multiple settings
    for read in config.reads:
        if len(read.settings_dict) > 1:
            print(f"\nRead file with multiple settings:")
            print(f"  {read.settings_dict}")
            print("\n  Both validation_level AND check_invalid_chars are applied!")
            break


def example_settings_precedence():
    """Example: Demonstrate settings precedence rules."""
    print("\n" + "=" * 70)
    print("Example 6: Settings Precedence Rules")
    print("=" * 70)

    config = ConfigManager.load("examples/data/config_with_validation_levels.json")

    print("\nPrecedence: Config < Code")
    print("\nScenario 1: No code settings")
    print("  → Uses config settings")
    print(f"  → validation_level = {config.ref_genome.settings_dict.get('validation_level')}")

    print("\nScenario 2: Code provides different validation_level")
    print("  → Code overrides config")
    code_settings = GenomeValidator.Settings(validation_level='minimal')
    print(f"  → validation_level = minimal (from code)")

    print("\nScenario 3: Code provides some settings, config provides others")
    print("  → Code settings override matching config settings")
    print("  → Non-matching config settings still apply")
    code_settings = GenomeValidator.Settings(validation_level='minimal')
    print(f"  → validation_level = minimal (from code)")
    print(f"  → min_sequence_length = {config.ref_genome.settings_dict.get('min_sequence_length')} (from config)")


def main():
    """Run all examples."""
    print("\n" + "=" * 70)
    print("CONFIG-LEVEL VALIDATION LEVELS - COMPREHENSIVE EXAMPLES")
    print("=" * 70)

    print("\nThis example demonstrates the new config-level settings feature.")
    print("You can now specify validation_level and other settings in config.json!")

    try:
        example_basic_config_settings()
        example_override_config_with_code()
        example_partial_override()
        example_directory_reads_inherit_settings()
        example_multiple_settings_per_file()
        example_settings_precedence()

        print("\n" + "=" * 70)
        print("✅ ALL EXAMPLES COMPLETED SUCCESSFULLY!")
        print("=" * 70)

        print("\nKey Takeaways:")
        print("  1. Specify any validator setting in config.json")
        print("  2. Settings are automatically applied by functional API")
        print("  3. Code settings override config settings (code wins)")
        print("  4. Directory reads inherit all settings")
        print("  5. Mix and match config and code settings as needed")

    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
