#!/usr/bin/env python3
"""
Example 5: Custom Feature Validation Settings

This example demonstrates advanced feature file validation:
- Coordinate validation
- Position-based sorting
- Sequence name replacement
- Format conversion (BED to GFF)
- Validation levels

This is useful for processing annotation files from different sources
and ensuring they are properly formatted and sorted.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from validation_pkg import ConfigManager, FeatureValidator, validate_feature

def main():
    # Path to the configuration file
    config_path = Path(__file__).parent / "data" / "config.json"

    # Output directory
    output_dir = Path(__file__).parent / "output" / "custom_features"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Custom Feature Validation Settings Example")
    print("=" * 70)
    print()

    # Load configuration
    config = ConfigManager.load(str(config_path))

    # Example 1: Strict validation with sorting
    print("1. STRICT validation with position-based sorting")
    print("-" * 70)

    settings_sorted = FeatureValidator.Settings()
    settings_sorted = settings_sorted.update(
        validation_level='strict',
        sort_by_position=True,
        check_coordinates=True,
        allow_zero_length=False,
        coding_type='gz',
        output_filename_suffix='sorted'
    )

    print(f"Settings: {settings_sorted}")
    print()

    if config.ref_feature:
        validate_feature(config.ref_feature, output_dir, settings_sorted)
        print(f"Validated and sorted: {config.ref_feature.filename}")
    print()

    # Example 2: Sequence name replacement
    print("2. Replace sequence names in features")
    print("-" * 70)

    settings_renamed = FeatureValidator.Settings()
    settings_renamed = settings_renamed.update(
        validation_level='strict',
        replace_id_with='chromosome',
        sort_by_position=True,
        coding_type='gz',
        output_filename_suffix='renamed'
    )

    print(f"Settings: {settings_renamed}")
    print("All sequence names will be replaced with 'chromosome' prefix")
    print()

    if config.ref_feature:
        validate_feature(config.ref_feature, output_dir, settings_renamed)
        print(f"Validated with renamed sequences: {config.ref_feature.filename}")
    print()

    # Example 3: Trust mode for fast processing
    print("3. TRUST validation for fast processing")
    print("-" * 70)

    settings_trust = FeatureValidator.Settings()
    settings_trust = settings_trust.update(
        validation_level='trust',
        coding_type='gz',
        output_filename_suffix='trust'
    )

    print(f"Settings: {settings_trust}")
    print()

    if config.mod_feature:
        validate_feature(config.mod_feature, output_dir, settings_trust)
        print(f"Validated (trust mode): {config.mod_feature.filename}")
    print()

    # Example 4: Coordinate validation without sorting
    print("4. Coordinate validation WITHOUT sorting")
    print("-" * 70)

    settings_no_sort = FeatureValidator.Settings()
    settings_no_sort = settings_no_sort.update(
        validation_level='strict',
        sort_by_position=False,
        check_coordinates=True,
        coding_type='gz',
        output_filename_suffix='validated_unsorted'
    )

    print(f"Settings: {settings_no_sort}")
    print("Features will be validated but maintain their original order")
    print()

    if config.ref_feature:
        validate_feature(config.ref_feature, output_dir, settings_no_sort)
        print(f"Validated without sorting: {config.ref_feature.filename}")
    print()

    print("=" * 70)
    print(f"All validated files saved to: {output_dir}")
    print("=" * 70)

if __name__ == "__main__":
    main()
