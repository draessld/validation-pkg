#!/usr/bin/env python3
"""
Example 9: Settings Patterns

This example demonstrates the immutable settings pattern used
throughout the validation package:
- Creating settings objects
- Updating settings (immutable pattern)
- Converting to/from dictionaries
- Inspecting settings
- Reusing settings across validators

Understanding this pattern is crucial for working effectively
with the package.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from validation_pkg.validators.genome_validator import GenomeValidator
from validation_pkg.validators.read_validator import ReadValidator
from validation_pkg.validators.feature_validator import FeatureValidator

def main():
    print("=" * 70)
    print("Settings Patterns Example")
    print("=" * 70)
    print()

    # Example 1: Creating and updating settings (IMMUTABLE pattern)
    print("1. Immutable settings pattern (IMPORTANT!)")
    print("-" * 70)

    # CORRECT way - assign the result of update()
    settings = GenomeValidator.Settings()
    print(f"Default settings:\n{settings}\n")

    # Update returns a NEW settings object - you must assign it!
    new_settings = settings.update(
        validation_level='trust',
        min_sequence_length=500,
        coding_type='gz'
    )

    print(f"Original settings (unchanged):\n{settings}\n")
    print(f"New settings (with updates):\n{new_settings}\n")

    # WRONG way - changes are lost!
    print("WRONG: Don't do this!")
    settings.update(min_sequence_length=1000)  # This does nothing!
    print(f"Settings still unchanged: min_sequence_length = {settings.min_sequence_length}\n")
    print()

    # Example 2: Chaining updates
    print("2. Chaining multiple updates")
    print("-" * 70)

    settings = GenomeValidator.Settings()
    settings = settings.update(
        validation_level='strict',
        min_sequence_length=100
    ).update(
        replace_id_with='chr',
        coding_type='gz'
    ).update(
        plasmid_split=True
    )

    print(f"Chained updates result:\n{settings}\n")
    print()

    # Example 3: Converting to/from dictionaries
    print("3. Converting settings to/from dictionaries")
    print("-" * 70)

    settings = ReadValidator.Settings()
    settings = settings.update(
        validation_level='trust',
        outdir_by_ngs_type=True,
        coding_type='gz'
    )

    # Convert to dictionary
    settings_dict = settings.to_dict()
    print(f"Settings as dictionary:")
    for key, value in settings_dict.items():
        print(f"  {key}: {value}")
    print()

    # Create settings from dictionary
    new_settings = ReadValidator.Settings.from_dict(settings_dict)
    print(f"Recreated from dict:\n{new_settings}\n")
    print()

    # Example 4: Creating reusable settings templates
    print("4. Creating reusable settings templates")
    print("-" * 70)

    # Define common settings templates
    STRICT_QC_TEMPLATE = GenomeValidator.Settings().update(
        validation_level='strict',
        min_sequence_length=500,
        coding_type='gz'
    )

    FAST_PRODUCTION_TEMPLATE = GenomeValidator.Settings().update(
        validation_level='trust',
        coding_type='gz'
    )

    ARCHIVE_TEMPLATE = GenomeValidator.Settings().update(
        validation_level='minimal'
    )

    print("Defined templates:")
    print(f"  STRICT_QC: {STRICT_QC_TEMPLATE.validation_level}, "
          f"min_len={STRICT_QC_TEMPLATE.min_sequence_length}")
    print(f"  FAST_PRODUCTION: {FAST_PRODUCTION_TEMPLATE.validation_level}")
    print(f"  ARCHIVE: {ARCHIVE_TEMPLATE.validation_level}")
    print()

    # Use templates with customization
    custom_settings = STRICT_QC_TEMPLATE.update(
        replace_id_with='custom_prefix',
        output_filename_suffix='processed'
    )

    print(f"Customized from template:\n{custom_settings}\n")
    print()

    # Example 5: Settings for different validators
    print("5. Settings for different validator types")
    print("-" * 70)

    # Each validator has its own Settings class with specific options

    # Genome settings
    genome_settings = GenomeValidator.Settings()
    genome_settings = genome_settings.update(
        plasmid_split=True,  # Specific to genome
        min_sequence_length=100
    )
    print(f"Genome settings:\n{genome_settings}\n")

    # Read settings
    read_settings = ReadValidator.Settings()
    read_settings = read_settings.update(
        outdir_by_ngs_type=True,  # Specific to reads
        check_invalid_chars=True
    )
    print(f"Read settings:\n{read_settings}\n")

    # Feature settings
    feature_settings = FeatureValidator.Settings()
    feature_settings = feature_settings.update(
        sort_by_position=True,  # Specific to features
        check_coordinates=True
    )
    print(f"Feature settings:\n{feature_settings}\n")
    print()

    # Example 6: Inspecting default values
    print("6. Inspecting default settings values")
    print("-" * 70)

    defaults = GenomeValidator.Settings()
    print("GenomeValidator default settings:")
    print(f"  validation_level: {defaults.validation_level}")
    print(f"  min_sequence_length: {defaults.min_sequence_length}")
    print(f"  allow_empty_sequences: {defaults.allow_empty_sequences}")
    print(f"  plasmid_split: {defaults.plasmid_split}")
    print(f"  coding_type: {defaults.coding_type}")
    print()

    # Example 7: Common patterns for different use cases
    print("7. Common settings patterns for different use cases")
    print("-" * 70)

    print("A. Development/QC (thorough validation):")
    dev_settings = GenomeValidator.Settings().update(
        validation_level='strict',
        min_sequence_length=100,
        allow_empty_sequences=False,
        coding_type='gz'
    )
    print(f"   {dev_settings.validation_level}, min_len={dev_settings.min_sequence_length}\n")

    print("B. Production (fast processing):")
    prod_settings = GenomeValidator.Settings().update(
        validation_level='trust',
        coding_type='gz'
    )
    print(f"   {prod_settings.validation_level}\n")

    print("C. Data archiving (no validation):")
    archive_settings = GenomeValidator.Settings().update(
        validation_level='minimal'
    )
    print(f"   {archive_settings.validation_level}\n")

    print("D. Format conversion (trust + compression change):")
    convert_settings = ReadValidator.Settings().update(
        validation_level='trust',
        coding_type='bz2'
    )
    print(f"   {convert_settings.validation_level}, output={convert_settings.coding_type}\n")

    print("=" * 70)
    print("Key Takeaways:")
    print("=" * 70)
    print("1. Settings are IMMUTABLE - always assign the result of update()")
    print("2. Use update() to modify settings, returns a NEW object")
    print("3. Can convert to/from dictionaries for serialization")
    print("4. Create reusable templates for common scenarios")
    print("5. Each validator type has its own Settings class")
    print("6. Settings objects have nice string representations")
    print("=" * 70)

if __name__ == "__main__":
    main()
