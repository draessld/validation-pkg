#!/usr/bin/env python3
"""
Example 4: Custom Read Validation Settings

This example demonstrates advanced read validation with custom settings:
- Validation levels for different NGS types
- Automatic organization by NGS type
- Quality checks (invalid characters, duplicate IDs)
- Compression format conversion
- BAM file handling

This is useful for processing sequencing data from different platforms
with different quality requirements.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from validation_pkg import ConfigManager, ReadValidator, validate_read, validate_reads

def main():
    # Path to the configuration file
    config_path = Path(__file__).parent / "data" / "config.json"

    # Output directory
    output_dir = Path(__file__).parent / "output" / "custom_reads"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Custom Read Validation Settings Example")
    print("=" * 70)
    print()

    # Load configuration
    config = ConfigManager.load(str(config_path))

    # Example 1: Strict validation with quality checks
    print("1. STRICT validation with quality checks")
    print("-" * 70)

    settings_strict = ReadValidator.Settings()
    settings_strict = settings_strict.update(
        validation_level='strict',
        check_invalid_chars=True,
        allow_duplicate_ids=False,
        allow_empty_id=False,
        coding_type='gz'
    )

    print(f"Settings: {settings_strict}")
    print()

    if config.reads and len(config.reads) > 0:
        validate_read(config.reads[0], output_dir, settings_strict)
        print(f"Validated: {config.reads[0].filename}")
    print()

    # Example 2: Trust mode with automatic organization by NGS type
    print("2. TRUST validation with automatic NGS type organization")
    print("-" * 70)

    settings_trust = ReadValidator.Settings()
    settings_trust = settings_trust.update(
        validation_level='trust',
        outdir_by_ngs_type=True,
        coding_type='gz'
    )

    print(f"Settings: {settings_trust}")
    print("Files will be organized into subdirectories by NGS type (illumina, ont, pacbio)")
    print()

    if config.reads:
        validate_reads(config.reads, output_dir, settings_trust)
        print(f"Validated {len(config.reads)} read file(s)")
        print("Check subdirectories: illumina/, ont/, pacbio/")
    print()

    # Example 3: Different settings for different NGS types
    print("3. Custom settings per NGS type")
    print("-" * 70)

    # Illumina settings - trust mode is sufficient
    illumina_settings = ReadValidator.Settings()
    illumina_settings = illumina_settings.update(
        validation_level='trust',
        check_invalid_chars=False,
        coding_type='gz'
    )

    # ONT/PacBio settings - strict validation for long reads
    long_read_settings = ReadValidator.Settings()
    long_read_settings = long_read_settings.update(
        validation_level='strict',
        check_invalid_chars=True,
        coding_type='gz'
    )

    if config.reads:
        for read_config in config.reads:
            if read_config.ngs_type == 'illumina':
                settings = illumina_settings
                print(f"  Using TRUST mode for Illumina: {read_config.filename}")
            else:
                settings = long_read_settings
                print(f"  Using STRICT mode for {read_config.ngs_type}: {read_config.filename}")

            validate_read(read_config, output_dir, settings)
    print()

    # Example 4: Fast compression conversion
    print("4. Fast compression conversion (TRUST mode)")
    print("-" * 70)

    settings_convert = ReadValidator.Settings()
    settings_convert = settings_convert.update(
        validation_level='trust',
        coding_type='bz2',  # Convert to bzip2
        output_filename_suffix='converted'
    )

    print(f"Settings: {settings_convert}")
    print("Converting all reads to bzip2 compression")
    print()

    if config.reads:
        validate_reads(config.reads, output_dir, settings_convert)
        print(f"Converted {len(config.reads)} file(s) to .bz2 format")
    print()

    print("=" * 70)
    print(f"All validated files saved to: {output_dir}")
    print("=" * 70)
    print()

    # Tip about config-level settings
    print("💡 TIP: You can also specify these settings in config.json!")
    print("Example config.json:")
    print('{')
    print('  "reads": [')
    print('    {')
    print('      "filename": "reads.fastq.gz",')
    print('      "ngs_type": "illumina",')
    print('      "validation_level": "trust",')
    print('      "check_invalid_chars": false,')
    print('      "outdir_by_ngs_type": true')
    print('    }')
    print('  ]')
    print('}')
    print()
    print("See example 12_config_validation_levels.py for details on config-level settings.")

if __name__ == "__main__":
    main()
