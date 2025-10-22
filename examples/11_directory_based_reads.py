#!/usr/bin/env python3
"""
Example 11: Directory-Based Read Validation

This example demonstrates how to validate all read files in a directory:
- Using the "directory" option instead of individual "filename" entries
- Applying the same settings to all files in a directory
- Different settings for different sequencing platforms
- Automatic file discovery

This is very useful when you have many read files from the same sequencing run
and want to process them all with the same settings.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from validation_pkg import ConfigManager, validate_read, ReadValidator

def main():
    # Path to the configuration file with directory entries
    config_path = Path(__file__).parent / "data" / "config_with_directories.json"

    # Output directory
    output_dir = Path(__file__).parent / "output" / "directory_reads"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Directory-Based Read Validation Example")
    print("=" * 70)
    print(f"Config file: {config_path}")
    print(f"Output directory: {output_dir}")
    print()

    # Load configuration
    print("Loading configuration...")
    config = ConfigManager.load(str(config_path))
    print()

    # Show what was loaded
    print("Configuration Summary:")
    print("-" * 70)
    print(f"Total read entries: {len(config.reads)}")
    for idx, read in enumerate(config.reads, 1):
        print(f"\nRead entry {idx}:")
        print(f"  Filename: {read.filename}")
        print(f"  NGS type: {read.ngs_type}")
        print(f"  Format: {read.detected_format}")
        print(f"  Compression: {read.coding_type}")
    print()

    # Example 1: Validate all reads with default settings
    print("=" * 70)
    print("Example 1: Validating all reads with default settings")
    print("=" * 70)
    print()

    for read in config.reads:
        validate_read(read, output_dir)
        print(f"✓ Validated: {read.filename}")

    print()
    print(f"All {len(config.reads)} read file(s) validated")
    print()

    # Example 2: Validate with custom settings per NGS type
    print("=" * 70)
    print("Example 2: Custom settings per NGS type")
    print("=" * 70)
    print()

    # Create different settings for different platforms
    illumina_settings = ReadValidator.Settings()
    illumina_settings = illumina_settings.update(
        validation_level='trust',  # Fast validation
        coding_type='gz',           # Convert to gzip
        outdir_by_ngs_type=True,   # Organize by platform
        output_filename_suffix='processed'
    )

    ont_settings = ReadValidator.Settings()
    ont_settings = ont_settings.update(
        validation_level='strict',     # Thorough validation
        check_invalid_chars=True,      # Check for invalid bases
        coding_type='gz',              # Convert to gzip
        outdir_by_ngs_type=True,      # Organize by platform
        output_filename_suffix='qc_passed'
    )

    output_dir2 = Path(__file__).parent / "output" / "directory_reads_custom"
    output_dir2.mkdir(parents=True, exist_ok=True)

    # Process with appropriate settings
    illumina_count = 0
    ont_count = 0

    for read in config.reads:
        if read.ngs_type == 'illumina':
            validate_read(read, output_dir2, illumina_settings)
            print(f"✓ Illumina (trust mode): {read.filename}")
            illumina_count += 1
        elif read.ngs_type == 'ont':
            validate_read(read, output_dir2, ont_settings)
            print(f"✓ ONT (strict mode): {read.filename}")
            ont_count += 1
        else:
            validate_read(read, output_dir2)
            print(f"✓ Other: {read.filename}")

    print()
    print(f"Processed {illumina_count} Illumina file(s) and {ont_count} ONT file(s)")
    print()

    # Summary
    print("=" * 70)
    print("DIRECTORY-BASED VALIDATION ADVANTAGES")
    print("=" * 70)
    print()
    print("Benefits of using directory option:")
    print("  1. No need to list every file individually in config")
    print("  2. Automatically processes all files in directory")
    print("  3. Apply same settings to all files in directory")
    print("  4. Easy to add new files - just drop them in the directory")
    print("  5. Ideal for batch processing of sequencing runs")
    print()
    print("Config file structure:")
    print('  {')
    print('    "reads": [')
    print('      {')
    print('        "directory": "illumina_reads/",')
    print('        "ngs_type": "illumina",')
    print('        "validation_level": "trust"')
    print('      },')
    print('      {')
    print('        "directory": "ont_reads/",')
    print('        "ngs_type": "ont",')
    print('        "validation_level": "strict"')
    print('      }')
    print('    ]')
    print('  }')
    print()
    print("=" * 70)
    print(f"Output saved to:")
    print(f"  - {output_dir}")
    print(f"  - {output_dir2}")
    print("=" * 70)

if __name__ == "__main__":
    main()
