#!/usr/bin/env python3
"""
Example 1: Simple Validation Using Functional API

This example demonstrates the simplest way to use the validation package:
- Load a configuration file
- Validate all files specified in the config using the functional API
- Check output files

This is the recommended approach for most use cases.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from validation_pkg import (
    ConfigManager,
    validate_genome,
    validate_read,
    validate_feature
)

def main():
    # Path to the configuration file
    config_path = Path(__file__).parent / "data" / "config.json"

    # Output directory for validated files
    output_dir = Path(__file__).parent / "output" / "simple_validation"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Simple Validation Example")
    print("=" * 70)
    print(f"Config file: {config_path}")
    print(f"Output directory: {output_dir}")
    print()

    # Load configuration
    print("Loading configuration...")
    config = ConfigManager.load(str(config_path))
    print()

    # Validate genomes
    print("Validating genome files...")
    if config.ref_genome:
        validate_genome(config.ref_genome, output_dir)
        print(f"  ✓ Validated: {config.ref_genome.filename}")

    if config.mod_genome:
        validate_genome(config.mod_genome, output_dir)
        print(f"  ✓ Validated: {config.mod_genome.filename}")

    if config.ref_plasmid:
        validate_genome(config.ref_plasmid, output_dir)
        print(f"  ✓ Validated: {config.ref_plasmid.filename}")
    print()

    # Validate reads
    print("Validating read files...")
    if config.reads:
        for read in config.reads:
            validate_read(read, output_dir)
            print(f"  ✓ Validated: {read.filename}")
    else:
        print("  No read files to validate")
    print()

    # Validate features
    print("Validating feature files...")
    if config.ref_feature:
        validate_feature(config.ref_feature, output_dir)
        print(f"  ✓ Validated: {config.ref_feature.filename}")

    if config.mod_feature:
        validate_feature(config.mod_feature, output_dir)
        print(f"  ✓ Validated: {config.mod_feature.filename}")
    print()

    # Summary
    print("=" * 70)
    print("VALIDATION COMPLETE")
    print("=" * 70)
    print("All files have been validated successfully!")
    print(f"\nOutput files saved to: {output_dir}")
    print("\nNext steps:")
    print("  - Check the output directory for validated files")
    print("  - Review other examples for custom settings")
    print("=" * 70)

if __name__ == "__main__":
    main()
