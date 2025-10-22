#!/usr/bin/env python3
"""
Example 2: Selective Validation

This example demonstrates how to validate only specific file types:
- Validate only genome files
- Validate only read files
- Validate only feature files

This is useful when you want to process different file types separately
or validate only a subset of your data.
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
    output_dir = Path(__file__).parent / "output" / "selective_validation"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Selective Validation Example")
    print("=" * 70)
    print()

    # Load configuration
    config = ConfigManager.load(str(config_path))

    # Counter for validated files
    genome_count = 0
    read_count = 0
    feature_count = 0

    # Validate only genomes
    print("1. Validating GENOME files...")
    print("-" * 70)
    if config.ref_genome:
        validate_genome(config.ref_genome, output_dir)
        print(f"  ✓ {config.ref_genome.filename}")
        genome_count += 1
    if config.mod_genome:
        validate_genome(config.mod_genome, output_dir)
        print(f"  ✓ {config.mod_genome.filename}")
        genome_count += 1
    if config.ref_plasmid:
        validate_genome(config.ref_plasmid, output_dir)
        print(f"  ✓ {config.ref_plasmid.filename}")
        genome_count += 1
    print(f"Validated {genome_count} genome file(s)")
    print()

    # Validate only reads
    print("2. Validating READ files...")
    print("-" * 70)
    if config.reads:
        for read in config.reads:
            validate_read(read, output_dir)
            print(f"  ✓ {read.filename}")
            read_count += 1
        print(f"Validated {read_count} read file(s)")
    else:
        print("No read files to validate")
    print()

    # Validate only features
    print("3. Validating FEATURE files...")
    print("-" * 70)
    if config.ref_feature:
        validate_feature(config.ref_feature, output_dir)
        print(f"  ✓ {config.ref_feature.filename}")
        feature_count += 1
    if config.mod_feature:
        validate_feature(config.mod_feature, output_dir)
        print(f"  ✓ {config.mod_feature.filename}")
        feature_count += 1
    print(f"Validated {feature_count} feature file(s)")
    print()

    # Summary
    print("=" * 70)
    print("OVERALL SUMMARY")
    print("=" * 70)
    total_files = genome_count + read_count + feature_count

    print(f"Total files validated: {total_files}")
    print(f"  - Genomes: {genome_count}")
    print(f"  - Reads: {read_count}")
    print(f"  - Features: {feature_count}")
    print(f"\nOutput directory: {output_dir}")

if __name__ == "__main__":
    main()
