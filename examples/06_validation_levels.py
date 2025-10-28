#!/usr/bin/env python3
"""
Example 6: Understanding Validation Levels

This example demonstrates the three validation levels available:
- STRICT: Full parsing and validation (slowest, most thorough)
- TRUST: Selective validation (fast, for trusted data)
- MINIMAL: No validation (fastest, for archiving)

Each level has different use cases:
- STRICT: Development, QC, untrusted data
- TRUST: Production pipelines with previously validated data
- MINIMAL: Fast archiving/staging where validation is not needed
"""

import sys
from pathlib import Path
import time

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from validation_pkg import ConfigManager, GenomeValidator, validate_genome

def main():
    # Path to the configuration file
    config_path = Path(__file__).parent / "data" / "config.json"

    # Output directory
    output_dir = Path(__file__).parent / "output" / "validation_levels"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Validation Levels Comparison")
    print("=" * 70)
    print()

    # Load configuration
    config = ConfigManager.load(str(config_path))

    if not config.ref_genome:
        print("No genome file found in config")
        return

    # Level 1: STRICT
    print("1. STRICT VALIDATION LEVEL")
    print("-" * 70)
    print("Behavior:")
    print("  - Parse all sequences")
    print("  - Validate all sequences")
    print("  - Apply all edits")
    print("  - BioPython write")
    print("Use case: Development, QC, untrusted data")
    print()

    settings_strict = GenomeValidator.Settings()
    settings_strict = settings_strict.update(
        validation_level='strict',
        min_sequence_length=50,
        coding_type='gz',
        output_filename_suffix='strict'
    )

    start = time.time()
    validate_genome(config.ref_genome, output_dir, settings_strict)
    elapsed = time.time() - start

    print(f"Completed in {elapsed:.3f} seconds")
    print(f"Output: {config.ref_genome.filename}_strict.fasta.gz")
    print()

    # Level 2: TRUST
    print("2. TRUST VALIDATION LEVEL")
    print("-" * 70)
    print("Behavior:")
    print("  - Parse all sequences")
    print("  - Validate FIRST sequence only")
    print("  - Apply all edits")
    print("  - BioPython write")
    print("Use case: Production pipelines, trusted data, fast processing")
    print()

    settings_trust = GenomeValidator.Settings()
    settings_trust = settings_trust.update(
        validation_level='trust',
        coding_type='gz',
        output_filename_suffix='trust'
    )

    start = time.time()
    validate_genome(config.ref_genome, output_dir, settings_trust)
    elapsed = time.time() - start

    print(f"Completed in {elapsed:.3f} seconds")
    print(f"Output: {config.ref_genome.filename}_trust.fasta.gz")
    print()

    # Level 3: MINIMAL
    print("3. MINIMAL VALIDATION LEVEL")
    print("-" * 70)
    print("Behavior:")
    print("  - No parsing")
    print("  - No validation")
    print("  - No edits")
    print("  - File copy as-is")
    print("Requirements: Input/output format and compression must match")
    print("Use case: Fast archiving, staging, when validation is not needed")
    print()

    settings_minimal = GenomeValidator.Settings()
    settings_minimal = settings_minimal.update(
        validation_level='minimal',
        # Note: For minimal mode to work, the coding_type must match the input
        # Since our input is uncompressed, we don't set coding_type (defaults to None)
        output_filename_suffix='minimal'
    )

    start = time.time()
    validate_genome(config.ref_genome, output_dir, settings_minimal)
    elapsed = time.time() - start

    print(f"Completed in {elapsed:.3f} seconds")
    print(f"Output: {config.ref_genome.filename}_minimal.fasta")
    print()

    # Comparison
    print("=" * 70)
    print("PERFORMANCE COMPARISON")
    print("=" * 70)
    print("For large files (>10GB FASTQ):")
    print("  - TRUST is typically 10-15x faster than STRICT")
    print("  - MINIMAL is essentially instant (file copy only)")
    print()
    print("Recommendations:")
    print("  - Use STRICT for new/untrusted data and during development")
    print("  - Use TRUST for production with previously validated data")
    print("  - Use MINIMAL only for archiving when validation is not needed")
    print("=" * 70)
    print()

    # Tip about config-level settings
    print("💡 TIP: You can also specify validation_level in config.json!")
    print("Example config.json:")
    print('{')
    print('  "ref_genome_filename": {')
    print('    "filename": "genome.fasta",')
    print('    "validation_level": "trust"')
    print('  },')
    print('  "reads": [')
    print('    {')
    print('      "filename": "reads.fastq.gz",')
    print('      "ngs_type": "illumina",')
    print('      "validation_level": "minimal"')
    print('    }')
    print('  ]')
    print('}')
    print()
    print("See example 12_config_validation_levels.py for config-based approach.")

if __name__ == "__main__":
    main()
