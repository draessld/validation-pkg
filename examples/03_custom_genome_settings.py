#!/usr/bin/env python3
"""
Example 3: Custom Genome Validation Settings

This example demonstrates advanced genome validation with custom settings:
- Setting validation levels (strict, trust, minimal)
- Plasmid handling (split, merge)
- Sequence filtering (minimum length)
- ID replacement (adding prefixes)
- Output format control (compression)

This is useful for production pipelines where you need fine control
over the validation process.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from validation_pkg import ConfigManager, GenomeValidator, validate_genome

def main():
    # Path to the configuration file
    config_path = Path(__file__).parent / "data" / "config.json"

    # Output directory
    output_dir = Path(__file__).parent / "output" / "custom_genome"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Custom Genome Validation Settings Example")
    print("=" * 70)
    print()

    # Load configuration
    config = ConfigManager.load(str(config_path))

    # Example 1: Strict validation with sequence filtering
    print("1. STRICT validation with minimum sequence length filter")
    print("-" * 70)

    settings_strict = GenomeValidator.Settings()
    settings_strict = settings_strict.update(
        validation_level='strict',
        min_sequence_length=100,
        replace_id_with='chr',
        coding_type='gz',
        output_filename_suffix='strict_filtered'
    )

    print(f"Settings: {settings_strict}")
    print()

    if config.ref_genome:
        validate_genome(config.ref_genome, output_dir, settings_strict)
        print(f"Validated: {config.ref_genome.filename}")
    print()

    # Example 2: Trust mode for fast processing
    print("2. TRUST validation for fast processing")
    print("-" * 70)

    settings_trust = GenomeValidator.Settings()
    settings_trust = settings_trust.update(
        validation_level='trust',
        coding_type='gz',
        output_filename_suffix='trust_fast'
    )

    print(f"Settings: {settings_trust}")
    print()

    if config.mod_genome:
        validate_genome(config.mod_genome, output_dir, settings_trust)
        print(f"Validated: {config.mod_genome.filename}")
    print()

    # Example 3: Plasmid splitting
    print("3. Plasmid handling with SPLIT mode")
    print("-" * 70)

    settings_plasmid = GenomeValidator.Settings()
    settings_plasmid = settings_plasmid.update(
        validation_level='strict',
        plasmid_split=True,
        coding_type='gz',
        output_subdir_name='plasmids'
    )

    print(f"Settings: {settings_plasmid}")
    print()

    if config.ref_plasmid:
        validate_genome(config.ref_plasmid, output_dir, settings_plasmid)
        print(f"Validated: {config.ref_plasmid.filename}")
        print("Each plasmid saved to separate file in 'plasmids' subdirectory")
    print()

    # Example 4: Minimal mode for archiving
    print("4. MINIMAL validation for fast archiving")
    print("-" * 70)

    settings_minimal = GenomeValidator.Settings()
    settings_minimal = settings_minimal.update(
        validation_level='minimal',
        coding_type='gz'
    )

    print(f"Settings: {settings_minimal}")
    print("Note: Minimal mode just copies the file - input/output compression must match")
    print()

    print("=" * 70)
    print(f"All validated files saved to: {output_dir}")
    print("=" * 70)

if __name__ == "__main__":
    main()
