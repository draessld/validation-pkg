#!/usr/bin/env python3
"""
Example 10: Complete Validation Pipeline

This example demonstrates a complete, production-ready validation pipeline
that combines all the concepts from previous examples:
- Configuration loading
- Custom settings for different file types
- Error handling
- Logging
- Progress reporting
- Output organization

This serves as a template for building your own validation pipelines.
"""

import sys
from pathlib import Path
import time

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from validation_pkg import (
    ValidationCoordinator,
    ConfigManager,
    GenomeValidator,
    ReadValidator,
    FeatureValidator,
    validate_genome,
    validate_reads,
    validate_feature
)
from validation_pkg.exceptions import ValidationError

def main():
    # Configuration
    config_path = Path(__file__).parent / "data" / "config.json"
    output_base = Path(__file__).parent / "output" / "complete_pipeline"
    output_base.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("COMPLETE VALIDATION PIPELINE")
    print("=" * 70)
    print(f"Config: {config_path}")
    print(f"Output: {output_base}")
    print("=" * 70)
    print()

    start_time = time.time()

    try:
        # Load configuration
        print("Loading configuration...")
        config = ConfigManager.load(str(config_path))
        print(f"  Found {len([x for x in [config.ref_genome, config.mod_genome] if x])} genome(s)")
        print(f"  Found {len(config.reads)} read file(s)")
        print(f"  Found {len([x for x in [config.ref_feature, config.mod_feature] if x])} feature file(s)")
        print()

        # Stage 1: Validate genomes
        print("-" * 70)
        print("STAGE 1: Genome Validation")
        print("-" * 70)

        genome_output = output_base / "genomes"
        genome_output.mkdir(exist_ok=True)

        # Reference genome: strict validation with QC
        if config.ref_genome:
            print(f"Validating reference genome: {config.ref_genome.filename}")
            ref_settings = GenomeValidator.Settings().update(
                validation_level='strict',
                min_sequence_length=100,
                replace_id_with='ref',
                coding_type='gz',
                output_filename_suffix='validated'
            )
            validate_genome(config.ref_genome, genome_output, ref_settings)
            print("  ✓ Reference genome validated")

        # Modified genome: trust mode for faster processing
        if config.mod_genome:
            print(f"Validating modified genome: {config.mod_genome.filename}")
            mod_settings = GenomeValidator.Settings().update(
                validation_level='trust',
                coding_type='gz',
                output_filename_suffix='validated'
            )
            validate_genome(config.mod_genome, genome_output, mod_settings)
            print("  ✓ Modified genome validated")

        # Plasmid: split into separate files
        if config.ref_plasmid:
            print(f"Validating plasmid: {config.ref_plasmid.filename}")
            plasmid_settings = GenomeValidator.Settings().update(
                validation_level='strict',
                plasmid_split=True,
                coding_type='gz',
                output_subdir_name='plasmids'
            )
            validate_genome(config.ref_plasmid, genome_output, plasmid_settings)
            print("  ✓ Plasmid validated and split")

        print()

        # Stage 2: Validate reads
        print("-" * 70)
        print("STAGE 2: Read Validation")
        print("-" * 70)

        read_output = output_base / "reads"
        read_output.mkdir(exist_ok=True)

        if config.reads:
            # Different settings based on NGS type
            illumina_settings = ReadValidator.Settings().update(
                validation_level='trust',
                outdir_by_ngs_type=True,
                coding_type='gz'
            )

            long_read_settings = ReadValidator.Settings().update(
                validation_level='strict',
                check_invalid_chars=True,
                outdir_by_ngs_type=True,
                coding_type='gz'
            )

            for read_config in config.reads:
                print(f"Validating {read_config.ngs_type} reads: {read_config.filename}")

                # Choose settings based on NGS type
                if read_config.ngs_type == 'illumina':
                    settings = illumina_settings
                    mode = "trust (fast)"
                else:
                    settings = long_read_settings
                    mode = "strict (thorough)"

                validate_reads([read_config], read_output, settings)
                print(f"  ✓ Validated using {mode} mode")
        else:
            print("No read files to validate")

        print()

        # Stage 3: Validate features
        print("-" * 70)
        print("STAGE 3: Feature Validation")
        print("-" * 70)

        feature_output = output_base / "features"
        feature_output.mkdir(exist_ok=True)

        # Reference features: strict validation with sorting
        if config.ref_feature:
            print(f"Validating reference features: {config.ref_feature.filename}")
            ref_feat_settings = FeatureValidator.Settings().update(
                validation_level='strict',
                sort_by_position=True,
                check_coordinates=True,
                replace_id_with='ref',
                coding_type='gz',
                output_filename_suffix='validated'
            )
            validate_feature(config.ref_feature, feature_output, ref_feat_settings)
            print("  ✓ Reference features validated and sorted")

        # Modified features: trust mode
        if config.mod_feature:
            print(f"Validating modified features: {config.mod_feature.filename}")
            mod_feat_settings = FeatureValidator.Settings().update(
                validation_level='trust',
                sort_by_position=True,
                coding_type='gz',
                output_filename_suffix='validated'
            )
            validate_feature(config.mod_feature, feature_output, mod_feat_settings)
            print("  ✓ Modified features validated")

        print()

        # Summary
        elapsed = time.time() - start_time
        print("=" * 70)
        print("PIPELINE COMPLETED SUCCESSFULLY")
        print("=" * 70)
        print(f"Total time: {elapsed:.2f} seconds")
        print()
        print("Output structure:")
        print(f"  {output_base}/")
        print(f"    ├── genomes/")
        print(f"    │   ├── validated genome files")
        print(f"    │   └── plasmids/ (if plasmid_split enabled)")
        print(f"    ├── reads/")
        print(f"    │   ├── illumina/ (if outdir_by_ngs_type enabled)")
        print(f"    │   ├── ont/")
        print(f"    │   └── pacbio/")
        print(f"    └── features/")
        print(f"        └── validated feature files")
        print()
        print("All files have been:")
        print("  - Validated according to their file type")
        print("  - Converted to gzip compression (.gz)")
        print("  - Organized into appropriate subdirectories")
        print("  - Named with '_validated' suffix")
        print("=" * 70)

        return 0

    except ValidationError as e:
        print()
        print("=" * 70)
        print("PIPELINE FAILED")
        print("=" * 70)
        print(f"Error: {e}")
        print("=" * 70)
        return 1

    except Exception as e:
        print()
        print("=" * 70)
        print("UNEXPECTED ERROR")
        print("=" * 70)
        print(f"Error: {type(e).__name__}: {e}")
        print("=" * 70)
        return 1

if __name__ == "__main__":
    sys.exit(main())
