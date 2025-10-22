#!/usr/bin/env python3
"""
Example 8: Error Handling

This example demonstrates how to handle errors and exceptions
that may occur during validation:
- Catching specific validation errors
- Handling file format errors
- Dealing with configuration errors
- Proper exception handling patterns

This is important for building robust pipelines that can
gracefully handle errors.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from validation_pkg import ConfigManager, validate_genome, validate_read, validate_feature
from validation_pkg.exceptions import (
    ValidationError,
    ConfigurationError,
    FileFormatError,
    GenomeValidationError,
    ReadValidationError,
    FeatureValidationError
)

def main():
    output_dir = Path(__file__).parent / "output" / "error_handling"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Error Handling Example")
    print("=" * 70)
    print()

    # Example 1: Handling missing config file
    print("1. Handling missing configuration file")
    print("-" * 70)

    try:
        config = ConfigManager.load("nonexistent_config.json")
    except (ConfigurationError, FileNotFoundError) as e:
        print(f"Caught expected error: {type(e).__name__}")
        print(f"Message: {e}")
    print()

    # Example 2: Handling validation errors with error tracking
    print("2. Handling validation errors with error tracking")
    print("-" * 70)

    config_path = Path(__file__).parent / "data" / "config.json"

    try:
        config = ConfigManager.load(str(config_path))

        errors = []
        validated_files = []

        # Validate genomes with error tracking
        if config.ref_genome:
            try:
                validate_genome(config.ref_genome, output_dir)
                validated_files.append(config.ref_genome.filename)
                print(f"  ✓ {config.ref_genome.filename}")
            except ValidationError as e:
                errors.append(f"{config.ref_genome.filename}: {e}")
                print(f"  ✗ {config.ref_genome.filename}: {e}")

        # Validate reads with error tracking
        if config.reads:
            for read in config.reads:
                try:
                    validate_read(read, output_dir)
                    validated_files.append(read.filename)
                    print(f"  ✓ {read.filename}")
                except ValidationError as e:
                    errors.append(f"{read.filename}: {e}")
                    print(f"  ✗ {read.filename}: {e}")

        # Validate features with error tracking
        if config.ref_feature:
            try:
                validate_feature(config.ref_feature, output_dir)
                validated_files.append(config.ref_feature.filename)
                print(f"  ✓ {config.ref_feature.filename}")
            except ValidationError as e:
                errors.append(f"{config.ref_feature.filename}: {e}")
                print(f"  ✗ {config.ref_feature.filename}: {e}")

        print()
        if errors:
            print(f"Validation completed with {len(errors)} error(s):")
            for error in errors:
                print(f"  - {error}")
        else:
            print(f"All {len(validated_files)} file(s) validated successfully!")

    except ValidationError as e:
        print(f"Validation error occurred: {e}")
    except Exception as e:
        print(f"Unexpected error: {type(e).__name__}: {e}")
    print()

    # Example 3: Catching specific validation errors
    print("3. Catching specific validator errors")
    print("-" * 70)

    try:
        config = ConfigManager.load(str(config_path))

        # Try to validate genome with specific error handling
        if config.ref_genome:
            try:
                validate_genome(config.ref_genome, output_dir)
                print("Genome validation successful")
            except GenomeValidationError as e:
                print(f"Genome validation error: {e}")

    except FileFormatError as e:
        print(f"File format error: {e}")
    except ConfigurationError as e:
        print(f"Configuration error: {e}")
    except ValidationError as e:
        print(f"General validation error: {e}")
    print()

    # Example 4: Creating a robust validation function
    print("4. Robust validation with comprehensive error handling")
    print("-" * 70)

    def robust_validate(config_path, output_dir):
        """
        A robust validation function with comprehensive error handling.
        Returns: (success: bool, message: str, details: dict)
        """
        try:
            # Validate inputs
            config_path = Path(config_path)
            if not config_path.exists():
                return False, f"Config file not found: {config_path}", {}

            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            # Load configuration
            config = ConfigManager.load(str(config_path))

            # Track validation results
            validated = []
            failed = []

            # Validate genomes
            for genome in [config.ref_genome, config.mod_genome, config.ref_plasmid, config.mod_plasmid]:
                if genome:
                    try:
                        validate_genome(genome, output_dir)
                        validated.append(genome.filename)
                    except ValidationError as e:
                        failed.append((genome.filename, str(e)))

            # Validate reads
            if config.reads:
                for read in config.reads:
                    try:
                        validate_read(read, output_dir)
                        validated.append(read.filename)
                    except ValidationError as e:
                        failed.append((read.filename, str(e)))

            # Validate features
            for feature in [config.ref_feature, config.mod_feature]:
                if feature:
                    try:
                        validate_feature(feature, output_dir)
                        validated.append(feature.filename)
                    except ValidationError as e:
                        failed.append((feature.filename, str(e)))

            # Generate result
            if failed:
                error_msg = f"Validation failed: {len(failed)} error(s)"
                details = {
                    'validated': validated,
                    'failed': failed,
                    'total': len(validated) + len(failed)
                }
                return False, error_msg, details
            else:
                success_msg = f"Successfully validated {len(validated)} file(s)"
                details = {
                    'validated': validated,
                    'failed': [],
                    'total': len(validated)
                }
                return True, success_msg, details

        except ConfigurationError as e:
            return False, f"Configuration error: {e}", {}
        except FileFormatError as e:
            return False, f"File format error: {e}", {}
        except ValidationError as e:
            return False, f"Validation error: {e}", {}
        except Exception as e:
            return False, f"Unexpected error: {type(e).__name__}: {e}", {}

    # Use the robust function
    success, message, details = robust_validate(config_path, output_dir)
    print(f"Result: {'SUCCESS' if success else 'FAILURE'}")
    print(f"Message: {message}")
    if details:
        print(f"Details:")
        print(f"  - Total files: {details.get('total', 0)}")
        print(f"  - Validated: {len(details.get('validated', []))}")
        print(f"  - Failed: {len(details.get('failed', []))}")
    print()

    # Example 5: Exception hierarchy
    print("5. Understanding the exception hierarchy")
    print("-" * 70)
    print("ValidationError (base)")
    print("├── ConfigurationError")
    print("├── FileNotFoundError")
    print("├── FileFormatError")
    print("│   ├── FastaFormatError")
    print("│   ├── GenBankFormatError")
    print("│   ├── FastqFormatError")
    print("│   ├── BamFormatError")
    print("│   ├── GffFormatError")
    print("│   └── BedFormatError")
    print("├── CompressionError")
    print("├── GenomeValidationError")
    print("├── ReadValidationError")
    print("├── FeatureValidationError")
    print("└── InterFileValidationError")
    print()

    print("=" * 70)
    print("Best Practices:")
    print("=" * 70)
    print("1. Always catch specific exceptions before general ones")
    print("2. Track validated and failed files separately")
    print("3. Log errors for debugging and monitoring")
    print("4. Validate inputs before processing")
    print("5. Use try-except blocks around each validation call")
    print("6. Return structured results (success, message, details)")
    print("=" * 70)

if __name__ == "__main__":
    main()
