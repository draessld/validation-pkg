"""
Command-line interface for validation package.

Provides CLI entry point for running validations from the command line.

Usage:
    python -m validation_pkg validate config.json
    python -m validation_pkg validate config.json --only genomes
    python -m validation_pkg validate config.json --output-dir ./output
"""

import sys
import argparse
from pathlib import Path

from validation_pkg.coordinator import ValidationCoordinator
from validation_pkg.logger import setup_logging


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Bioinformatics Validation Package - Validate genomic data files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate all files
  python -m validation_pkg validate config.json

  # Validate only genomes
  python -m validation_pkg validate config.json --only genomes

  # Validate only reads
  python -m validation_pkg validate config.json --only reads

  # Validate only features
  python -m validation_pkg validate config.json --only features

  # Specify output directory
  python -m validation_pkg validate config.json --output-dir ./output

  # Enable verbose logging
  python -m validation_pkg validate config.json --verbose
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Command to execute')

    # Validate command
    validate_parser = subparsers.add_parser(
        'validate',
        help='Validate genomic data files'
    )
    validate_parser.add_argument(
        'config',
        type=str,
        help='Path to configuration JSON file'
    )
    validate_parser.add_argument(
        '--only',
        type=str,
        choices=['genomes', 'reads', 'features'],
        help='Validate only specific file type'
    )
    validate_parser.add_argument(
        '--output-dir',
        type=str,
        help='Output directory for validated files (overrides config)'
    )
    validate_parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    validate_parser.add_argument(
        '--log-file',
        type=str,
        help='Write logs to file'
    )

    # Parse arguments
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    # Setup logging
    log_level = 'DEBUG' if args.verbose else 'INFO'
    log_file = Path(args.log_file) if args.log_file else None
    setup_logging(level=log_level, log_file=log_file)

    # Execute command
    if args.command == 'validate':
        return validate_command(args)

    return 0


def validate_command(args):
    """Execute validation command."""
    try:
        # Check config file exists
        config_path = Path(args.config)
        if not config_path.exists():
            print(f"Error: Config file not found: {config_path}", file=sys.stderr)
            return 1

        # Create coordinator
        coordinator = ValidationCoordinator(
            str(config_path),
            output_dir=args.output_dir
        )

        # Run validation based on --only flag
        if args.only == 'genomes':
            print("Validating genomes only...")
            report = coordinator.validate_genomes()
        elif args.only == 'reads':
            print("Validating reads only...")
            report = coordinator.validate_reads()
        elif args.only == 'features':
            print("Validating features only...")
            report = coordinator.validate_features()
        else:
            print("Validating all files...")
            report = coordinator.validate_all()

        # Print report
        print()
        print(report.summary())

        # Return exit code
        return 0 if report.passed else 1

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
