#!/usr/bin/env python3
"""
Main entry point for bioinformatics validation package.

Usage:
    python main.py config.json
    python main.py config.json --verbose
    python main.py config.json --log-dir logs/
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime

from validation_pkg.logger import setup_logging, get_logger
from validation_pkg.coordinator import Coordinator
from validation_pkg.exceptions import ValidationError, ConfigurationError, FileNotFoundError


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Bioinformatics File Validator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py config.json
  python main.py config.json --verbose
  python main.py config.json --log-dir results/logs/
  python main.py config.json --no-report
        """
    )
    
    parser.add_argument(
        'config',
        type=str,
        help='Path to configuration file (config.json)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output (DEBUG level)'
    )
    
    parser.add_argument(
        '--log-dir',
        type=str,
        default='logs',
        help='Directory for log files (default: logs/)'
    )
    
    parser.add_argument(
        '--no-report',
        action='store_true',
        help='Do not generate validation report'
    )
    
    parser.add_argument(
        '--no-log-file',
        action='store_true',
        help='Do not create log file (console only)'
    )
    
    return parser.parse_args()


def run_validation(config_path: str, args):
    """
    Run the validation pipeline.
    
    Args:
        config_path: Path to config.json
        args: Command line arguments
        
    Returns:
        int: Exit code (0 = success, 1 = failure)
    """
    # Setup logging
    log_dir = Path(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = None if args.no_log_file else log_dir / f"validation_{timestamp}.log"
    report_file = None if args.no_report else log_dir / f"report_{timestamp}.txt"
    
    console_level = "DEBUG" if args.verbose else "INFO"
    
    logger = setup_logging(
        console_level=console_level,
        log_file=log_file,
        report_file=report_file
    )
    
    logger.info("="*70)
    logger.info("BIOINFORMATICS VALIDATION PIPELINE")
    logger.info("="*70)
    logger.info(f"Configuration file: {config_path}")
    logger.info(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("="*70)
    
    try:
        # Step 1: Load and validate configuration
        logger.info("\n[1/4] Loading configuration...")
        config = Coordinator.load(config_path)
        logger.info("✓ Configuration loaded successfully")
        
        # Step 2: Validate genome files
        logger.info("\n[2/4] Validating genome files...")
        logger.info(f"  - Reference genome: {config.ref_genome.filename}")
        logger.info(f"  - Modified genome: {config.mod_genome.filename}")
        
        if config.ref_plasmid:
            logger.info(f"  - Reference plasmid: {config.ref_plasmid.filename}")
        if config.mod_plasmid:
            logger.info(f"  - Modified plasmid: {config.mod_plasmid.filename}")
        
        # Import and use GenomeValidator
        from validators.genome_validator import GenomeValidator
        
        output_dir = config.output_dir
        
        # Validate reference genome
        ref_genome_validator = GenomeValidator(
            config.ref_genome,
            config.config_dir,
            output_dir
        )
        ref_genome_validator.validate()
        
        # Validate modified genome
        mod_genome_validator = GenomeValidator(
            config.mod_genome,
            config.config_dir,
            output_dir
        )
        mod_genome_validator.validate()
        
        # Validate plasmids if present
        if config.ref_plasmid:
            ref_plasmid_validator = GenomeValidator(
                config.ref_plasmid,
                config.config_dir,
                output_dir / "plasmids"
            )
            ref_plasmid_validator.validate()
        
        if config.mod_plasmid:
            mod_plasmid_validator = GenomeValidator(
                config.mod_plasmid,
                config.config_dir,
                output_dir / "plasmids"
            )
            mod_plasmid_validator.validate()
        
        logger.info("✓ Genome validation completed")
        
        # Step 3: Validate feature files
        logger.info("\n[3/4] Validating feature files...")
        
        if config.ref_feature:
            logger.info(f"  - Reference features: {config.ref_feature.filename}")
        if config.mod_feature:
            logger.info(f"  - Modified features: {config.mod_feature.filename}")
        
        if not config.ref_feature and not config.mod_feature:
            logger.info("  - No feature files specified (skipping)")
        else:
            # TODO: Add feature validation here
            # from validators.feature_validator import FeatureValidator
            # if config.ref_feature:
            #     feature_validator = FeatureValidator(config.ref_feature)
            #     feature_validator.validate()
            
            logger.info("✓ Feature validation completed")
        
        # Step 4: Validate read files
        logger.info("\n[4/4] Validating read files...")
        logger.info(f"  - Number of read entries: {len(config.reads)}")
        
        for idx, read in enumerate(config.reads):
            if read.filename:
                logger.info(f"  - Read {idx+1}: {read.filename} ({read.ngs_type})")
            else:
                logger.info(f"  - Read {idx+1}: {read.directory}/ ({read.ngs_type})")
        
        # TODO: Add read validation here
        # from validators.read_validator import ReadValidator
        # for read in config.reads:
        #     read_validator = ReadValidator(read)
        #     read_validator.validate()
        
        logger.info("✓ Read validation completed")
        
        # Step 5: Inter-file validation
        logger.info("\n[5/5] Performing inter-file validation...")
        # TODO: Add inter-file validation here
        # - Check if features are within genome bounds
        # - Check sequence IDs match
        # etc.
        logger.info("✓ Inter-file validation completed")
        
        # Generate final report
        logger.info("\n" + "="*70)
        logger.info("VALIDATION COMPLETE")
        logger.info("="*70)
        
        summary = logger.get_summary()
        
        if summary['passed']:
            logger.info("✓ All validations passed successfully!")
            logger.info(f"  Total issues: {summary['total_issues']}")
            logger.info(f"  Errors: {summary['errors']}")
            logger.info(f"  Warnings: {summary['warnings']}")
        else:
            logger.error("✗ Validation failed with errors")
            logger.error(f"  Total issues: {summary['total_issues']}")
            logger.error(f"  Errors: {summary['errors']}")
            logger.error(f"  Warnings: {summary['warnings']}")
        
        if not args.no_report:
            logger.info("\n" + "="*70)
            logger.generate_report()
        
        logger.info("="*70)
        logger.info(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info("="*70)
        
        # Return exit code
        return 0 if summary['passed'] else 1
        
    except ConfigurationError as e:
        logger.error("\n" + "="*70)
        logger.error("CONFIGURATION ERROR")
        logger.error("="*70)
        logger.error(str(e))
        logger.error("\nPlease check your configuration file and try again.")
        
        if not args.no_report:
            logger.generate_report()
        
        return 1
        
    except FileNotFoundError as e:
        logger.error("\n" + "="*70)
        logger.error("FILE NOT FOUND ERROR")
        logger.error("="*70)
        logger.error(str(e))
        logger.error("\nPlease ensure all referenced files exist.")
        
        if not args.no_report:
            logger.generate_report()
        
        return 1
        
    except ValidationError as e:
        logger.error("\n" + "="*70)
        logger.error("VALIDATION ERROR")
        logger.error("="*70)
        logger.error(str(e))
        
        if not args.no_report:
            logger.generate_report()
        
        return 1
        
    except Exception as e:
        logger.critical("\n" + "="*70)
        logger.critical("UNEXPECTED ERROR")
        logger.critical("="*70)
        logger.critical(f"An unexpected error occurred: {type(e).__name__}: {e}")
        logger.critical("\nPlease report this issue to the developers.")
        
        if not args.no_report:
            logger.generate_report()
        
        # Print full traceback in debug mode
        if args.verbose:
            import traceback
            logger.critical("\nFull traceback:")
            logger.critical(traceback.format_exc())
        
        return 1


def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Check if config file exists
    if not Path(args.config).exists():
        print(f"Error: Configuration file not found: {args.config}")
        print("\nUsage: python main.py <config.json>")
        print("Run 'python main.py --help' for more options.")
        return 1
    
    # Run validation
    exit_code = run_validation(args.config, args)
    
    return exit_code


if __name__ == "__main__":
    sys.exit(main())