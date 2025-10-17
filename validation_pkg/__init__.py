"""
Bioinformatics Validation Package
==================================

A comprehensive validation package for genomic data files used in bioinformatics
pipelines. Validates and standardizes genome files, sequencing reads, and
feature annotations with detailed error reporting and quality checks.

Supported File Types
-------------------
- **Genome files:** FASTA (.fasta, .fa, .fna), GenBank (.gb, .gbk, .genbank)
- **Read files:** FASTQ (.fastq, .fq), BAM (.bam)
- **Feature files:** GFF (.gff, .gff3, .gtf), BED (.bed)
- **Compression:** gzip (.gz), bzip2 (.bz2)

Features
--------
- Format detection and validation
- Automatic compression/decompression
- Sequence quality checks (duplicate IDs, invalid characters, empty sequences)
- Coordinate validation for features
- Format conversion (e.g., GenBank to FASTA, BED to GFF)
- Plasmid splitting for bacterial genomes
- BAM to FASTQ conversion
- Detailed validation reports with structured error logging
- Flexible configuration via JSON
- Command-line interface

Quick Start
-----------

**Using the API:**

>>> from validation_pkg import ValidationCoordinator
>>> coordinator = ValidationCoordinator("config.json")
>>> report = coordinator.validate_all()
>>> print(report.summary())

**Using the CLI:**

.. code-block:: bash

    # Validate all files in config
    python -m validation_pkg validate config.json

    # Validate only genomes
    python -m validation_pkg validate config.json --only genomes

    # Enable verbose output
    python -m validation_pkg validate config.json --verbose

Configuration Example
--------------------

Create a config.json file:

.. code-block:: json

    {
      "output_dir": "./output",
      "ref_genome_filename": {"filename": "reference.fasta"},
      "mod_genome_filename": {"filename": "modified.fasta"},
      "reads": [
        {"filename": "sample1.fastq.gz", "ngs_type": "illumina"},
        {"filename": "sample2.fastq.gz", "ngs_type": "illumina"}
      ],
      "ref_feature_filename": {"filename": "annotations.gff"}
    }

See docs/config_guide.md for detailed configuration options.

Package Structure
----------------
- validation_pkg.config_manager: Configuration loading and parsing
- validation_pkg.coordinator: Workflow orchestration
- validation_pkg.validators: Individual file validators
    - genome_validator: Genome file validation
    - read_validator: Read file validation
    - feature_validator: Feature annotation validation
- validation_pkg.utils: Utility functions and formats
- validation_pkg.logger: Structured logging system

Workflow
--------
1. **Load configuration:** Parse JSON config and detect file formats
2. **Validate genomes:** Check reference and modified genome files
3. **Validate features:** Validate feature annotations (optional)
4. **Validate reads:** Process all sequencing read files
5. **Generate report:** Create detailed validation report with statistics

Error Handling
-------------
The package uses a hierarchical exception system:

- ValidationError: Base exception for all validation errors
    - GenomeValidationError: Genome-specific errors
        - FastaFormatError
        - GenBankFormatError
    - ReadValidationError: Read-specific errors
        - FastqFormatError
        - BamFormatError
    - FeatureValidationError: Feature-specific errors
    - ConfigurationError: Config file errors
    - FileNotFoundError: Missing file errors

All errors are logged with structured details for debugging.

Version Information
------------------
"""

__version__ = "0.1.0"
__author__ = "Dominika Bohuslavova"
__license__ = "MIT"

# Public API exports
from validation_pkg.coordinator import ValidationCoordinator, Coordinator, ValidationReport
from validation_pkg.config_manager import ConfigManager, Config
from validation_pkg.validators.genome_validator import GenomeValidator
from validation_pkg.validators.read_validator import ReadValidator
from validation_pkg.validators.feature_validator import FeatureValidator
from validation_pkg.logger import setup_logging, get_logger

# Functional API imports
from pathlib import Path
from typing import Optional, List, Union


# ============================================================================
# Functional API - Simplified wrapper functions
# ============================================================================

def validate_genome(
    genome_config,
    output_dir: Union[str, Path],
    settings: Optional[GenomeValidator.Settings] = None
) -> dict:
    """
    Validate a genome file with optional custom settings.

    This is a simplified wrapper around GenomeValidator for easier usage.

    Args:
        genome_config: GenomeConfig object (from ConfigManager)
        output_dir: Directory for output files
        settings: Optional GenomeValidator.Settings object (uses defaults if None)

    Returns:
        Dictionary with validation statistics

    Example:
        >>> config = ConfigManager.load("config.json")
        >>> settings = GenomeValidator.Settings()
        >>> settings = settings.update(coding_type='gz', plasmid_split=True)
        >>> stats = validate_genome(config.ref_genome, config.output_dir, settings)
        >>> print(f"Validated {stats['total_sequences']} sequences")
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    validator = GenomeValidator(genome_config, output_path, settings)
    validator.validate()

    return validator.get_statistics()


def validate_read(
    read_config,
    output_dir: Union[str, Path],
    settings: Optional[ReadValidator.Settings] = None
) -> dict:
    """
    Validate a single read file with optional custom settings.

    This is a simplified wrapper around ReadValidator for easier usage.

    Args:
        read_config: ReadConfig object (from ConfigManager)
        output_dir: Directory for output files
        settings: Optional ReadValidator.Settings object (uses defaults if None)

    Returns:
        Dictionary with validation statistics

    Example:
        >>> config = ConfigManager.load("config.json")
        >>> settings = ReadValidator.Settings()
        >>> settings = settings.update(coding_type='gz')
        >>> stats = validate_read(config.reads[0], config.output_dir, settings)
        >>> print(f"Validated {stats['total_reads']} reads")
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    validator = ReadValidator(read_config, output_path, settings)
    validator.validate()

    return validator.get_statistics()


def validate_reads(
    read_configs: List,
    output_dir: Union[str, Path],
    settings: Optional[ReadValidator.Settings] = None
) -> List[dict]:
    """
    Validate multiple read files with optional custom settings.

    This is a simplified wrapper for validating all reads in a config.

    Args:
        read_configs: List of ReadConfig objects (from ConfigManager)
        output_dir: Directory for output files
        settings: Optional ReadValidator.Settings object (uses defaults if None)

    Returns:
        List of dictionaries with validation statistics for each file

    Example:
        >>> config = ConfigManager.load("config.json")
        >>> settings = ReadValidator.Settings()
        >>> settings = settings.update(coding_type='gz')
        >>> stats_list = validate_reads(config.reads, config.output_dir, settings)
        >>> for idx, stats in enumerate(stats_list, 1):
        ...     print(f"Read file {idx}: {stats['total_reads']} reads")
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    results = []
    for read_config in read_configs:
        validator = ReadValidator(read_config, output_path, settings)
        validator.validate()
        results.append(validator.get_statistics())

    return results


def validate_feature(
    feature_config,
    output_dir: Union[str, Path],
    settings: Optional[FeatureValidator.Settings] = None
) -> dict:
    """
    Validate a feature annotation file with optional custom settings.

    This is a simplified wrapper around FeatureValidator for easier usage.

    Args:
        feature_config: FeatureConfig object (from ConfigManager)
        output_dir: Directory for output files
        settings: Optional FeatureValidator.Settings object (uses defaults if None)

    Returns:
        Dictionary with validation statistics

    Example:
        >>> config = ConfigManager.load("config.json")
        >>> settings = FeatureValidator.Settings()
        >>> settings = settings.update(coding_type='gz', sort_by_position=True)
        >>> stats = validate_feature(config.ref_feature, config.output_dir, settings)
        >>> print(f"Validated {stats['total_features']} features")
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    validator = FeatureValidator(feature_config, output_path, settings)
    validator.validate()

    return validator.get_statistics()


# Alias for backward compatibility and convenience
validate_features = validate_feature


__all__ = [
    # Main classes
    'ValidationCoordinator',
    'Coordinator',  # Backward compatibility alias
    'ValidationReport',
    'ConfigManager',
    'Config',

    # Validators
    'GenomeValidator',
    'ReadValidator',
    'FeatureValidator',

    # Functional API
    'validate_genome',
    'validate_read',
    'validate_reads',
    'validate_feature',
    'validate_features',

    # Logging
    'setup_logging',
    'get_logger',

    # Version info
    '__version__',
    '__author__',
    '__license__',
]
