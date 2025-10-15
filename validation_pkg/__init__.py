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

    # Logging
    'setup_logging',
    'get_logger',

    # Version info
    '__version__',
    '__author__',
    '__license__',
]
