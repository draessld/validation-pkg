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

Quick Start
-----------

**Using the Functional API:**

>>> from validation_pkg import ConfigManager, validate_genome, validate_read, validate_feature
>>> config = ConfigManager.load("config.json")
>>>
>>> # Validate genome files
>>> if config.ref_genome:
...     validate_genome(config.ref_genome, config.output_dir)
>>>
>>> # Validate read files
>>> if config.reads:
...     for read in config.reads:
...         validate_read(read, config.output_dir)
>>>
>>> # Validate feature files
>>> if config.ref_feature:
...     validate_feature(config.ref_feature, config.output_dir)

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
from validation_pkg.config_manager import ConfigManager, Config
from validation_pkg.validators.genome_validator import GenomeValidator
from validation_pkg.validators.read_validator import ReadValidator
from validation_pkg.validators.feature_validator import FeatureValidator
from validation_pkg.logger import setup_logging, get_logger

# Functional API imports
from pathlib import Path
from typing import Optional, List, Union
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing


# ============================================================================
# Functional API - Simplified wrapper functions
# ============================================================================

# ============================================================================
# Parallel Processing Helper Functions
# ============================================================================

def _validate_single_file(args):
    """
    Helper function for parallel validation of a single file.

    This function is designed to be called by ProcessPoolExecutor.
    It must be a module-level function (not nested) to be picklable.

    Args:
        args: Tuple of (validator_type, config, output_dir, settings_dict, config_threads, worker_id)

    Returns:
        dict: Result with 'success', 'filename', 'error' keys
    """
    validator_type, config, output_dir, settings_dict, config_threads, worker_id = args

    try:
        # Import validators and logger inside worker process
        from validation_pkg.validators.genome_validator import GenomeValidator
        from validation_pkg.validators.read_validator import ReadValidator
        from validation_pkg.validators.feature_validator import FeatureValidator
        from validation_pkg.logger import get_logger

        # Bind worker context to logger for this process
        logger = get_logger()
        filename_short = Path(config.filename).name if hasattr(config, 'filename') else 'unknown'
        logger.bind_worker_context(worker_id=worker_id, file_context=filename_short)

        # Reconstruct settings from dict
        if validator_type == 'genome':
            if settings_dict:
                settings = GenomeValidator.Settings.from_dict(settings_dict)
            else:
                settings = GenomeValidator.Settings()
            validator = GenomeValidator(config, Path(output_dir), settings)
        elif validator_type == 'read':
            if settings_dict:
                settings = ReadValidator.Settings.from_dict(settings_dict)
            else:
                settings = ReadValidator.Settings()
            validator = ReadValidator(config, Path(output_dir), settings)
        elif validator_type == 'feature':
            if settings_dict:
                settings = FeatureValidator.Settings.from_dict(settings_dict)
            else:
                settings = FeatureValidator.Settings()
            validator = FeatureValidator(config, Path(output_dir), settings)
        else:
            raise ValueError(f"Unknown validator type: {validator_type}")

        # Run validation (logs will now include worker context)
        validator.validate()

        return {
            'success': True,
            'filename': config.filename,
            'error': None
        }

    except Exception as e:
        return {
            'success': False,
            'filename': config.filename if hasattr(config, 'filename') else 'unknown',
            'error': str(e)
        }

def validate_genome(
    genome_config,
    output_dir: Union[str, Path],
    settings: Optional[GenomeValidator.Settings] = None,
    config_threads: Optional[int] = None
) -> dict:
    """
    Validate a genome file with optional custom settings.

    This is a simplified wrapper around GenomeValidator for easier usage.

    Args:
        genome_config: GenomeConfig object (from ConfigManager)
        output_dir: Directory for output files
        settings: Optional GenomeValidator.Settings object (uses defaults if None)
        config_threads: Thread count from config.get_threads() (for internal use)

    Example:
        >>> config = ConfigManager.load("config.json")
        >>> settings = GenomeValidator.Settings()
        >>> settings = settings.update(coding_type='gz', plasmid_split=True)
        >>> stats = validate_genome(config.ref_genome, config.output_dir, settings)
        >>> print(f"Validated {stats['total_sequences']} sequences")
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Apply settings from config if available
    if settings is None:
        # No user settings provided - use config settings if available
        if genome_config.settings_dict:
            settings = GenomeValidator.Settings.from_dict(genome_config.settings_dict)
        else:
            settings = GenomeValidator.Settings()
    else:
        # User settings provided - merge with config settings (user takes precedence)
        if genome_config.settings_dict:
            # Start with config settings, then apply user settings on top
            merged_settings = genome_config.settings_dict.copy()
            merged_settings.update(settings.to_dict())
            settings = GenomeValidator.Settings.from_dict(merged_settings)

    # If threads specified in config but not in settings, apply it
    if config_threads is not None:
        if settings.compression_threads is None:
            settings = settings.update(compression_threads=config_threads)

    validator = GenomeValidator(genome_config, output_path, settings)
    validator.validate()


def validate_read(
    read_config,
    output_dir: Union[str, Path],
    settings: Optional[ReadValidator.Settings] = None,
    config_threads: Optional[int] = None
) -> dict:
    """
    Validate a single read file with optional custom settings.

    This is a simplified wrapper around ReadValidator for easier usage.

    Args:
        read_config: ReadConfig object (from ConfigManager)
        output_dir: Directory for output files
        settings: Optional ReadValidator.Settings object (uses defaults if None)
        config_threads: Thread count from config.get_threads() (for internal use)


    Example:
        >>> config = ConfigManager.load("config.json")
        >>> settings = ReadValidator.Settings()
        >>> settings = settings.update(coding_type='gz')
        >>> stats = validate_read(config.reads[0], config.output_dir, settings)
        >>> print(f"Validated {stats['total_reads']} reads")
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Apply settings from config if available
    if settings is None:
        # No user settings provided - use config settings if available
        if read_config.settings_dict:
            settings = ReadValidator.Settings.from_dict(read_config.settings_dict)
        else:
            settings = ReadValidator.Settings()
    else:
        # User settings provided - merge with config settings (user takes precedence)
        if read_config.settings_dict:
            # Start with config settings, then apply user settings on top
            merged_settings = read_config.settings_dict.copy()
            merged_settings.update(settings.to_dict())
            settings = ReadValidator.Settings.from_dict(merged_settings)

    # If threads specified in config but not in settings, apply it
    if config_threads is not None:
        if settings.compression_threads is None:
            settings = settings.update(compression_threads=config_threads)

    validator = ReadValidator(read_config, output_path, settings)
    validator.validate()


def validate_reads(
    read_configs: List,
    output_dir: Union[str, Path],
    settings: Optional[ReadValidator.Settings] = None,
    config_threads: Optional[int] = None
) -> List[dict]:
    """
    Validate multiple read files with optional custom settings.

    Supports parallel processing if settings.max_workers is set.

    Args:
        read_configs: List of ReadConfig objects (from ConfigManager)
        output_dir: Directory for output files
        settings: Optional ReadValidator.Settings object (uses defaults if None)
        config_threads: Thread count from config.get_threads() (for internal use)

    Returns:
        List of validation results (dicts with 'success', 'filename', 'error' keys)

    Example:
        >>> config = ConfigManager.load("config.json")
        >>> settings = ReadValidator.Settings()
        >>> settings = settings.update(coding_type='gz', max_workers=4)
        >>> results = validate_reads(config.reads, config.output_dir, settings)
        >>> for result in results:
        ...     if result['success']:
        ...         print(f"✓ {result['filename']}")
        ...     else:
        ...         print(f"✗ {result['filename']}: {result['error']}")
    """
    logger = get_logger()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Prepare settings for each file
    tasks = []
    for read_config in read_configs:
        # Each read file gets its own settings (config-specific + user-provided)
        file_settings = None

        if settings is None:
            # No user settings - use per-file config settings if available
            if read_config.settings_dict:
                file_settings = ReadValidator.Settings.from_dict(read_config.settings_dict)
            else:
                file_settings = ReadValidator.Settings()
        else:
            # User settings provided - merge with per-file config settings
            if read_config.settings_dict:
                merged_settings = read_config.settings_dict.copy()
                merged_settings.update(settings.to_dict())
                file_settings = ReadValidator.Settings.from_dict(merged_settings)
            else:
                file_settings = settings.copy()

        # Apply threads from config
        if config_threads is not None:
            if file_settings.compression_threads is None:
                file_settings = file_settings.update(compression_threads=config_threads)

        # Set output subdirectory by ngs_type
        file_settings = file_settings.update(output_subdir_name=read_config.ngs_type)

        tasks.append((read_config, file_settings))

    # Apply smart thread splitting if unified 'threads' parameter is provided
    # Only apply if explicit max_workers/compression_threads are not set
    if settings and settings.threads is not None:
        if settings.max_workers is None and settings.compression_threads is None:
            # Import here to avoid circular dependency
            from validation_pkg.utils.file_handler import calculate_thread_distribution

            # Calculate optimal distribution
            max_workers_calc, compression_threads_calc = calculate_thread_distribution(
                settings.threads, len(tasks)
            )

            logger.debug(
                f"Smart thread splitting: {settings.threads} threads → "
                f"{max_workers_calc} workers × {compression_threads_calc} compression threads"
            )

            # Apply to all task settings
            updated_tasks = []
            for read_config, file_settings in tasks:
                file_settings = file_settings.update(
                    max_workers=max_workers_calc,
                    compression_threads=compression_threads_calc
                )
                updated_tasks.append((read_config, file_settings))
            tasks = updated_tasks

            # Use calculated max_workers for parallel execution decision
            max_workers = max_workers_calc
        else:
            # Explicit settings provided, use them
            max_workers = settings.max_workers
    else:
        # Check if parallel processing is requested (old behavior)
        max_workers = settings.max_workers if settings else None

    if max_workers and max_workers > 1 and len(tasks) > 1:
        # Parallel execution
        logger.info(f"Processing {len(tasks)} read files in parallel with {max_workers} workers")

        results = []
        args_list = [
            ('read', read_config, str(output_path), file_settings.to_dict(), config_threads, worker_id)
            for worker_id, (read_config, file_settings) in enumerate(tasks, start=1)
        ]

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_filename = {
                executor.submit(_validate_single_file, args): args[1].filename
                for args in args_list
            }

            # Collect results as they complete
            for future in as_completed(future_to_filename):
                filename = future_to_filename[future]
                try:
                    result = future.result()
                    results.append(result)

                    if result['success']:
                        logger.info(f"✓ Completed: {filename}")
                    else:
                        logger.error(f"✗ Failed: {filename} - {result['error']}")

                except Exception as e:
                    error_result = {
                        'success': False,
                        'filename': filename,
                        'error': str(e)
                    }
                    results.append(error_result)
                    logger.error(f"✗ Exception: {filename} - {e}")

        # Check for failures
        failures = [r for r in results if not r['success']]
        if failures:
            logger.error(f"Validation failed for {len(failures)}/{len(tasks)} files")
            for failure in failures:
                logger.error(f"  - {failure['filename']}: {failure['error']}")
        else:
            logger.info(f"✓ All {len(tasks)} read files validated successfully")

        return results

    else:
        # Sequential execution (original behavior)
        if max_workers:
            logger.debug(f"Sequential processing (max_workers={max_workers} but only {len(tasks)} file(s))")
        else:
            logger.debug("Sequential processing (max_workers not set)")

        results = []
        for read_config, file_settings in tasks:
            try:
                validator = ReadValidator(read_config, output_path, file_settings)
                validator.validate()
                results.append({
                    'success': True,
                    'filename': read_config.filename,
                    'error': None
                })
            except Exception as e:
                logger.error(f"Validation failed for {read_config.filename}: {e}")
                results.append({
                    'success': False,
                    'filename': read_config.filename,
                    'error': str(e)
                })

        return results


def validate_genomes(
    genome_configs: List,
    output_dir: Union[str, Path],
    settings: Optional[GenomeValidator.Settings] = None,
    config_threads: Optional[int] = None
) -> List[dict]:
    """
    Validate multiple genome files with optional custom settings.

    Supports parallel processing if settings.max_workers is set.

    Args:
        genome_configs: List of GenomeConfig objects (from ConfigManager)
        output_dir: Directory for output files
        settings: Optional GenomeValidator.Settings object (uses defaults if None)
        config_threads: Thread count from config.get_threads() (for internal use)

    Returns:
        List of validation results (dicts with 'success', 'filename', 'error' keys)

    Example:
        >>> config = ConfigManager.load("config.json")
        >>> settings = GenomeValidator.Settings()
        >>> settings = settings.update(coding_type='gz', max_workers=2)
        >>> genome_list = [config.ref_genome, config.mod_genome]
        >>> results = validate_genomes(genome_list, config.output_dir, settings)
    """
    logger = get_logger()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Prepare settings for each file
    tasks = []
    for genome_config in genome_configs:
        file_settings = None

        if settings is None:
            if genome_config.settings_dict:
                file_settings = GenomeValidator.Settings.from_dict(genome_config.settings_dict)
            else:
                file_settings = GenomeValidator.Settings()
        else:
            if genome_config.settings_dict:
                merged_settings = genome_config.settings_dict.copy()
                merged_settings.update(settings.to_dict())
                file_settings = GenomeValidator.Settings.from_dict(merged_settings)
            else:
                file_settings = settings.copy()

        if config_threads is not None:
            if file_settings.compression_threads is None:
                file_settings = file_settings.update(compression_threads=config_threads)

        tasks.append((genome_config, file_settings))

    # Apply smart thread splitting if unified 'threads' parameter is provided
    # Only apply if explicit max_workers/compression_threads are not set
    if settings and settings.threads is not None:
        if settings.max_workers is None and settings.compression_threads is None:
            # Import here to avoid circular dependency
            from validation_pkg.utils.file_handler import calculate_thread_distribution

            # Calculate optimal distribution
            max_workers_calc, compression_threads_calc = calculate_thread_distribution(
                settings.threads, len(tasks)
            )

            logger.debug(
                f"Smart thread splitting: {settings.threads} threads → "
                f"{max_workers_calc} workers × {compression_threads_calc} compression threads"
            )

            # Apply to all task settings
            updated_tasks = []
            for genome_config, file_settings in tasks:
                file_settings = file_settings.update(
                    max_workers=max_workers_calc,
                    compression_threads=compression_threads_calc
                )
                updated_tasks.append((genome_config, file_settings))
            tasks = updated_tasks

            # Use calculated max_workers for parallel execution decision
            max_workers = max_workers_calc
        else:
            # Explicit settings provided, use them
            max_workers = settings.max_workers
    else:
        # Check if parallel processing is requested (old behavior)
        max_workers = settings.max_workers if settings else None

    if max_workers and max_workers > 1 and len(tasks) > 1:
        logger.info(f"Processing {len(tasks)} genome files in parallel with {max_workers} workers")

        results = []
        args_list = [
            ('genome', genome_config, str(output_path), file_settings.to_dict(), config_threads, worker_id)
            for worker_id, (genome_config, file_settings) in enumerate(tasks, start=1)
        ]

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_filename = {
                executor.submit(_validate_single_file, args): args[1].filename
                for args in args_list
            }

            for future in as_completed(future_to_filename):
                filename = future_to_filename[future]
                try:
                    result = future.result()
                    results.append(result)

                    if result['success']:
                        logger.info(f"✓ Completed: {filename}")
                    else:
                        logger.error(f"✗ Failed: {filename} - {result['error']}")

                except Exception as e:
                    error_result = {
                        'success': False,
                        'filename': filename,
                        'error': str(e)
                    }
                    results.append(error_result)
                    logger.error(f"✗ Exception: {filename} - {e}")

        failures = [r for r in results if not r['success']]
        if failures:
            logger.error(f"Validation failed for {len(failures)}/{len(tasks)} files")
            for failure in failures:
                logger.error(f"  - {failure['filename']}: {failure['error']}")
        else:
            logger.info(f"✓ All {len(tasks)} genome files validated successfully")

        return results

    else:
        if max_workers:
            logger.debug(f"Sequential processing (max_workers={max_workers} but only {len(tasks)} file(s))")
        else:
            logger.debug("Sequential processing (max_workers not set)")

        results = []
        for genome_config, file_settings in tasks:
            try:
                validator = GenomeValidator(genome_config, output_path, file_settings)
                validator.validate()
                results.append({
                    'success': True,
                    'filename': genome_config.filename,
                    'error': None
                })
            except Exception as e:
                logger.error(f"Validation failed for {genome_config.filename}: {e}")
                results.append({
                    'success': False,
                    'filename': genome_config.filename,
                    'error': str(e)
                })

        return results


def validate_feature(
    feature_config,
    output_dir: Union[str, Path],
    settings: Optional[FeatureValidator.Settings] = None,
    config_threads: Optional[int] = None
) -> dict:
    """
    Validate a feature annotation file with optional custom settings.

    This is a simplified wrapper around FeatureValidator for easier usage.

    Args:
        feature_config: FeatureConfig object (from ConfigManager)
        output_dir: Directory for output files
        settings: Optional FeatureValidator.Settings object (uses defaults if None)
        config_threads: Thread count from config.get_threads() (for internal use)


    Example:
        >>> config = ConfigManager.load("config.json")
        >>> settings = FeatureValidator.Settings()
        >>> settings = settings.update(coding_type='gz', sort_by_position=True)
        >>> stats = validate_feature(config.ref_feature, config.output_dir, settings)
        >>> print(f"Validated {stats['total_features']} features")
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Apply settings from config if available
    if settings is None:
        # No user settings provided - use config settings if available
        if feature_config.settings_dict:
            settings = FeatureValidator.Settings.from_dict(feature_config.settings_dict)
        else:
            settings = FeatureValidator.Settings()
    else:
        # User settings provided - merge with config settings (user takes precedence)
        if feature_config.settings_dict:
            # Start with config settings, then apply user settings on top
            merged_settings = feature_config.settings_dict.copy()
            merged_settings.update(settings.to_dict())
            settings = FeatureValidator.Settings.from_dict(merged_settings)

    # If threads specified in config but not in settings, apply it
    if config_threads is not None:
        if settings.compression_threads is None:
            settings = settings.update(compression_threads=config_threads)

    validator = FeatureValidator(feature_config, output_path, settings)
    validator.validate()


def validate_features_list(
    feature_configs: List,
    output_dir: Union[str, Path],
    settings: Optional[FeatureValidator.Settings] = None,
    config_threads: Optional[int] = None
) -> List[dict]:
    """
    Validate multiple feature annotation files with optional custom settings.

    Supports parallel processing if settings.max_workers is set.

    Args:
        feature_configs: List of FeatureConfig objects (from ConfigManager)
        output_dir: Directory for output files
        settings: Optional FeatureValidator.Settings object (uses defaults if None)
        config_threads: Thread count from config.get_threads() (for internal use)

    Returns:
        List of validation results (dicts with 'success', 'filename', 'error' keys)

    Example:
        >>> config = ConfigManager.load("config.json")
        >>> settings = FeatureValidator.Settings()
        >>> settings = settings.update(coding_type='gz', max_workers=2)
        >>> feature_list = [config.ref_feature, config.mod_feature]
        >>> results = validate_features_list(feature_list, config.output_dir, settings)
    """
    logger = get_logger()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Prepare settings for each file
    tasks = []
    for feature_config in feature_configs:
        file_settings = None

        if settings is None:
            if feature_config.settings_dict:
                file_settings = FeatureValidator.Settings.from_dict(feature_config.settings_dict)
            else:
                file_settings = FeatureValidator.Settings()
        else:
            if feature_config.settings_dict:
                merged_settings = feature_config.settings_dict.copy()
                merged_settings.update(settings.to_dict())
                file_settings = FeatureValidator.Settings.from_dict(merged_settings)
            else:
                file_settings = settings.copy()

        if config_threads is not None:
            if file_settings.compression_threads is None:
                file_settings = file_settings.update(compression_threads=config_threads)

        tasks.append((feature_config, file_settings))

    # Apply smart thread splitting if unified 'threads' parameter is provided
    # Only apply if explicit max_workers/compression_threads are not set
    if settings and settings.threads is not None:
        if settings.max_workers is None and settings.compression_threads is None:
            # Import here to avoid circular dependency
            from validation_pkg.utils.file_handler import calculate_thread_distribution

            # Calculate optimal distribution
            max_workers_calc, compression_threads_calc = calculate_thread_distribution(
                settings.threads, len(tasks)
            )

            logger.debug(
                f"Smart thread splitting: {settings.threads} threads → "
                f"{max_workers_calc} workers × {compression_threads_calc} compression threads"
            )

            # Apply to all task settings
            updated_tasks = []
            for feature_config, file_settings in tasks:
                file_settings = file_settings.update(
                    max_workers=max_workers_calc,
                    compression_threads=compression_threads_calc
                )
                updated_tasks.append((feature_config, file_settings))
            tasks = updated_tasks

            # Use calculated max_workers for parallel execution decision
            max_workers = max_workers_calc
        else:
            # Explicit settings provided, use them
            max_workers = settings.max_workers
    else:
        # Check if parallel processing is requested (old behavior)
        max_workers = settings.max_workers if settings else None

    if max_workers and max_workers > 1 and len(tasks) > 1:
        logger.info(f"Processing {len(tasks)} feature files in parallel with {max_workers} workers")

        results = []
        args_list = [
            ('feature', feature_config, str(output_path), file_settings.to_dict(), config_threads, worker_id)
            for worker_id, (feature_config, file_settings) in enumerate(tasks, start=1)
        ]

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_filename = {
                executor.submit(_validate_single_file, args): args[1].filename
                for args in args_list
            }

            for future in as_completed(future_to_filename):
                filename = future_to_filename[future]
                try:
                    result = future.result()
                    results.append(result)

                    if result['success']:
                        logger.info(f"✓ Completed: {filename}")
                    else:
                        logger.error(f"✗ Failed: {filename} - {result['error']}")

                except Exception as e:
                    error_result = {
                        'success': False,
                        'filename': filename,
                        'error': str(e)
                    }
                    results.append(error_result)
                    logger.error(f"✗ Exception: {filename} - {e}")

        failures = [r for r in results if not r['success']]
        if failures:
            logger.error(f"Validation failed for {len(failures)}/{len(tasks)} files")
            for failure in failures:
                logger.error(f"  - {failure['filename']}: {failure['error']}")
        else:
            logger.info(f"✓ All {len(tasks)} feature files validated successfully")

        return results

    else:
        if max_workers:
            logger.debug(f"Sequential processing (max_workers={max_workers} but only {len(tasks)} file(s))")
        else:
            logger.debug("Sequential processing (max_workers not set)")

        results = []
        for feature_config, file_settings in tasks:
            try:
                validator = FeatureValidator(feature_config, output_path, file_settings)
                validator.validate()
                results.append({
                    'success': True,
                    'filename': feature_config.filename,
                    'error': None
                })
            except Exception as e:
                logger.error(f"Validation failed for {feature_config.filename}: {e}")
                results.append({
                    'success': False,
                    'filename': feature_config.filename,
                    'error': str(e)
                })

        return results

# Alias for backward compatibility and convenience
validate_features = validate_feature


__all__ = [
    # Configuration
    'ConfigManager',
    'Config',

    # Validators
    'GenomeValidator',
    'ReadValidator',
    'FeatureValidator',

    # Functional API (Primary Interface)
    'validate_genome',
    'validate_genomes',  # New: parallel genome validation
    'validate_read',
    'validate_reads',
    'validate_feature',
    'validate_features',
    'validate_features_list',  # New: parallel feature validation

    # Logging
    'setup_logging',
    'get_logger',

    # Version info
    '__version__',
    '__author__',
    '__license__',
]
