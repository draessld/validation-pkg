"""
Configuration management for bioinformatics file validation.

This module provides classes and utilities for loading and validating JSON configuration
files that specify input files, output directories, and validation settings.

Main Components:
    - GenomeConfig: Configuration for genome/plasmid files
    - ReadConfig: Configuration for sequencing read files
    - FeatureConfig: Configuration for feature annotation files
    - Config: Main configuration container
    - ConfigManager: JSON parser and validator

Key Features:
    - Automatic path resolution relative to config file location
    - Compression and format detection
    - Directory-based read loading
    - Validator-specific settings extraction from config
    - Thread count configuration for parallel compression
"""

import json
from pathlib import Path
from typing import List, Optional, Dict, Any, Union, Type, Tuple
from dataclasses import dataclass

from validation_pkg.utils.formats import CodingType, GenomeFormat, ReadFormat, FeatureFormat
from validation_pkg.logger import get_logger
from validation_pkg.exceptions import (
    ConfigurationError,
    FileNotFoundError as ValidationFileNotFoundError
)

@dataclass
class GenomeConfig:
    """
    Configuration for genome or plasmid files.

    Attributes:
        filename: Original filename (relative to config)
        filepath: Absolute resolved path
        coding_type: Compression format (GZIP, BZIP2, or NONE)
        detected_format: File format (FASTA or GENBANK)
        output_dir: Base output directory from config
        extra: Additional non-validator fields from config
        settings_dict: Validator-specific settings from config
    """
    filename: str
    filepath: Path
    coding_type: CodingType = None
    detected_format: GenomeFormat = None
    output_dir: Path = None
    extra: Dict[str, Any] = None
    settings_dict: Dict[str, Any] = None

    def __post_init__(self):
        if self.extra is None:
            self.extra = {}
        if self.settings_dict is None:
            self.settings_dict = {}

@dataclass
class ReadConfig:
    """
    Configuration for sequencing read files.

    Attributes:
        filename: Original filename (relative to config)
        filepath: Absolute resolved path
        ngs_type: Sequencing platform ("illumina", "ont", or "pacbio")
        coding_type: Compression format (GZIP, BZIP2, or NONE)
        detected_format: File format (FASTQ or BAM)
        output_dir: Base output directory from config
        extra: Additional non-validator fields from config
        settings_dict: Validator-specific settings from config
    """
    filename: str
    filepath: Path
    ngs_type: str = None
    coding_type: CodingType = None
    detected_format: ReadFormat = None
    output_dir: Path = None
    extra: Dict[str, Any] = None
    settings_dict: Dict[str, Any] = None

    def __post_init__(self):
        if self.ngs_type not in ["illumina", "ont", "pacbio"]:
            raise ValueError(f"Invalid ngs_type: {self.ngs_type}")
        if self.extra is None:
            self.extra = {}
        if self.settings_dict is None:
            self.settings_dict = {}

@dataclass
class FeatureConfig:
    """
    Configuration for feature annotation files.

    Attributes:
        filename: Original filename (relative to config)
        filepath: Absolute resolved path
        coding_type: Compression format (GZIP, BZIP2, or NONE)
        detected_format: File format (GFF, GTF, or BED)
        output_dir: Base output directory from config
        extra: Additional non-validator fields from config
        settings_dict: Validator-specific settings from config
    """
    filename: str
    filepath: Path
    coding_type: CodingType = None
    detected_format: FeatureFormat = None
    output_dir: Path = None
    extra: Dict[str, Any] = None
    settings_dict: Dict[str, Any] = None

    def __post_init__(self):
        if self.extra is None:
            self.extra = {}
        if self.settings_dict is None:
            self.settings_dict = {}


class Config:
    """
    Main configuration container for validation workflow.

    Stores validated configuration data loaded from JSON file, including file paths,
    compression/format details, and options. Provides helper methods for path resolution.

    Attributes:
        ref_genome: Reference genome configuration (required)
        mod_genome: Modified genome configuration (required)
        reads: List of read file configurations (required, non-empty)
        ref_plasmid: Reference plasmid configuration (optional)
        mod_plasmid: Modified plasmid configuration (optional)
        ref_feature: Reference feature configuration (optional)
        mod_feature: Modified feature configuration (optional)
        options: Additional options (e.g., thread count)
        config_dir: Directory containing config file
        output_dir: Base output directory for validated files
    """

    def __init__(self):
        # Required fields
        self.ref_genome: Optional[GenomeConfig] = None
        self.mod_genome: Optional[GenomeConfig] = None
        self.reads: List[ReadConfig] = []

        # Optional fields
        self.ref_plasmid: Optional[GenomeConfig] = None
        self.mod_plasmid: Optional[GenomeConfig] = None
        self.ref_feature: Optional[FeatureConfig] = None
        self.mod_feature: Optional[FeatureConfig] = None
        self.options: Dict[str, Any] = {}

        # Internal fields
        self.config_dir: Optional[Path] = None
        self.output_dir: Optional[Path] = None

    def get_threads(self) -> Optional[int]:
        """
        Get thread count from options, or None for auto-detection.

        Returns:
            Thread count if specified in options, None for auto-detection

        Example:
            >>> config = Config()
            >>> config.options = {'threads': 4}
            >>> config.get_threads()
            4
            >>> config.options = {}
            >>> config.get_threads()
            None
        """
        return self.options.get('threads')
    
    def resolve_path(self, relative_path: str) -> Path:
        """
        Resolve a path relative to config directory.
        
        Args:
            relative_path: Path relative to config file
            
        Returns:
            Absolute Path object
        """
        if self.config_dir is None:
            return Path(relative_path)
        return self.config_dir / relative_path
    
    def __repr__(self):
        return (
            f"Config(\n"
            f"  ref_genome={self.ref_genome},\n"
            f"  mod_genome={self.mod_genome},\n"
            f"  reads={len(self.reads)} file(s),\n"
            f"  ref_plasmid={self.ref_plasmid},\n"
            f"  mod_plasmid={self.mod_plasmid},\n"
            f"  ref_feature={self.ref_feature},\n"
            f"  mod_feature={self.mod_feature}\n"
            f")"
        )


class ConfigManager:
    """
    Configuration file parser and validator.

    Responsible for loading JSON configuration files, validating their structure,
    resolving file paths, detecting formats/compression, and extracting validator
    settings.

    Primary usage:
        config = ConfigManager.load("config.json")
    """
    
    @staticmethod
    def load(config_path: str) -> Config:
        """
        Load and validate configuration from JSON file.
        
        Args:
            config_path: Path to config.json
            
        Returns:
            Validated Config object
            
        Raises:
            FileNotFoundError: If config file doesn't exist
            ValueError: If configuration is invalid
            json.JSONDecodeError: If JSON is malformed
        """
        logger = get_logger()
        logger.info(f"Loading configuration from: {config_path}")
        
        config_file = Path(config_path)
        
        if not config_file.exists():
            error_msg = f"Configuration file not found: {config_path}"
            logger.error(error_msg)
            raise ValidationFileNotFoundError(error_msg)

        # Load JSON
        try:
            logger.debug(f"Reading configuration file...")
            with open(config_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            logger.debug(f"Configuration file parsed successfully")
        except json.JSONDecodeError as e:
            error_msg = f"Invalid JSON in configuration file: {e}"
            logger.error(error_msg)
            logger.add_validation_issue(
                level='ERROR',
                category='configuration',
                message='Malformed JSON in config file',
                details={'file': str(config_path), 'error': str(e)}
            )
            raise ConfigurationError(error_msg) from e
        
        # Create Config object
        config = Config()
        config.config_dir = config_file.parent
        
        # Validate and parse
        try:
            logger.debug("Validating configuration structure...")
            ConfigManager._validate_required_fields(data)
            logger.debug("Setup base for outputs")
            ConfigManager._setup_output_directory(config)
            logger.debug("Parsing genome configurations...")
            ConfigManager._parse_genome_configs(data, config)
            logger.debug("Parsing reads configurations...")
            ConfigManager._parse_reads_configs(data, config)
            logger.debug("Parsing feature configurations...")
            ConfigManager._parse_feature_configs(data, config)
            logger.debug("Parsing options...")
            ConfigManager._parse_options(data, config)

            logger.info("✓ Configuration loaded and validated successfully")
            return config
            
        except (ValueError, KeyError) as e:
            error_msg = f"Configuration validation failed: {e}"
            logger.error(error_msg)
            raise ConfigurationError(error_msg) from e
        
    
    @staticmethod
    def _validate_required_fields(data: dict):
        """Validate that required top-level fields exist."""
        required = ['ref_genome_filename', 'mod_genome_filename', 'reads']
        
        for field in required:
            if field not in data:
                raise ValueError(f"Missing required field: {field}")
        
        # Reads must be non-empty
        if not data['reads']:
            raise ValueError("'reads' must be a non-empty list")
    
    @staticmethod
    def _parse_genome_configs(data: dict, config: Config):
        """Parse genome and plasmid configurations."""
        # Required genomes
        config.ref_genome = ConfigManager._parse_genome_config(
            data['ref_genome_filename'], 'ref_genome_filename', config.config_dir, config.output_dir
        )
        config.mod_genome = ConfigManager._parse_genome_config(
            data['mod_genome_filename'], 'mod_genome_filename', config.config_dir, config.output_dir 
        )
        
        # Optional plasmids
        if 'ref_plasmid_filename' in data and data['ref_plasmid_filename']:
            config.ref_plasmid = ConfigManager._parse_genome_config(
                data['ref_plasmid_filename'], 'ref_plasmid_filename', config.config_dir, config.output_dir
            )
        
        if 'mod_plasmid_filename' in data and data['mod_plasmid_filename']:
            config.mod_plasmid = ConfigManager._parse_genome_config(
                data['mod_plasmid_filename'], 'mod_plasmid_filename', config.config_dir, config.output_dir
            )
    
    @staticmethod
    def _parse_genome_config(value: Any, field_name: str, config_dir: Path, output_dir: Path) -> GenomeConfig:
        """
        Parse a genome configuration entry.

        Args:
            value: Config value (dict with 'filename' or string)
            field_name: Name of field for error messages
            config_dir: Base directory for resolving relative paths
            output_dir: Output directory for validated files

        Returns:
            GenomeConfig with resolved paths and detected format/compression
        """
        from validation_pkg.validators.genome_validator import GenomeValidator
        logger = get_logger()

        # Parse filename and extra fields using unified utility
        filename, extra = ConfigManager._parse_config_file_value(value, field_name)

        # Extract validator-specific settings from extra fields
        settings_dict, remaining_extra = ConfigManager._extract_validator_settings(extra, GenomeValidator)

        # Log extracted settings
        if settings_dict:
            logger.debug(f"{field_name}: Extracted validator settings: {list(settings_dict.keys())}")

        # Warn about remaining extra fields that are not validator settings
        if remaining_extra:
            logger.warning(
                f"{field_name}: Found extra fields that will be ignored: {list(remaining_extra.keys())}"
            )

        # Resolve absolute path
        filepath = config_dir / filename

        if not filepath.exists():
            raise ValidationFileNotFoundError(
                f"The following file was not found: {filepath}\n")

        # Detect compression and format using unified utilities
        coding_type = ConfigManager._detect_compression_type(filepath)
        detected_format = ConfigManager._detect_file_format(filepath, GenomeFormat)

        return GenomeConfig(
            filename=filepath.name,
            filepath=filepath,
            coding_type=coding_type,
            detected_format=detected_format,
            extra=remaining_extra,
            settings_dict=settings_dict,
            output_dir=output_dir
        )
    
    @staticmethod
    def _parse_reads_configs(data: dict, config: Config):
        """Parse reads configuration."""
        logger = get_logger()
        reads_data = data['reads']

        if not isinstance(reads_data, list):
            raise ValueError("'reads' must be a list")

        for idx, read_entry in enumerate(reads_data):
            if not isinstance(read_entry, dict):
                raise ValueError(f"reads[{idx}] must be a dict")
            
            try:
                # Get filename or directory
                directory = read_entry.get('directory')
                filename = read_entry.get('filename')

                if not (filename or directory):
                    raise ValueError("'filename' or 'directory' field is required")

                if filename and directory:
                    raise ValueError("Cannot specify both 'filename' and 'directory' - use one or the other")

                if filename:
                    config.reads.append(ConfigManager._parse_read_config(read_entry,'filename',config.config_dir, config.output_dir))

                if directory:
                    #   Resolve absolute path directory
                    dirpath = config.config_dir / directory
                    if not dirpath.exists():
                        raise ValidationFileNotFoundError(
                            f"Directory not found: {dirpath}"
                        )

                    if not dirpath.is_dir():
                        raise ValueError(f"Path is not a directory: {dirpath}")

                    files = [Path(f) for f in dirpath.iterdir()]
                    if not files:
                        logger.warning(
                            f"No valid read files found in directory: {dirpath}. "
                        )
                        continue

                    logger.info(f"Found {len(files)} read file(s) in directory: {dirpath}")

                    for file in files:
                        # Create a dict with filename and ngs_type from directory entry
                        file_entry = {
                            'filename': str(Path(directory) / file.name),
                            'ngs_type': read_entry.get('ngs_type')
                        }
                        # Include any other extra fields from the directory entry
                        for key, value in read_entry.items():
                            if key not in ['directory', 'filename', 'ngs_type']:
                                file_entry[key] = value

                        config.reads.append(ConfigManager._parse_read_config(file_entry, "", config.config_dir, config.output_dir))
            
            except ValueError as e:
                raise ValueError(f"Invalid reads[{idx}]: {e}")
    
    @staticmethod
    def _parse_read_config(value: Any, field_name: str, config_dir: Path, output_dir: Path) -> ReadConfig:
        """
        Parse a read configuration entry.

        Args:
            value: Config value (dict with 'filename' and 'ngs_type' or string)
            field_name: Name of field for error messages
            config_dir: Base directory for resolving relative paths
            output_dir: Output directory for validated files

        Returns:
            ReadConfig with resolved paths and detected format/compression
        """
        from validation_pkg.validators.read_validator import ReadValidator
        logger = get_logger()

        # Parse filename and extra fields using unified utility
        filename, extra = ConfigManager._parse_config_file_value(value, field_name)

        # Extract ngs_type first (it's a required field, not a validator setting)
        ngs_type = extra.pop('ngs_type', 'illumina')

        if not ngs_type:
            raise ValueError("Missing required 'ngs_type'")

        # Extract validator-specific settings from remaining extra fields
        settings_dict, remaining_extra = ConfigManager._extract_validator_settings(extra, ReadValidator)

        # Log extracted settings
        if settings_dict:
            logger.debug(f"{field_name}: Extracted validator settings: {list(settings_dict.keys())}")

        # Warn about remaining extra fields that are not validator settings
        if remaining_extra:
            logger.warning(
                f"{field_name}: Found extra fields that will be ignored: {list(remaining_extra.keys())}"
            )

        # Resolve absolute path
        filepath = config_dir / filename

        if not filepath.exists():
            raise ValidationFileNotFoundError(
                f"The following file was not found: {filepath}\n")

        # Detect compression and format using unified utilities
        coding_type = ConfigManager._detect_compression_type(filepath)
        detected_format = ConfigManager._detect_file_format(filepath, ReadFormat)

        return ReadConfig(
            filename=filepath.name,
            filepath=filepath,
            ngs_type=ngs_type,
            coding_type=coding_type,
            detected_format=detected_format,
            output_dir=output_dir,
            extra=remaining_extra,
            settings_dict=settings_dict
        )

    @staticmethod
    def _parse_feature_configs(data: dict, config: Config):
        """Parse feature configurations."""
        if 'ref_feature_filename' in data and data['ref_feature_filename']:
            config.ref_feature = ConfigManager._parse_feature_config(
                data['ref_feature_filename'], 'ref_feature_filename', config.config_dir, config.output_dir
            )
        
        if 'mod_feature_filename' in data and data['mod_feature_filename']:
            config.mod_feature = ConfigManager._parse_feature_config(
                data['mod_feature_filename'], 'mod_feature_filename', config.config_dir, config.output_dir
            )
    
    @staticmethod
    def _parse_feature_config(value: Any, field_name: str, config_dir: Path, output_dir : Path) -> FeatureConfig:
        """
        Parse a single feature config entry and resolve paths.

        Uses unified utilities from file_handler for consistent parsing.

        Args:
            value: Config value (dict with 'filename' or string)
            field_name: Name of field for error messages
            config_dir: Base directory for resolving relative paths

        Returns:
            FeatureConfig with resolved paths and detected format/compression

        Raises:
            ValueError: If value format is invalid
        """
        from validation_pkg.validators.feature_validator import FeatureValidator
        logger = get_logger()

        # Parse filename and extra fields using unified utility
        filename, extra = ConfigManager._parse_config_file_value(value, field_name)

        # Extract validator-specific settings from extra fields
        settings_dict, remaining_extra = ConfigManager._extract_validator_settings(extra, FeatureValidator)

        # Log extracted settings
        if settings_dict:
            logger.debug(f"{field_name}: Extracted validator settings: {list(settings_dict.keys())}")

        # Warn about remaining extra fields that are not validator settings
        if remaining_extra:
            logger.warning(
                f"{field_name}: Found extra fields that will be ignored: {list(remaining_extra.keys())}"
            )

        # Resolve absolute path and file existence
        filepath = config_dir / filename
        if not filepath.exists():
            raise ValidationFileNotFoundError(
                f"The following file was not found: {filepath}\n")

        # Detect compression and format using unified utilities
        coding_type = ConfigManager._detect_compression_type(filepath)
        detected_format = ConfigManager._detect_file_format(filepath, FeatureFormat)

        return FeatureConfig(
            filename=filepath.name,
            filepath=filepath,
            coding_type=coding_type,
            detected_format=detected_format,
            output_dir=output_dir,
            extra=remaining_extra,
            settings_dict=settings_dict
        )
        
    @staticmethod
    def _parse_options(data: dict, config: Config):
        """
        Parse optional configuration options.

        Supported options:
        - threads: Number of threads for parallel compression (positive integer)

        Args:
            data: Configuration dictionary
            config: Config object to update

        Raises:
            ValueError: If threads is invalid (not positive integer)
        """
        logger = get_logger()

        if 'options' not in data:
            logger.debug("No options specified in config - using defaults")
            return

        options = data['options']

        if not isinstance(options, dict):
            raise ValueError("'options' must be a dictionary")

        # Parse threads option
        if 'threads' in options:
            threads = options['threads']

            # Allow null/None to mean auto-detect
            if threads is None:
                logger.info("Thread count: auto-detect (null specified)")
                config.options['threads'] = None
            elif isinstance(threads, int):
                # Validate positive integer
                if threads <= 0:
                    raise ValueError(f"'threads' must be a positive integer, got {threads}")

                # Warn if excessive (diminishing returns beyond 16)
                if threads > 16:
                    logger.warning(
                        f"Thread count {threads} is high - diminishing returns beyond 16 threads. "
                        f"Consider using 4-8 threads for optimal performance."
                    )

                logger.info(f"Using {threads} threads for parallel compression")
                config.options['threads'] = threads
            else:
                raise ValueError(
                    f"'threads' must be an integer or null, got {type(threads).__name__}: {threads}"
                )
        else:
            logger.info("Thread count: auto-detect (not specified)")

        # Store any other options for future use
        for key, value in options.items():
            if key not in ['threads']:
                config.options[key] = value
                logger.debug(f"Stored option '{key}': {value}")

    @staticmethod
    def _extract_validator_settings(
        extra: Dict[str, Any],
        validator_class: Type
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """
        Extract validator-specific settings from config extra fields.

        This method separates settings that belong to a validator's Settings class
        from other config fields. It uses the validator's Settings class field names
        to determine which keys are valid validator settings.

        Args:
            extra: Dictionary of extra fields from config.json
            validator_class: Validator class (GenomeValidator, ReadValidator, or FeatureValidator)

        Returns:
            Tuple of (settings_dict, remaining_extra):
                - settings_dict: Dict with valid validator settings
                - remaining_extra: Dict with non-validator fields

        Example:
            >>> extra = {'validation_level': 'trust', 'ngs_type': 'illumina', 'custom_field': 'value'}
            >>> settings, remaining = _extract_validator_settings(extra, ReadValidator)
            >>> settings
            {'validation_level': 'trust'}
            >>> remaining
            {'ngs_type': 'illumina', 'custom_field': 'value'}
        """
        from dataclasses import fields as dataclass_fields

        # Get valid field names from the validator's Settings class
        valid_settings = {f.name for f in dataclass_fields(validator_class.Settings)}

        # Separate validator settings from other fields
        settings_dict = {}
        remaining_extra = {}

        for key, value in extra.items():
            if key in valid_settings:
                settings_dict[key] = value
            else:
                remaining_extra[key] = value

        return settings_dict, remaining_extra

    @staticmethod
    def _setup_output_directory(config: Config):
        """
        Set up output directory for validation results.

        Creates the output directory at config_dir/output/ if it doesn't exist.
        Sets config.output_dir to the created directory path.

        Note: Individual validators can use output_subdir_name in their settings
        to create subdirectories within this base output directory.

        Args:
            config: Config object to update with output directory

        Raises:
            OSError: If directory creation fails
        """
        logger = get_logger()

        # Create output directory path
        output_dir = config.config_dir.parent / "valid"

        # Create directory if it doesn't exist
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            config.output_dir = output_dir
            logger.info(f"Output directory: {output_dir}")
            logger.debug(f"Output directory created/verified: {output_dir}")
        except OSError as e:
            error_msg = f"Failed to create output directory: {e}"
            logger.error(error_msg)
            raise ConfigurationError(error_msg) from e
        
    @staticmethod
    def _detect_compression_type(filepath: Path) -> CodingType:
        """
        Detect compression type from file path and return CodingType enum.

        Checks the last extension to determine compression:
        - .gz or .gzip → CodingType.GZIP
        - .bz2 or .bzip2 → CodingType.BZIP2
        - other → CodingType.NONE

        Args:
            filepath: Path to file

        Returns:
            CodingType enum indicating compression

        Examples:
            >>> detect_compression_type(Path('genome.fasta'))
            CodingType.NONE
            >>> detect_compression_type(Path('genome.fasta.gz'))
            CodingType.GZIP
            >>> detect_compression_type(Path('genome.fasta.gzip'))
            CodingType.GZIP
            >>> detect_compression_type(Path('reads.fastq.bz2'))
            CodingType.BZIP2

        Note:
            The function will detect .gz from file.tar.gz, which will
            cause issues since validators do not handle TAR extraction.
            Users must extract TAR archives before processing.
        """
        suffixes = Path(filepath).suffixes

        if not suffixes:
            return CodingType.NONE

        # Check last extension for compression
        last_ext = suffixes[-1].lower()

        if last_ext in ['.gz', '.gzip']:
            return CodingType.GZIP
        elif last_ext in ['.bz2', '.bzip2']:
            return CodingType.BZIP2
        else:
            return CodingType.NONE

    @staticmethod
    def _detect_file_format(filepath: Path, format_enum: Type[Union[GenomeFormat, ReadFormat, FeatureFormat]]) -> Union[GenomeFormat, ReadFormat, FeatureFormat]:
        """
        Detect file format from file path and return appropriate enum.

        Checks the first extension (before compression) to determine format.

        Args:
            filepath: Path to file
            format_enum: Format enum class (GenomeFormat, ReadFormat, or FeatureFormat)

        Returns:
            Format enum member

        Raises:
            ValueError: If format cannot be determined

        Examples:
            >>> detect_file_format(Path('genome.fasta'), GenomeFormat)
            GenomeFormat.FASTA
            >>> detect_file_format(Path('genome.fasta.gz'), GenomeFormat)
            GenomeFormat.FASTA
            >>> detect_file_format(Path('reads.fastq.bz2'), ReadFormat)
            ReadFormat.FASTQ
            >>> detect_file_format(Path('features.gff'), FeatureFormat)
            FeatureFormat.GFF
        """
        suffixes = Path(filepath).suffixes

        if not suffixes:
            raise ValueError(f"Cannot determine format: no extension found in {filepath.name}")

        # Get first extension (the format, before compression)
        format_ext = suffixes[0]

        # Let the format enum's _missing_ method handle the conversion
        return format_enum(format_ext)

    @staticmethod
    def _parse_config_file_value(
        value: Any,
        field_name: str
    ) -> Tuple[str, Dict[str, Any]]:
        """
        Parse file configuration value from JSON config.

        Handles both dict and string formats:
        - Dict: {"filename": "genome.fasta", "extra_key": "value"}
        - String: "genome.fasta" (backwards compatibility)

        Args:
            value: Configuration value (dict or string)
            field_name: Name of configuration field (for error messages)

        Returns:
            Tuple of (filename, extra_fields_dict)

        Raises:
            ValueError: If value is invalid format

        Examples:
            >>> _parse_config_file_value({"filename": "genome.fasta"}, "ref_genome")
            ("genome.fasta", {})
            >>> _parse_config_file_value({"filename": "reads.fastq", "ngs_type": "illumina"}, "reads")
            ("reads.fastq", {"ngs_type": "illumina"})
            >>> _parse_config_file_value("genome.fasta", "ref_genome")
            ("genome.fasta", {})
        """
        filename = None
        extra = {}

        if isinstance(value, dict):
            if 'filename' not in value:
                raise ValueError(f"{field_name} must contain 'filename' field")
            filename = value['filename']

            # Extract extra keys (excluding 'filename')
            extra = {k: v for k, v in value.items() if k != 'filename'}
                    # Warn about extra fields

        elif isinstance(value, str):
            # Accept plain string for backwards compatibility
            filename = value
        else:
            raise ValueError(f"{field_name} must be a dict or string, got {type(value).__name__}")

        return filename, extra