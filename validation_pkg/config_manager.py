"""
Configuration and ConfigManager classes for bioinformatics file validation.

Contains:
- GenomeConfig, ReadConfig, FeatureConfig: Dataclasses for file configurations
- Config: Main configuration container
- ConfigManager: Configuration file parser and validator

Note: ConfigManager class may be renamed to ConfigManager in future.
"""

import json
from pathlib import Path
from typing import List, Optional, Dict, Any, Union, Type, Tuple
from dataclasses import dataclass

# Import logging and exceptions
from validation_pkg.utils.formats import CodingType, GenomeFormat, ReadFormat, FeatureFormat
from validation_pkg.logger import get_logger
from validation_pkg.exceptions import (
    ConfigurationError,
    FileNotFoundError as ValidationFileNotFoundError
)

#   TODO: rescrict user of the package to call anything else except ConfigManager.load()
#   TODO[OPTIONAL]: rename filename to file?

@dataclass
class GenomeConfig:
    """Configuration for genome or plasmid files."""
    filename: str  # Original filename
    filepath: Path  # Absolute resolved path
    coding_type: CodingType = None  # Coding type applied on file
    detected_format: GenomeFormat = None  # Genome file format
    output_dir: Path = None  # Output folder base - heritage from config
    
    # Store any additional keys from config
    extra: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.extra is None:
            self.extra = {}

@dataclass
class ReadConfig:
    """Configuration for read files."""
    filename: str  # Original filename
    filepath: Path  # Absolute resolved path
    ngs_type: str = "illumina"  # "illumina", "ont", "pacbio"
    coding_type: CodingType = None # Coding type applied on file
    detected_format: ReadFormat = None  # Read file format
    output_dir: Path = None  # Output folder base - heritage from config

    # Store any additional keys from config
    extra: Dict[str, Any] = None

    def __post_init__(self):
        if self.ngs_type not in ["illumina", "ont", "pacbio"]:
            raise ValueError(f"Invalid ngs_type: {self.ngs_type}")
        if self.extra is None:
            self.extra = {}

@dataclass
class FeatureConfig:
    """Configuration for feature files."""
    filename: str  # Original relative filename
    filepath: Path  # Absolute resolved path
    coding_type: CodingType = None # Coding type applied on file
    detected_format: FeatureFormat = None  # Genome file format
    output_dir: Path = None  # Output folder base - heritage from config
    
    # Store any additional keys from config
    extra: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.extra is None:
            self.extra = {}


class Config:
    """
    Data holder for configuration.
    Stores validated configuration data without business logic.
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
        self.options: Dict[str, Any] = {}   #   TODO: implement release/debug mode and threads to run
        
        # Store config directory for relative path resolution
        self.config_dir: Optional[Path] = None
        
        # Output directory (set to config_dir.parent / "output")
        self.output_dir: Optional[Path] = None
    
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
    Reads and validates configuration files.
    Responsible for parsing JSON and validating keys/values.
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
            logger.debug("Parsing read configurations...")
            ConfigManager._parse_read_configs(data, config)
            logger.debug("Parsing feature configurations...")
            ConfigManager._parse_feature_configs(data, config)
            logger.debug("Parsing options...")
            ConfigManager._parse_options(data, config)
            logger.debug("Validating file existence...")
            ConfigManager._validate_file_existence(config)

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
            data['ref_genome_filename'], 'ref_genome_filename', config.config_dir
        )
        config.mod_genome = ConfigManager._parse_genome_config(
            data['mod_genome_filename'], 'mod_genome_filename', config.config_dir
        )
        
        # Optional plasmids
        if 'ref_plasmid_filename' in data and data['ref_plasmid_filename']:
            config.ref_plasmid = ConfigManager._parse_genome_config(
                data['ref_plasmid_filename'], 'ref_plasmid_filename', config.config_dir
            )
        
        if 'mod_plasmid_filename' in data and data['mod_plasmid_filename']:
            config.mod_plasmid = ConfigManager._parse_genome_config(
                data['mod_plasmid_filename'], 'mod_plasmid_filename', config.config_dir
            )
    
    @staticmethod
    def _parse_genome_config(value: Any, field_name: str, config_dir: Path) -> GenomeConfig:
        """
        Parse a genome condiguration.
        """
        logger = get_logger()

        # Parse filename and extra fields using unified utility
        filename, extra = ConfigManager._parse_config_file_value(value, field_name)

        # Warn about extra fields - we are not able to proccess them yet
        if extra:
            logger.warning(
                f"{field_name}: Found extra fields that will be ignored: {list(extra.keys())}"
            )

        # Resolve absolute path
        filepath = config_dir / filename
        output_dir_base = config_dir.parent / "valid"

        # Detect compression and format using unified utilities
        coding_type = ConfigManager._detect_compression_type(filepath)
        detected_format = ConfigManager._detect_file_format(filepath, GenomeFormat)

        return GenomeConfig(
            filename=filepath.name,
            filepath=filepath,
            coding_type=coding_type,
            detected_format=detected_format,
            extra=extra,
            output_dir=output_dir_base
        )
    
    @staticmethod
    def _parse_read_configs(data: dict, config: Config):
        """Parse reads configuration."""
        logger = get_logger()
        reads_data = data['reads']

        if not isinstance(reads_data, list):
            raise ValueError("'reads' must be a list")

        output_base_dir = config.config_dir.parent / "valid"
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

                ngs_type = read_entry.get('ngs_type', 'illumina')

                # Extract extra keys
                extra = {k: v for k, v in read_entry.items()
                        if k not in ['filename', 'directory', 'ngs_type']}

                # Warn about extra fields
                if extra:
                    logger.warning(
                        f"reads[{idx}]: Found extra fields that will be ignored: {list(extra.keys())}"
                    )

                if filename:
                    # Resolve absolute path
                    filepath = config.config_dir / filename

                    # Detect compression and format using unified utilities
                    coding_type = ConfigManager._detect_compression_type(filepath)
                    detected_format = ConfigManager._detect_file_format(filepath, ReadFormat)

                    read_config = ReadConfig(
                        filename=filepath.name,
                        filepath=filepath,
                        ngs_type=ngs_type,
                        coding_type=coding_type,
                        detected_format=detected_format,
                        output_dir=output_base_dir,
                        extra=extra
                    )
                    config.reads.append(read_config)
                
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
                        try:
                            filepath = Path(file)

                            # Detect compression and format using unified utilities
                            coding_type = ConfigManager._detect_compression_type(filepath)
                            detected_format = ConfigManager._detect_file_format(filepath, ReadFormat)
                            
                            #   NOTE: BAM file will be ignore in this version of package
                            if detected_format == ReadFormat.BAM:
                                logger.warning(
                                    f"BAM file will be ignored: {file}. "
                                )
                                continue

                            read_config = ReadConfig(
                                filename=filepath.name,
                                filepath=filepath,
                                ngs_type=ngs_type,
                                coding_type=coding_type,
                                detected_format=detected_format,
                                extra=extra
                            )
                            config.reads.append(read_config)
                        except ValueError as file_error:
                            logger.warning(f"Skipping file {file.name}: {file_error}")
            
            except ValueError as e:
                raise ValueError(f"Invalid reads[{idx}]: {e}")
    
    @staticmethod
    def _parse_feature_configs(data: dict, config: Config):
        """Parse feature configurations."""
        if 'ref_feature_filename' in data and data['ref_feature_filename']:
            config.ref_feature = ConfigManager._parse_feature_config(
                data['ref_feature_filename'], 'ref_feature_filename', config.config_dir
            )
        
        if 'mod_feature_filename' in data and data['mod_feature_filename']:
            config.mod_feature = ConfigManager._parse_feature_config(
                data['mod_feature_filename'], 'mod_feature_filename', config.config_dir
            )
    
    @staticmethod
    def _parse_feature_config(value: Any, field_name: str, config_dir: Path) -> FeatureConfig:
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
        logger = get_logger()

        # Parse filename and extra fields using unified utility
        filename, extra = ConfigManager._parse_config_file_value(value, field_name)

        # Warn about extra fields
        if extra:
            logger.warning(
                f"{field_name}: Found extra fields that will be ignored: {list(extra.keys())}"
            )

        # Resolve absolute path
        filepath = config_dir / filename

        # Detect compression and format using unified utilities
        coding_type = ConfigManager._detect_compression_type(filepath)
        detected_format = ConfigManager._detect_file_format(filepath, FeatureFormat)

        return FeatureConfig(
            filename=filepath.name,
            filepath=filepath,
            coding_type=coding_type,
            detected_format=detected_format,
            extra=extra
        )
        
    @staticmethod
    def _parse_options(data: dict, config: Config):
        """
        Parse optional configuration options.

        Note: Options functionality not yet fully implemented.
        Reserved for future configuration extensions.
        """
        if 'options' in data:
            config.options = data['options']
    
    @staticmethod
    def _validate_file_existence(config: Config):
        """Validate that all specified files exist."""
        files_to_check = []

        # Genomes
        if config.ref_genome:
            files_to_check.append(('ref_genome', config.ref_genome.filepath))
        if config.mod_genome:
            files_to_check.append(('mod_genome', config.mod_genome.filepath))

        # Plasmids
        if config.ref_plasmid:
            files_to_check.append(('ref_plasmid', config.ref_plasmid.filepath))
        if config.mod_plasmid:
            files_to_check.append(('mod_plasmid', config.mod_plasmid.filepath))

        # Features
        if config.ref_feature:
            files_to_check.append(('ref_feature', config.ref_feature.filepath))
        if config.mod_feature:
            files_to_check.append(('mod_feature', config.mod_feature.filepath))

        # Reads
        for idx, read in enumerate(config.reads):
            if read.filepath:
                files_to_check.append((f'reads[{idx}]', read.filepath))

        # Check all files exist
        missing_files = []
        for name, filepath in files_to_check:
            if not filepath.exists():
                missing_files.append(f"{name}: {filepath}")

        if missing_files:
            raise ValidationFileNotFoundError(
                f"The following files were not found:\n" +
                "\n".join(f"  - {f}" for f in missing_files)
            )
    
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

        elif isinstance(value, str):
            # Accept plain string for backwards compatibility
            filename = value
        else:
            raise ValueError(f"{field_name} must be a dict or string, got {type(value).__name__}")

        return filename, extra

    @staticmethod
    def _validate_file_exists(filepath: Path, field_name: str) -> None:
        """
        Validate that file exists and is readable.

        Args:
            filepath: Path to file
            field_name: Name of configuration field (for error messages)

        Raises:
            FileNotFoundError: If file does not exist
            ValueError: If path is not a file
            PermissionError: If file is not readable

        Examples:
            >>> validate_file_exists(Path('genome.fasta'), 'ref_genome')
            # Raises FileNotFoundError if file doesn't exist
        """
        if not filepath.exists():
            raise FileNotFoundError(
                f"{field_name}: File not found: {filepath}\n"
                f"  Make sure the file exists relative to the config file location."
            )

        if not filepath.is_file():
            raise ValueError(f"{field_name}: Path is not a file: {filepath}")

        # Check if file is readable
        try:
            with open(filepath, 'rb') as f:
                f.read(1)
        except PermissionError as e:
            raise PermissionError(f"{field_name}: File is not readable: {filepath}") from e