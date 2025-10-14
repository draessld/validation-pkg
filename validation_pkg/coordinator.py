"""
Configuration and Coordinator classes for bioinformatics file validation.
"""

import json
from pathlib import Path
from typing import List, Optional, Dict, Any
from dataclasses import dataclass

# Import logging and exceptions
from validation_pkg.utils import settings
from validation_pkg.logger import get_logger
from validation_pkg.exceptions import (
    ConfigurationError,
    FileNotFoundError as ValidationFileNotFoundError
)



@dataclass
class GenomeConfig:
    """Configuration for genome or plasmid files."""
    filename: str  # Original relative filename
    filepath: Path  # Absolute resolved path
    
    # Store any additional keys from config
    extra: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.extra is None:
            self.extra = {}

@dataclass
class ReadConfig:
    """Configuration for read files."""
    filename: Optional[str] = None  # Original relative filename
    directory: Optional[str] = None  # Original relative directory
    filepath: Optional[Path] = None  # Absolute resolved path (for file)
    dirpath: Optional[Path] = None  # Absolute resolved path (for directory)
    ngs_type: str = "illumina"  # "illumina", "ont", "pacbio"
    
    # Store any additional keys from config
    extra: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.filename is None and self.directory is None:
            raise ValueError("Either filename or directory must be provided")
        if self.filename and self.directory:
            raise ValueError("Provide either filename or directory, not both")
        if self.ngs_type not in ["illumina", "ont", "pacbio"]:
            raise ValueError(f"Invalid ngs_type: {self.ngs_type}")
        if self.extra is None:
            self.extra = {}

@dataclass
class FeatureConfig:
    """Configuration for feature files."""
    filename: str  # Original relative filename
    filepath: Path  # Absolute resolved path
    
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
        self.options: Dict[str, Any] = {}
        
        # Store config directory for relative path resolution
        self.config_dir: Optional[Path] = None
        
        # Output directory (set to config_dir / "output")
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


class Coordinator:
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
            # Validate and parse
            logger.debug("Validating configuration structure...")
            Coordinator._validate_required_fields(data)
            logger.debug("Parsing genome configurations...")
            Coordinator._parse_genome_configs(data, config)
            logger.debug("Parsing read configurations...")
            Coordinator._parse_read_configs(data, config)
            logger.debug("Parsing feature configurations...")
            Coordinator._parse_feature_configs(data, config)
            logger.debug("Parsing options...")
            Coordinator._parse_options(data, config)
            logger.debug("Validating file existence...")
            Coordinator._validate_file_existence(config)
            logger.debug("Setup output directory...")
            
            logger.info("✓ Configuration loaded and validated successfully")
            return config
            
        except (ValueError, KeyError) as e:
            error_msg = f"Configuration validation failed: {e}"
            logger.error(error_msg)
            raise ConfigurationError(error_msg) from e
        
        
        return config
    
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
        config.ref_genome = Coordinator._parse_genome_config(
            data['ref_genome_filename'], 'ref_genome_filename', config.config_dir
        )
        config.mod_genome = Coordinator._parse_genome_config(
            data['mod_genome_filename'], 'mod_genome_filename', config.config_dir
        )
        
        # Optional plasmids
        if 'ref_plasmid_filename' in data and data['ref_plasmid_filename']:
            config.ref_plasmid = Coordinator._parse_genome_config(
                data['ref_plasmid_filename'], 'ref_plasmid_filename', config.config_dir
            )
        
        if 'mod_plasmid_filename' in data and data['mod_plasmid_filename']:
            config.mod_plasmid = Coordinator._parse_genome_config(
                data['mod_plasmid_filename'], 'mod_plasmid_filename', config.config_dir
            )
    
    @staticmethod
    def _parse_genome_config(value: Any, field_name: str, config_dir: Path) -> GenomeConfig:
        """Parse a single genome config entry and resolve paths."""
        filename = None
        extra = {}  
        if isinstance(value, dict):
            if 'filename' not in value:
                raise ValueError(f"{field_name} must contain 'filename' field")
            filename = value['filename']
            
            # Store any extra keys
            extra = {k: v for k, v in value.items() if k not in ['filename']}   
        elif isinstance(value, str):
            # Accept plain string for backwards compatibility
            filename = value
        else:
            raise ValueError(f"{field_name} must be a dict or string")
        
        # Resolve absolute path
        filepath = config_dir / filename
        
        return GenomeConfig(
            filename=filename,
            filepath=filepath,
            extra=extra)
    
    @staticmethod
    def _parse_read_configs(data: dict, config: Config):
        """Parse reads configuration."""
        reads_data = data['reads']
        
        if not isinstance(reads_data, list):
            raise ValueError("'reads' must be a list")
        
        for idx, read_entry in enumerate(reads_data):
            if not isinstance(read_entry, dict):
                raise ValueError(f"reads[{idx}] must be a dict")
            
            try:
                filename = read_entry.get('filename')
                directory = read_entry.get('directory')
                ngs_type = read_entry.get('ngs_type', 'illumina')
                
                # Extract extra keys
                extra = {k: v for k, v in read_entry.items() 
                        if k not in ['filename', 'directory', 'ngs_type']}
                
                # Resolve paths
                filepath = config.config_dir / filename if filename else None
                dirpath = config.config_dir / directory if directory else None
                
                read_config = ReadConfig(
                    filename=filename,
                    directory=directory,
                    filepath=filepath,
                    dirpath=dirpath,
                    ngs_type=ngs_type,
                    extra=extra
                )
                config.reads.append(read_config)
            except ValueError as e:
                raise ValueError(f"Invalid reads[{idx}]: {e}")
    
    @staticmethod
    def _parse_feature_configs(data: dict, config: Config):
        """Parse feature configurations."""
        if 'ref_feature_filename' in data and data['ref_feature_filename']:
            config.ref_feature = Coordinator._parse_feature_config(
                data['ref_feature_filename'], 'ref_feature_filename', config.config_dir
            )
        
        if 'mod_feature_filename' in data and data['mod_feature_filename']:
            config.mod_feature = Coordinator._parse_feature_config(
                data['mod_feature_filename'], 'mod_feature_filename', config.config_dir
            )
    
    @staticmethod
    def _parse_feature_config(value: Any, field_name: str, config_dir: Path) -> FeatureConfig:
        """Parse a single feature config entry and resolve paths."""
        filename = None
        extra = {}
        
        if isinstance(value, dict):
            if 'filename' not in value:
                raise ValueError(f"{field_name} must contain 'filename' field")
            filename = value['filename']
            
            # Store any extra keys
            extra = {k: v for k, v in value.items() 
                    if k not in ['filename']}
            
        elif isinstance(value, str):
            filename = value
        else:
            raise ValueError(f"{field_name} must be a dict or string")
        
        # Resolve absolute path
        filepath = config_dir / filename
        
        return FeatureConfig(
            filename=filename,
            filepath=filepath,
            extra=extra
        )
        
    @staticmethod
    def _parse_options(data: dict, config: Config):
        """Parse optional configuration options."""
        if 'options' in data:
            config.options = data['options']
    
    @staticmethod
    def _validate_file_existence(config: Config):
        """Validate that all specified files exist."""
        files_to_check = []
        
        # Genomes
        if config.ref_genome:
            files_to_check.append(('ref_genome', config.ref_genome.filename))
        if config.mod_genome:
            files_to_check.append(('mod_genome', config.mod_genome.filename))
        
        # Plasmids
        if config.ref_plasmid:
            files_to_check.append(('ref_plasmid', config.ref_plasmid.filename))
        if config.mod_plasmid:
            files_to_check.append(('mod_plasmid', config.mod_plasmid.filename))
        
        # Features
        if config.ref_feature:
            files_to_check.append(('ref_feature', config.ref_feature.filename))
        if config.mod_feature:
            files_to_check.append(('mod_feature', config.mod_feature.filename))
        
        # Reads
        for idx, read in enumerate(config.reads):
            if read.filename:
                files_to_check.append((f'reads[{idx}]', read.filename))
            elif read.directory:
                # Check directory exists
                dir_path = config.resolve_path(read.directory)
                if not dir_path.exists():
                    raise FileNotFoundError(f"reads[{idx}] directory not found: {read.directory}")
                if not dir_path.is_dir():
                    raise ValueError(f"reads[{idx}] path is not a directory: {read.directory}")
        
        # Check all files exist
        missing_files = []
        for name, filepath in files_to_check:
            full_path = config.resolve_path(filepath)
            if not full_path.exists():
                missing_files.append(f"{name}: {filepath}")
        
        if missing_files:
            raise FileNotFoundError(
                f"The following files were not found:\n" +
                "\n".join(f"  - {f}" for f in missing_files)
            )
    
    @staticmethod
    def _validate_file_format(filename: str, allowed_extensions: List[str]) -> bool:
        """
        Validate file format based on extension.
        
        Args:
            filename: File name to check
            allowed_extensions: List of allowed extensions (e.g., ['.fa', '.fasta'])
            
        Returns:
            True if valid, False otherwise
        """
        # Remove compression extensions first
        name = filename.lower()
        if name.endswith('.gz'):
            name = name[:-3]
        elif name.endswith('.bz2'):
            name = name[:-4]
        
        return any(name.endswith(ext) for ext in allowed_extensions)
    
        """
        Setup output path by default on tmp path
        """
        config.output_dir = config.config_dir / dirname