"""
Utility functions for file handling with compression support.

Provides unified utilities for:
- File opening with automatic decompression
- Compression detection
- Format detection
- File validation
- ConfigManager parsing helpers
"""

import gzip
import bz2
from pathlib import Path
from typing import Union, TextIO, Tuple, Dict, Any, Type
from enum import Enum

from validation_pkg.utils.formats import CodingType, GenomeFormat, ReadFormat, FeatureFormat


def open_file(filepath: Union[str, Path], mode: str = 'rt') -> TextIO:
    """
    Open a file with automatic decompression.
    
    Supports:
    - Plain files
    - Gzip compressed files (.gz)
    - Bzip2 compressed files (.bz2)
    
    Args:
        filepath: Path to file
        mode: Opening mode (default: 'rt' for text read)
        
    Returns:
        File handle
        
    Example:
        with open_file('genome.fasta.gz') as f:
            content = f.read()
    """
    filepath = Path(filepath)
    
    if filepath.suffix == '.gz':
        return gzip.open(filepath, mode)
    elif filepath.suffix == '.bz2':
        return bz2.open(filepath, mode)
    else:
        return open(filepath, mode)


def detect_compression(filepath: Union[str, Path]) -> str:
    """
    Detect compression type from file extension.
    
    Args:
        filepath: Path to file
        
    Returns:
        'gz', 'bz2', or 'none'
    """
    filepath = Path(filepath)
    suffix = filepath.suffix.lower()
    
    if suffix == '.gz':
        return 'gz'
    elif suffix == '.bz2':
        return 'bz2'
    else:
        return 'none'


def get_base_filename(filepath: Union[str, Path]) -> str:
    """
    Get base filename without compression extension.
    
    Args:
        filepath: Path to file
        
    Returns:
        Base filename
        
    Example:
        get_base_filename('genome.fasta.gz') -> 'genome.fasta'
        get_base_filename('reads.fastq.bz2') -> 'reads.fastq'
    """
    filepath = Path(filepath)
    
    # Remove compression extension if present
    if filepath.suffix in ['.gz', '.bz2']:
        return filepath.stem
    
    return filepath.name


def get_file_format(filepath: Union[str, Path]) -> str:
    """
    Detect file format from extension (ignoring compression).
    
    Args:
        filepath: Path to file
        
    Returns:
        File format: 'fasta', 'genbank', 'bed', 'gff', 'gtf', 'fastq', or 'unknown'
    """
    # Get filename without compression
    base_name = get_base_filename(filepath)
    base_path = Path(base_name)
    ext = base_path.suffix.lower()
    
    # Map extensions to formats
    format_map = {
        '.fa': 'fasta',
        '.fasta': 'fasta',
        '.fna': 'fasta',
        '.gb': 'genbank',
        '.gbk': 'genbank',
        '.genbank': 'genbank',
        '.bed': 'bed',
        '.gff': 'gff',
        '.gff3': 'gff',
        '.gtf': 'gtf',
        '.fastq': 'fastq',
        '.fq': 'fastq',
    }
    
    return format_map.get(ext, 'unknown')


# ===== Unified ConfigManager Utilities =====

def detect_compression_type(filepath: Path) -> CodingType:
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
        TAR archives (.tar.gz, .tgz) are NOT supported.
        The function will detect .gz from file.tar.gz, which will
        cause issues since validators do not handle TAR extraction.
        Users must extract TAR archives before processing.
    """
    suffixes = filepath.suffixes

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


def detect_file_format(
    filepath: Path,
    format_enum: Type[Union[GenomeFormat, ReadFormat, FeatureFormat]]
) -> Union[GenomeFormat, ReadFormat, FeatureFormat]:
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
    suffixes = filepath.suffixes

    if not suffixes:
        raise ValueError(f"Cannot determine format: no extension found in {filepath.name}")

    # Get first extension (the format, before compression)
    format_ext = suffixes[0]

    # Let the format enum's _missing_ method handle the conversion
    return format_enum(format_ext)


def parse_config_file_value(
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
        >>> parse_config_file_value({"filename": "genome.fasta"}, "ref_genome")
        ("genome.fasta", {})
        >>> parse_config_file_value({"filename": "reads.fastq", "ngs_type": "illumina"}, "reads")
        ("reads.fastq", {"ngs_type": "illumina"})
        >>> parse_config_file_value("genome.fasta", "ref_genome")
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


def validate_file_exists(filepath: Path, field_name: str) -> None:
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