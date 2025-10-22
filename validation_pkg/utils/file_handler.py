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
from typing import Union, TextIO

from validation_pkg.utils.formats import CodingType, GenomeFormat, ReadFormat, FeatureFormat
from validation_pkg.exceptions import CompressionError

import subprocess
import os


def open_file_with_coding_type(
    filepath: Union[str, Path],
    coding_type: CodingType,
    mode: str = 'rt'
) -> TextIO:
    """
    Open a file with automatic decompression based on CodingType enum.

    This function provides centralized file opening logic for all validators,
    eliminating code duplication and ensuring consistent compression handling.

    Args:
        filepath: Path to file
        coding_type: CodingType enum indicating compression
        mode: Opening mode (default: 'rt' for text read)

    Returns:
        File handle with automatic decompression

    Raises:
        CompressionError: If file cannot be opened or decompressed

    Example:
        >>> from validation_pkg.utils.formats import CodingType
        >>> from pathlib import Path
        >>> filepath = Path('genome.fasta.gz')
        >>> coding = CodingType.GZIP
        >>> with open_file_with_coding_type(filepath, coding) as f:
        ...     content = f.read()
    """
    filepath = Path(filepath)

    try:
        if coding_type == CodingType.GZIP:
            return gzip.open(filepath, mode)
        elif coding_type == CodingType.BZIP2:
            return bz2.open(filepath, mode)
        else:
            # No compression or unknown
            return open(filepath, mode)
    except Exception as e:
        raise CompressionError(f"Failed to open file {filepath}: {e}") from e


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


def gz_to_bz2(gz_file:Path, bz2_file:Path):

    subprocess.run(
        f"gzip -dc {gz_file} | bzip2 -c > {bz2_file}",
        shell=True, check=True
    )


def bz2_to_gz(bz2_file:Path, gz_file:Path):
    subprocess.run(
        f"bzip2 -dc {bz2_file} | gzip -c > {gz_file}",
        shell=True, check=True
    )

def none_to_gz(none_file:Path, gz_file:Path):
    subprocess.run(
        f"gzip -c {none_file} > {gz_file}",
        shell=True, check=True
    )

def gz_to_none(gz_file:Path, none_file:Path):
    """Decompress gzip file to uncompressed file."""
    subprocess.run(
        f"gzip -dc {gz_file} > {none_file}",
        shell=True, check=True
    )

def bz2_to_none(bz2_file:Path, none_file:Path):
    """Decompress bzip2 file to uncompressed file."""
    subprocess.run(
        f"bzip2 -dc {bz2_file} > {none_file}",
        shell=True, check=True
    )

def none_to_bz2(none_file:Path, bz2_file:Path):
    """Compress uncompressed file to bzip2."""
    subprocess.run(
        f"bzip2 -c {none_file} > {bz2_file}",
        shell=True, check=True
    )