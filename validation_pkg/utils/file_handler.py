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
import shutil


# Cache for tool availability checks to avoid repeated subprocess calls
_TOOL_CACHE = {}
_LOGGER_INITIALIZED = False


def check_compression_tool_available(tool_name: str) -> bool:
    """
    Check if a compression tool is available on the system.

    Results are cached to avoid repeated subprocess calls.

    Args:
        tool_name: Name of the tool to check ('pigz', 'pbzip2', 'gzip', 'bzip2')

    Returns:
        True if tool is available, False otherwise

    Example:
        >>> if check_compression_tool_available('pigz'):
        ...     print("pigz is available for faster compression")
    """
    if tool_name in _TOOL_CACHE:
        return _TOOL_CACHE[tool_name]

    # Check if tool exists using shutil.which
    available = shutil.which(tool_name) is not None
    _TOOL_CACHE[tool_name] = available

    return available


def get_optimal_thread_count() -> int:
    """
    Get optimal thread count for parallel compression.

    Auto-detects CPU cores and caps at 8 threads for balanced performance.
    More than 8 threads shows diminishing returns for compression.

    Returns:
        Number of threads to use (1-8)

    Example:
        >>> threads = get_optimal_thread_count()
        >>> print(f"Using {threads} threads for compression")
    """
    try:
        cpu_count = os.cpu_count() or 1
        # Cap at 8 threads - diminishing returns beyond this for compression
        return min(cpu_count, 8)
    except:
        return 4  # Safe default


def get_compression_command(coding_type: CodingType, mode: str = 'compress', threads: int = None) -> tuple:
    """
    Get the best available compression command for the given coding type.

    Prefers parallel tools (pigz, pbzip2) over standard tools (gzip, bzip2).
    Falls back gracefully if parallel tools are not available.

    Args:
        coding_type: CodingType enum (GZIP, BZIP2, or NONE)
        mode: 'compress' or 'decompress'
        threads: Number of threads to use (None = auto-detect)

    Returns:
        Tuple of (command_name, [args]) for subprocess

    Example:
        >>> cmd, args = get_compression_command(CodingType.GZIP, 'compress', threads=4)
        >>> # Returns ('pigz', ['-c', '-p', '4']) if pigz available
        >>> # Otherwise ('gzip', ['-c'])
    """
    global _LOGGER_INITIALIZED

    if threads is None:
        threads = get_optimal_thread_count()

    if coding_type == CodingType.GZIP:
        # Try pigz first (parallel gzip)
        if check_compression_tool_available('pigz'):
            if not _LOGGER_INITIALIZED:
                try:
                    from validation_pkg.logger import get_logger
                    logger = get_logger()
                    logger.info(f"Using pigz for gzip compression ({threads} threads)")
                    _LOGGER_INITIALIZED = True
                except:
                    pass

            if mode == 'decompress':
                return ('pigz', ['-dc', '-p', str(threads)])
            else:
                return ('pigz', ['-c', '-p', str(threads)])
        else:
            # Fallback to standard gzip
            if not _LOGGER_INITIALIZED:
                try:
                    from validation_pkg.logger import get_logger
                    logger = get_logger()
                    logger.info("Using standard gzip (install pigz for better performance: sudo apt-get install pigz)")
                    _LOGGER_INITIALIZED = True
                except:
                    pass

            if mode == 'decompress':
                return ('gzip', ['-dc'])
            else:
                return ('gzip', ['-c'])

    elif coding_type == CodingType.BZIP2:
        # Try pbzip2 first (parallel bzip2)
        if check_compression_tool_available('pbzip2'):
            if not _LOGGER_INITIALIZED:
                try:
                    from validation_pkg.logger import get_logger
                    logger = get_logger()
                    logger.info(f"Using pbzip2 for bzip2 compression ({threads} threads)")
                    _LOGGER_INITIALIZED = True
                except:
                    pass

            if mode == 'decompress':
                return ('pbzip2', ['-dc', '-p' + str(threads)])
            else:
                return ('pbzip2', ['-c', '-p' + str(threads)])
        else:
            # Fallback to standard bzip2
            if not _LOGGER_INITIALIZED:
                try:
                    from validation_pkg.logger import get_logger
                    logger = get_logger()
                    logger.info("Using standard bzip2 (install pbzip2 for better performance: sudo apt-get install pbzip2)")
                    _LOGGER_INITIALIZED = True
                except:
                    pass

            if mode == 'decompress':
                return ('bzip2', ['-dc'])
            else:
                return ('bzip2', ['-c'])

    else:
        # No compression
        return ('cat', [])


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

    # Remove compression extension if present (case-insensitive)
    if filepath.suffix.lower() in ['.gz', '.bz2']:
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


def gz_to_bz2(gz_file: Path, bz2_file: Path):
    """
    Convert gzip compressed file to bzip2.

    Uses parallel tools (pigz/pbzip2) if available, falls back to gzip/bzip2.

    Args:
        gz_file: Path to input .gz file
        bz2_file: Path to output .bz2 file
    """
    # Get best available compression commands
    decompress_cmd, decompress_args = get_compression_command(CodingType.GZIP, 'decompress')
    compress_cmd, compress_args = get_compression_command(CodingType.BZIP2, 'compress')

    # Build piped command
    cmd = f"set -o pipefail && {decompress_cmd} {' '.join(decompress_args)} {gz_file} | {compress_cmd} {' '.join(compress_args)} > {bz2_file}"

    subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')


def bz2_to_gz(bz2_file: Path, gz_file: Path):
    """
    Convert bzip2 compressed file to gzip.

    Uses parallel tools (pbzip2/pigz) if available, falls back to bzip2/gzip.

    Args:
        bz2_file: Path to input .bz2 file
        gz_file: Path to output .gz file
    """
    # Get best available compression commands
    decompress_cmd, decompress_args = get_compression_command(CodingType.BZIP2, 'decompress')
    compress_cmd, compress_args = get_compression_command(CodingType.GZIP, 'compress')

    # Build piped command
    cmd = f"set -o pipefail && {decompress_cmd} {' '.join(decompress_args)} {bz2_file} | {compress_cmd} {' '.join(compress_args)} > {gz_file}"

    subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')

def none_to_gz(none_file: Path, gz_file: Path):
    """
    Compress uncompressed file to gzip.

    Uses pigz if available for parallel compression, falls back to gzip.

    Args:
        none_file: Path to input uncompressed file
        gz_file: Path to output .gz file
    """
    compress_cmd, compress_args = get_compression_command(CodingType.GZIP, 'compress')

    cmd = f"{compress_cmd} {' '.join(compress_args)} {none_file} > {gz_file}"
    subprocess.run(cmd, shell=True, check=True)

def gz_to_none(gz_file: Path, none_file: Path):
    """
    Decompress gzip file to uncompressed file.

    Uses pigz if available for parallel decompression, falls back to gzip.

    Args:
        gz_file: Path to input .gz file
        none_file: Path to output uncompressed file
    """
    decompress_cmd, decompress_args = get_compression_command(CodingType.GZIP, 'decompress')

    cmd = f"{decompress_cmd} {' '.join(decompress_args)} {gz_file} > {none_file}"
    subprocess.run(cmd, shell=True, check=True)

def bz2_to_none(bz2_file: Path, none_file: Path):
    """
    Decompress bzip2 file to uncompressed file.

    Uses pbzip2 if available for parallel decompression, falls back to bzip2.

    Args:
        bz2_file: Path to input .bz2 file
        none_file: Path to output uncompressed file
    """
    decompress_cmd, decompress_args = get_compression_command(CodingType.BZIP2, 'decompress')

    cmd = f"{decompress_cmd} {' '.join(decompress_args)} {bz2_file} > {none_file}"
    subprocess.run(cmd, shell=True, check=True)

def none_to_bz2(none_file: Path, bz2_file: Path):
    """
    Compress uncompressed file to bzip2.

    Uses pbzip2 if available for parallel compression, falls back to bzip2.

    Args:
        none_file: Path to input uncompressed file
        bz2_file: Path to output .bz2 file
    """
    compress_cmd, compress_args = get_compression_command(CodingType.BZIP2, 'compress')

    cmd = f"{compress_cmd} {' '.join(compress_args)} {none_file} > {bz2_file}"
    subprocess.run(cmd, shell=True, check=True)


def open_compressed_writer(filepath: Union[str, Path], coding_type: CodingType, use_parallel: bool = True):
    """
    Open a file handle for writing compressed data efficiently.

    Prefers parallel compression tools (pigz/pbzip2) when use_parallel=True.
    Falls back to Python gzip/bz2 libraries for compatibility.

    This function returns a context manager that should be used with 'with' statement.

    Args:
        filepath: Path to output file
        coding_type: CodingType enum (GZIP, BZIP2, or NONE)
        use_parallel: If True, use parallel tools (pigz/pbzip2) when available

    Returns:
        File handle for writing (text mode)

    Example:
        >>> from validation_pkg.utils.formats import CodingType
        >>> from pathlib import Path
        >>> with open_compressed_writer('output.gz', CodingType.GZIP) as f:
        ...     f.write("Hello, world!\\n")

    Note:
        For parallel tools to work, pigz or pbzip2 must be installed.
        If not available or use_parallel=False, falls back to Python libraries.
    """
    filepath = Path(filepath)

    # If parallel compression is requested and available, use subprocess approach
    if use_parallel and coding_type in (CodingType.GZIP, CodingType.BZIP2):
        # Check if parallel tools are available
        use_subprocess = False

        if coding_type == CodingType.GZIP:
            use_subprocess = check_compression_tool_available('pigz')
        elif coding_type == CodingType.BZIP2:
            use_subprocess = check_compression_tool_available('pbzip2')

        if use_subprocess:
            # Use subprocess pipe for parallel compression
            compress_cmd, compress_args = get_compression_command(coding_type, 'compress')

            # Create subprocess with pipe to stdin
            process = subprocess.Popen(
                [compress_cmd] + compress_args,
                stdin=subprocess.PIPE,
                stdout=open(filepath, 'wb'),
                stderr=subprocess.PIPE,
                text=False  # Binary mode for pipe
            )

            # Wrap in a context manager that handles text encoding and cleanup
            class SubprocessWriter:
                def __init__(self, proc):
                    self.proc = proc
                    self.stdin = proc.stdin

                def write(self, text: str):
                    """Write text to subprocess stdin."""
                    self.stdin.write(text.encode('utf-8'))

                def __enter__(self):
                    return self

                def __exit__(self, exc_type, exc_val, exc_tb):
                    self.stdin.close()
                    self.proc.wait()
                    if self.proc.returncode != 0:
                        stderr_output = self.proc.stderr.read().decode('utf-8')
                        raise CompressionError(
                            f"Compression failed with return code {self.proc.returncode}: {stderr_output}"
                        )
                    return False

            return SubprocessWriter(process)

    # Fallback to Python compression libraries
    if coding_type == CodingType.GZIP:
        return gzip.open(filepath, 'wt')
    elif coding_type == CodingType.BZIP2:
        return bz2.open(filepath, 'wt')
    else:
        # No compression
        return open(filepath, 'w')