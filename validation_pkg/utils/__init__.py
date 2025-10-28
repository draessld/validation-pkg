"""
Utility modules for validation_pkg.

Provides:
    - file_handler: File I/O, compression handling, format detection
    - formats: File format and compression type enums
    - settings: Base settings class with immutable update pattern
"""

from .file_handler import (
    open_file,
    detect_compression,
    get_base_filename,
    get_file_format,
    calculate_thread_distribution
)

__all__ = [
    'open_file',
    'detect_compression',
    'get_base_filename',
    'get_file_format',
    'calculate_thread_distribution',
]