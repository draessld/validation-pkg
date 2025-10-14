"""
Utility functions for path resolution and validation.
"""

from pathlib import Path
from typing import Union


def resolve_path(relative_path: Union[str, Path], base_dir: Path) -> Path:
    """
    Resolve a relative path to an absolute path based on a base directory.
    
    Args:
        relative_path: Relative path to resolve
        base_dir: Base directory to resolve from
        
    Returns:
        Absolute Path object
        
    Example:
        >>> resolve_path("data/genome.fasta", Path("/home/user/project"))
        Path("/home/user/project/data/genome.fasta")
    """
    relative_path = Path(relative_path)
    
    # If already absolute, return as is
    if relative_path.is_absolute():
        return relative_path
    
    # Resolve relative to base directory
    return (base_dir / relative_path).resolve()


def validate_path_exists(path: Path, path_type: str = "file") -> bool:
    """
    Validate that a path exists and is of the correct type.
    
    Args:
        path: Path to validate
        path_type: Type of path ("file" or "directory")
        
    Returns:
        True if valid, False otherwise
    """
    if not path.exists():
        return False
    
    if path_type == "file" and not path.is_file():
        return False
    
    if path_type == "directory" and not path.is_dir():
        return False
    
    return True


def get_relative_path(absolute_path: Path, base_dir: Path) -> str:
    """
    Get relative path from absolute path based on base directory.
    
    Args:
        absolute_path: Absolute path
        base_dir: Base directory
        
    Returns:
        Relative path as string
        
    Example:
        >>> get_relative_path(Path("/home/user/project/data/file.txt"), Path("/home/user/project"))
        "data/file.txt"
    """
    try:
        return str(absolute_path.relative_to(base_dir))
    except ValueError:
        # If paths are not relative, return absolute
        return str(absolute_path)