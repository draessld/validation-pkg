from enum import Enum
from pathlib import Path

# ===== Enums  =====

class CodingType(Enum):
    GZIP = "gzip"
    BZIP2 = "bzip2"
    TGZ = "tgz"
    NONE = "none"
    

    @classmethod
    def _missing_(cls, value):
        """
        Called when enum lookup fails. Allows CodingType('.gz') or CodingType('gz')
        to work alongside CodingType.GZIP.
        """
        # Normalize the input
        value_lower = str(value).lower().strip()
        
        # Remove leading dot if present
        if value_lower.startswith('.'):
            value_lower = value_lower[1:]
        
        # Extension mapping
        extension_map = {
            'gz': cls.GZIP,
            'gzip': cls.GZIP,
            'bz2': cls.BZIP2,
            'bzip2': cls.BZIP2,
            'tgz': cls.TGZ,
            'tar.gz': cls.TGZ,
        }
        
        # Try to find matching enum
        if value_lower in extension_map:
            return extension_map[value_lower]
        
        # Check if it's a filename
        if '.' in value_lower:
            path = Path(value)
            filename = path.name.lower()
            
            # Check compound extensions
            if filename.endswith('.tar.gz'):
                return cls.TGZ
            
            # Check single extension
            ext = path.suffix.lower()
            if ext:
                return cls._missing_(ext)  # Recursive call with just extension
        
        # Default to NONE if nothing matches
        return cls.NONE


class GenomeFormat(Enum):
    FASTA = "fasta"
    GENBANK = "genbank"
    
    @classmethod
    def _missing_(cls, value):
        """
        Handle GenomeFormat('fasta'), GenomeFormat('.fa'), GenomeFormat('file.fasta')
        """
        value_lower = str(value).lower().strip()
        
        # Remove leading dot
        if value_lower.startswith('.'):
            value_lower = value_lower[1:]
        
        # Extension mapping
        extension_map = {
            'fa': cls.FASTA,
            'fasta': cls.FASTA,
            'fna': cls.FASTA,
            'faa': cls.FASTA,
            'genbank': cls.GENBANK,
            'gb': cls.GENBANK,
            'gbk': cls.GENBANK,
        }
        
        # Direct match
        if value_lower in extension_map:
            return extension_map[value_lower]
        
        # If it looks like a filename, extract extension
        if '.' in value_lower:
            ext = Path(value).suffix.lower()[1:]  # Remove dot
            if ext in extension_map:
                return extension_map[ext]
        
        raise ValueError(f"'{value}' is not a valid {cls.__name__}")
