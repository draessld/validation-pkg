from enum import Enum
from pathlib import Path

# ===== Enums  =====

class CodingType(Enum):
    """
    Supported compression types for genomic data files.

    Supported formats:
    - GZIP: .gz files (gzip compression)
    - BZIP2: .bz2 files (bzip2 compression)
    - NONE: Uncompressed files

    Note: TAR archives (.tar.gz, .tgz) are NOT supported.
    Please extract archives before processing.
    """
    GZIP = "gzip"
    BZIP2 = "bzip2"
    NONE = "none"

    @classmethod
    def _missing_(cls, value):
        """
        Called when enum lookup fails. Allows flexible input formats.

        Supports:
        - CodingType('gz') or CodingType('.gz') → GZIP
        - CodingType('bz2') or CodingType('.bz2') → BZIP2
        - CodingType('file.fasta.gz') → GZIP (extracts from filename)
        - CodingType('') or invalid → NONE

        Args:
            value: String value to convert (extension, filename, or format name)

        Returns:
            CodingType enum member
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
        }

        # Try to find matching enum
        if value_lower in extension_map:
            return extension_map[value_lower]

        # Check if it's a filename with extension
        if '.' in value_lower:
            path = Path(value)
            ext = path.suffix.lower()
            if ext:
                # Remove dot and try again
                return cls._missing_(ext)

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


class ReadFormat(Enum):
    FASTQ = "fastq"
    BAM = "bam"
    
    @classmethod
    def _missing_(cls, value):
        """
        Handle ReadFormat('fastq'), ReadFormat('.fq'), ReadFormat('file.fastq')
        """
        value_lower = str(value).lower().strip()
        
        # Remove leading dot
        if value_lower.startswith('.'):
            value_lower = value_lower[1:]
        
        # Extension mapping
        extension_map = {
            'fq': cls.FASTQ,
            'fastq': cls.FASTQ,
            'bam': cls.BAM,
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


class FeatureFormat(Enum):
    GFF = "gff"
    BED = "bed"
    
    @classmethod
    def _missing_(cls, value):
        """
        Handle ReadFormat('bed'), ReadFormat('.gtf'), ReadFormat('file.gff')
        """
        value_lower = str(value).lower().strip()
        
        # Remove leading dot
        if value_lower.startswith('.'):
            value_lower = value_lower[1:]
        
        # Extension mapping
        extension_map = {
            'gff': cls.GFF,
            'gff3': cls.GFF,
            'gtf': cls.GFF,
            'gff2': cls.GFF,
            'bed': cls.BED,
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