"""
Genome file validator and processor.

Handles FASTA and GenBank formats with compression support.
"""

from pathlib import Path
from typing import Optional, Dict, Any, List
from dataclasses import dataclass
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import bz2
import shutil

from validation_pkg.logger import get_logger
from validation_pkg.utils.settings import BaseSettings
from validation_pkg.exceptions import (
    FeatureValidationError,
    FileFormatError,
    BedFormatError,
    GffFormatError,
    CompressionError,
    FileNotFoundError as ValidationFileNotFoundError
)

class FeatureValidator:
    """
    Validates and processes genome files (FASTA and GenBank formats).

    Workflow:
    1. Detect and handle compression
    2. Detect file format
    3. Parse and validate using BioPython
    4. Apply editing specifications
    5. Convert to FASTA format
    6. Compress if requested
    7. Save to output directory
    """

    @dataclass
    class Settings(BaseSettings):
        """
        TODO: desribe
        """
        # Output format
        always_gff: bool = True
        coding_type: Optional[str] = None
        output_filename_suffix: Optional[str] = None
        output_subdir_name: Optional[str] = None

    # Supported file extensions
    GFF_EXTENSIONS = ['.gff', '.gff3', '.gtf', '.gff2']
    BED_EXTENSIONS = ['.bed',]
    COMPRESSION_EXTENSIONS = ['.gz', '.bz2']


    def __init__(self, feature_config, output_dir):
        """
        TODO: add description
        Initialize genome validator.

        Args:
            feature_config: GenomeConfig object
            output_dir: Directory for output files (Path object)
        """
        self.logger = get_logger()
        self.feature_config = feature_config
        self.output_dir = Path(output_dir)
        self.settings = self.Settings()

        # Log settings being used
        self.logger.debug(f"Initializing FeatureValidator with settings:\n{self.settings}")
        
        # Resolved paths
        self.input_path = feature_config.filepath
        
        # Parsed data
        self.data = []  # data from features
        
        # Statistics
        self.statistics = {
            #   TODO: think about later
        }
    
    def validate(self):
        """
        Main validation and processing workflow.
        
        Raises:
            FeatureValidationError: If validation fails
        """
        self.logger.info(f"Processing genome file: {self.feature_config.filename}")
        
        try:
            # Step 4: Parse and validate
            self._parse_file()
            
            # Step 5: Apply editing specifications
            self._apply_edits()
            
            # Step 6: Collect statistics
            self._collect_statistics()
            
            # Step 7: Convert to gff (if needed)
            self._convert_to_gff()
            
            # Step 8: Save to output directory
            output_path = self._save_output()
            
            self.logger.info(f"✓ Feature validation completed: {output_path.name}")
            self.logger.debug(f"Statistics: {self.statistics}")
            
        except Exception as e:
            self.logger.error(f"Feature validation failed: {e}")
            raise

    def _parse_file(self):
        pass

    def _apply_edits(self):
        pass

    def _collect_statistics(self):
        pass

    def _convert_to_gff(self):
        pass

    def _save_output(self):
        pass