"""
Genome file validator and processor.

Handles FASTA and GenBank formats with compression support.
"""

from pathlib import Path
from typing import Optional, Dict, Any, List
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import bz2
import shutil

from validation_pkg.logger import get_logger
from validation_pkg.exceptions import (
    GenomeValidationError,
    FileFormatError,
    FastaFormatError,
    GenBankFormatError,
    CompressionError,
    FileNotFoundError as ValidationFileNotFoundError
)


class GenomeValidator:
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
    
    # Supported file extensions
    FASTA_EXTENSIONS = ['.fa', '.fasta', '.fna']
    GENBANK_EXTENSIONS = ['.gb', '.gbk', '.genbank']
    COMPRESSION_EXTENSIONS = ['.gz', '.bz2']
    
    def __init__(self, genome_config, output_dir, settings: dict):
        """
        Initialize genome validator.
        
        Args:
            genome_config: GenomeConfig object
            config_dir: Directory containing the config file (for resolving paths)
            output_dir: Directory for output files
            settings: specifying the edit processes
        """
        self.logger = get_logger()
        self.genome_config = genome_config
        self.output_dir = output_dir
        self.settings = settings
        
        # Resolved paths
        self.input_path = genome_config.filepath
        
        # File properties (discovered during processing)
        self.is_compressed = False
        self.compression_type = None  # 'gz' or 'bz2'
        self.detected_format = None  # 'fasta' or 'genbank'
        
        # Parsed data
        self.sequences = []  # List of SeqRecord objects
        
        # Statistics
        self.statistics = {
            'num_sequences': 0,
            'total_length': 0,
            'gc_content': 0.0,
            'min_length': 0,
            'max_length': 0,
            'avg_length': 0.0
        }
    
    def validate(self):
        """
        Main validation and processing workflow.
        
        Raises:
            GenomeValidationError: If validation fails
        """
        self.logger.info(f"Processing genome file: {self.genome_config.filename}")
        
        try:
            # Step 1: Check file existence
            self._check_file_exists()
            
            # Step 2: Detect compression
            self._detect_compression()
            
            # Step 3: Detect file format
            self._detect_format()
            
            # Step 4: Parse and validate
            self._parse_file()
            
            # Step 5: Apply editing specifications
            self._apply_edits()
            
            # Step 6: Collect statistics
            self._collect_statistics()
            
            # Step 7: Convert to FASTA (if needed)
            self._convert_to_fasta()
            
            # Step 8: Save to output directory
            output_path = self._save_output()
            
            self.logger.info(f"✓ Genome validation completed: {output_path.name}")
            self.logger.debug(f"Statistics: {self.statistics}")
            
        except Exception as e:
            self.logger.error(f"Genome validation failed: {e}")
            raise
    
    def _check_file_exists(self):
        """Check if input file exists."""
        if not self.input_path.exists():
            error_msg = f"Genome file not found: {self.genome_config.filename}"
            self.logger.add_validation_issue(
                level='ERROR',
                category='genome',
                message=error_msg,
                details={'file': str(self.input_path)}
            )
            raise ValidationFileNotFoundError(error_msg)
        
        self.logger.debug(f"File exists: {self.input_path}")
    
    def _detect_compression(self):
        """Detect if file is compressed."""
        suffix = self.input_path.suffix.lower()
        
        if suffix == '.gz':
            self.is_compressed = True
            self.compression_type = 'gz'
            self.logger.debug("Detected gzip compression")
        elif suffix == '.bz2':
            self.is_compressed = True
            self.compression_type = 'bz2'
            self.logger.debug("Detected bzip2 compression")
        else:
            self.is_compressed = False
            self.logger.debug("No compression detected")
    
    def _detect_format(self):
        """Detect file format based on extension."""
        # Get filename without compression extension
        filename = str(self.input_path)
        if self.is_compressed:
            # Remove .gz or .bz2
            filename = filename.rsplit('.', 1)[0]
        
        filename_lower = filename.lower()
        
        # Check FASTA extensions
        for ext in self.FASTA_EXTENSIONS:
            if filename_lower.endswith(ext):
                self.detected_format = 'fasta'
                self.logger.debug(f"Detected FASTA format (extension: {ext})")
                return
        
        # Check GenBank extensions
        for ext in self.GENBANK_EXTENSIONS:
            if filename_lower.endswith(ext):
                self.detected_format = 'genbank'
                self.logger.debug(f"Detected GenBank format (extension: {ext})")
                return
        
        # Unknown format
        error_msg = f"Unknown file format: {self.input_path.name}"
        self.logger.add_validation_issue(
            level='ERROR',
            category='genome',
            message=error_msg,
            details={'file': self.genome_config.filename, 'supported': 'FASTA or GenBank'}
        )
        raise FileFormatError(error_msg)
    
    def _open_file(self, mode='rt'):
        """
        Open file with automatic decompression.
        
        Args:
            mode: File opening mode (default: 'rt' for text read)
            
        Returns:
            File handle
        """
        try:
            if self.compression_type == 'gz':
                return gzip.open(self.input_path, mode)
            elif self.compression_type == 'bz2':
                return bz2.open(self.input_path, mode)
            else:
                return open(self.input_path, mode)
        except Exception as e:
            error_msg = f"Failed to open file: {e}"
            self.logger.add_validation_issue(
                level='ERROR',
                category='genome',
                message=error_msg,
                details={'file': str(self.input_path), 'error': str(e)}
            )
            raise CompressionError(error_msg) from e
    
    def _parse_file(self):
        """Parse file using BioPython and validate format."""
        self.logger.debug(f"Parsing {self.detected_format} file...")
        
        try:
            with self._open_file() as handle:
                # Parse using BioPython
                self.sequences = list(SeqIO.parse(handle, self.detected_format))
            
            # Validate we got sequences
            if not self.sequences:
                error_msg = f"No sequences found in {self.detected_format} file"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='genome',
                    message=error_msg,
                    details={'file': self.genome_config.filename, 'format': self.detected_format}
                )
                raise GenomeValidationError(error_msg)
            
            self.logger.debug(f"Parsed {len(self.sequences)} sequence(s)")
            
            # Validate each sequence
            self._validate_sequences()
            
        except FileFormatError:
            raise
        except Exception as e:
            error_msg = f"Failed to parse {self.detected_format} file: {e}"
            
            if self.detected_format == 'fasta':
                exception_class = FastaFormatError
            elif self.detected_format == 'genbank':
                exception_class = GenBankFormatError
            else:
                exception_class = FileFormatError
            
            self.logger.add_validation_issue(
                level='ERROR',
                category='genome',
                message=error_msg,
                details={
                    'file': self.genome_config.filename,
                    'format': self.detected_format,
                    'error': str(e)
                }
            )
            raise exception_class(error_msg) from e
    
    def _validate_sequences(self):
        """Validate parsed sequences."""
        self.logger.debug("Validating sequences...")
        
        for idx, record in enumerate(self.sequences):
            # Check sequence ID
            if not record.id:
                self.logger.add_validation_issue(
                    level='WARNING',
                    category='genome',
                    message=f'Sequence {idx} has no ID',
                    details={'sequence_index': idx}
                )
            
            # Check sequence length
            if len(record.seq) == 0:
                error_msg = f"Sequence '{record.id}' has zero length"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='genome',
                    message=error_msg,
                    details={'sequence_id': record.id, 'index': idx}
                )
                raise GenomeValidationError(error_msg)
            
            # Check for valid nucleotides
            valid_chars = set('ATCGNatcgn')
            seq_str = str(record.seq)
            invalid_chars = set(seq_str) - valid_chars
            
            if invalid_chars:
                self.logger.add_validation_issue(
                    level='WARNING',
                    category='genome',
                    message=f"Sequence '{record.id}' contains non-standard characters",
                    details={
                        'sequence_id': record.id,
                        'invalid_chars': list(invalid_chars),
                        'count': sum(1 for c in seq_str if c in invalid_chars)
                    }
                )
        
        self.logger.debug("✓ Sequence validation passed")
    
    def _apply_edits(self):
        """
        Apply editing specifications to sequences.
        
        This is where you can add custom editing logic.
        Examples:
        - Remove sequences below certain length
        - Mask low-complexity regions
        - Trim sequences
        - Replace ambiguous nucleotides
        """
        self.logger.debug("Applying editing specifications...")
        
        # Example 1: Remove very short sequences
        min_length = 100
        original_count = len(self.sequences)
        self.sequences = [seq for seq in self.sequences if len(seq.seq) >= min_length]
        
        if len(self.sequences) < original_count:
            removed = original_count - len(self.sequences)
            self.logger.add_validation_issue(
                level='WARNING',
                category='genome',
                message=f'Removed {removed} sequence(s) shorter than {min_length}bp',
                details={'min_length': min_length, 'removed_count': removed}
            )
        
        # Example 2: Convert to uppercase
        for record in self.sequences:
            record.seq = record.seq.upper()
        
        # Example 3: Check for duplicate IDs
        seq_ids = [record.id for record in self.sequences]
        if len(seq_ids) != len(set(seq_ids)):
            duplicates = [sid for sid in seq_ids if seq_ids.count(sid) > 1]
            self.logger.add_validation_issue(
                level='WARNING',
                category='genome',
                message=f'Duplicate sequence IDs found',
                details={'duplicate_ids': list(set(duplicates))}
            )
        
        # TODO: Add more editing specifications as needed
        # - self._remove_ambiguous_nucleotides()
        # - self._trim_sequences()
        # - self._mask_repeats()
        
        self.logger.debug(f"✓ Edits applied, {len(self.sequences)} sequence(s) remaining")
    
    def _collect_statistics(self):
        """Collect statistics about the sequences."""
        self.logger.debug("Collecting statistics...")
        
        if not self.sequences:
            return
        
        lengths = [len(record.seq) for record in self.sequences]
        
        self.statistics = {
            'num_sequences': len(self.sequences),
            'total_length': sum(lengths),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'avg_length': sum(lengths) / len(lengths)
        }
        
        # Calculate GC content
        total_gc = 0
        total_bases = 0
        
        for record in self.sequences:
            seq_str = str(record.seq).upper()
            total_gc += seq_str.count('G') + seq_str.count('C')
            total_bases += len(seq_str)
        
        if total_bases > 0:
            self.statistics['gc_content'] = (total_gc / total_bases) * 100
        
        # Log unusual statistics
        if self.statistics['gc_content'] < 20 or self.statistics['gc_content'] > 80:
            self.logger.add_validation_issue(
                level='WARNING',
                category='genome',
                message=f"Unusual GC content: {self.statistics['gc_content']:.2f}%",
                details={'gc_content': self.statistics['gc_content']}
            )
        
        self.logger.debug(f"Statistics: {self.statistics}")
    
    def _convert_to_fasta(self):
        """Convert sequences to FASTA format (if not already)."""
        if self.detected_format == 'fasta':
            self.logger.debug("Already in FASTA format")
            return
        
        self.logger.debug(f"Converting from {self.detected_format} to FASTA...")
        
        # BioPython SeqRecord objects can be written as FASTA directly
        # No conversion needed, just change the output format
        
        self.logger.debug("✓ Ready for FASTA output")
    
    def _save_output(self) -> Path:
        """
        Save processed genome to output directory.
        
        Returns:
            Path to output file
        """
        self.logger.debug("Saving output file...")
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate output filename
        base_name = self.input_path.stem
        if self.is_compressed:
            # Remove compression extension from stem
            base_name = Path(base_name).stem
        
        # Always save as FASTA
        output_filename = f"{base_name}.fasta"
        
        # Add compression extension if requested
        if self.settings['output']['coding_type'] == 'gz':
            output_filename += '.gz'
        elif self.settings['output']['coding_type'] == 'bzip':
            output_filename += '.bz2'
        
        output_path = self.output_dir / output_filename
        
        # Write output
        if self.settings['output']['coding_type'] == 'gz':
            with gzip.open(output_path, 'wt') as handle:
                SeqIO.write(self.sequences, handle, 'fasta')
        elif self.settings['output']['coding_type'] == 'bzip':
            with bz2.open(output_path, 'wt') as handle:
                SeqIO.write(self.sequences, handle, 'fasta')
        else:
            with open(output_path, 'w') as handle:
                SeqIO.write(self.sequences, handle, 'fasta')
        
        self.logger.info(f"Output saved: {output_path}")
        
        return output_path
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        Get statistics about the processed genome.
        
        Returns:
            Dictionary of statistics
        """
        return self.statistics.copy()


# Convenience function for easy use
def validate_genome(genome_config, config_dir: Path, output_dir: Path) -> GenomeValidator:
    """
    Validate and process a genome file.
    
    Args:
        genome_config: GenomeConfig object
        config_dir: Directory containing config file
        output_dir: Output directory for processed files
        
    Returns:
        GenomeValidator instance with statistics
    """
    validator = GenomeValidator(genome_config, config_dir, output_dir)
    validator.validate()
    return validator