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
    GenomeValidationError,
    FileFormatError,
    FastaFormatError,
    GenBankFormatError,
    CompressionError
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

    @dataclass
    class Settings(BaseSettings):
        """
        Settings for genome validation and processing.

        Attributes:
            plasmid_split: Separate plasmid sequences into different file (default: True)
            sequence_prefix: Prefix to add to sequence IDs, None = no prefix (default: None)
            min_sequence_length: Minimum sequence length to keep in bp (default: 100)
            convert_to_uppercase: Convert all sequences to uppercase (default: False)
            check_duplicate_ids: Check for and warn about duplicate sequence IDs (default: True)
            warn_n_sequences: Warn if number of sequences exceeds this (default: 2)
            allow_duplicate_ids: Allow duplicate sequence IDs, just warn (default: True)
            always_fasta: Always output as FASTA format (default: True)
            coding_type: Output compression type: 'gz', 'bz2', or None (default: None)
            output_filename_suffix: Suffix to add to output filename (default: None)
            output_subdir_name: Subdirectory name for output files (default: None)
            line_width: Characters per line in FASTA output (default: 80)

        Example:
            >>> settings = GenomeValidator.Settings()
            >>> print(settings)  # View all settings
            >>> ref_settings = settings.update(sequence_prefix="chr1", output_filename_suffix="ref")
            >>> validator = GenomeValidator(genome_config, output_dir, ref_settings)
        """
        # Editing specifications
        plasmid_split: bool = True
        sequence_prefix: Optional[str] = None
        min_sequence_length: int = 100
        check_duplicate_ids: bool = True

        # Validation thresholds
        warn_n_sequences: int = 2
        allow_duplicate_ids: bool = True

        # Output format
        always_fasta: bool = True
        coding_type: Optional[str] = None
        output_filename_suffix: Optional[str] = None
        output_subdir_name: Optional[str] = None
        line_width: int = 80

    # Supported file extensions
    FASTA_EXTENSIONS = ['.fa', '.fasta', '.fna']
    GENBANK_EXTENSIONS = ['.gb', '.gbk', '.genbank']
    COMPRESSION_EXTENSIONS = ['.gz', '.bz2']

    def __init__(self, genome_config, output_dir, settings: Optional[Settings] = None):
        """
        Initialize genome validator.

        Args:
            genome_config: GenomeConfig object from ConfigManager with file info
            output_dir: Directory for output files (Path object)
            settings: Settings object with validation parameters (uses defaults if None)

        Note: In future versions, output_dir may be moved to settings.

        Example:
            >>> settings = GenomeValidator.Settings()
            >>> settings = settings.update(sequence_prefix="chr1")
            >>> validator = GenomeValidator(genome_config, output_dir, settings)
        """
        self.logger = get_logger()
        self.genome_config = genome_config
        self.output_dir = Path(output_dir)
        self.settings = settings if settings is not None else self.Settings()

        # Log settings being used
        self.logger.debug(f"Initializing GenomeValidator with settings:\n{self.settings}")

        # Resolved paths
        self.input_path = genome_config.filepath

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

        Uses genome_config data (format, compression) provided by ConfigManager.

        Raises:
            GenomeValidationError: If validation fails
        """
        self.logger.info(f"Processing genome file: {self.genome_config.filename}")
        self.logger.debug(f"Format: {self.genome_config.detected_format}, Compression: {self.genome_config.coding_type}")

        try:
            # File existence already checked by ConfigManager
            # Format and compression already detected in genome_config

            # Step 1: Parse and validate
            self._parse_file()

            # Step 2: Apply editing specifications
            self._apply_edits()

            # Step 3: Collect statistics
            self._collect_statistics()

            # Step 4: Convert to FASTA (if needed)
            self._convert_to_fasta()

            # Step 5: Save to output directory
            output_path = self._save_output()

            self.logger.info(f"✓ Genome validation completed: {output_path.name}")
            self.logger.debug(f"Statistics: {self.statistics}")

        except Exception as e:
            self.logger.error(f"Genome validation failed: {e}")
            raise

    def _open_file(self, mode='rt'):
        """
        Open file with automatic decompression based on genome_config.

        Note: File compression handling will be moved to utils/file_handler.py in future.

        Args:
            mode: File opening mode (default: 'rt' for text read)

        Returns:
            File handle

        Raises:
            CompressionError: If file cannot be opened or decompressed
        """
        try:
            coding_type_str = str(self.genome_config.coding_type)

            if coding_type_str == 'gz':
                return gzip.open(self.input_path, mode)
            elif coding_type_str == 'bz2':
                return bz2.open(self.input_path, mode)
            elif coding_type_str == '' or coding_type_str == 'CodingType()':
                # No compression
                return open(self.input_path, mode)
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
        """Parse file using BioPython and validate format from genome_config."""
        format_str = str(self.genome_config.detected_format)
        self.logger.debug(f"Parsing {format_str} file...")

        try:
            with self._open_file() as handle:
                # Parse using BioPython - convert GenomeFormat to BioPython format string
                biopython_format = format_str.lower().replace('genomeformat.', '')
                self.sequences = list(SeqIO.parse(handle, biopython_format))

            # Validate we got sequences
            if not self.sequences:
                error_msg = f"No sequences found in {format_str} file"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='genome',
                    message=error_msg,
                    details={'file': self.genome_config.filename, 'format': format_str}
                )
                raise GenomeValidationError(error_msg)

            self.logger.debug(f"Parsed {len(self.sequences)} sequence(s)")

            # Validate each sequence
            self._validate_sequences()

        except FileFormatError:
            raise
        except Exception as e:
            error_msg = f"Failed to parse {format_str} file: {e}"

            if 'fasta' in format_str.lower():
                exception_class = FastaFormatError
            elif 'genbank' in format_str.lower() or 'gb' in format_str.lower():
                exception_class = GenBankFormatError
            else:
                exception_class = FileFormatError

            self.logger.add_validation_issue(
                level='ERROR',
                category='genome',
                message=error_msg,
                details={
                    'file': self.genome_config.filename,
                    'format': format_str,
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
        Apply editing specifications to sequences based on settings.

        Applies the following edits (if enabled in settings):
        - Remove short sequences (min_sequence_length)
        - Add sequence prefix (sequence_prefix)
        - Check for duplicate IDs (check_duplicate_ids)
        """
        self.logger.debug("Applying editing specifications from settings...")

        # 1. Remove short sequences
        if self.settings.min_sequence_length > 0:
            min_length = self.settings.min_sequence_length
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

        # 2. Add sequence prefix
        if self.settings.sequence_prefix:
            prefix = self.settings.sequence_prefix
            for record in self.sequences:
                # Add prefix to sequence ID
                record.id = f"{prefix}_{record.id}"
                # Also update the name field to match
                record.name = record.id
            self.logger.debug(f"Added prefix '{prefix}' to all sequence IDs")

        # 3. Check for duplicate IDs
        if self.settings.check_duplicate_ids:
            seq_ids = [record.id for record in self.sequences]
            if len(seq_ids) != len(set(seq_ids)):
                duplicates = [sid for sid in seq_ids if seq_ids.count(sid) > 1]
                level = 'ERROR' if not self.settings.allow_duplicate_ids else 'WARNING'
                self.logger.add_validation_issue(
                    level=level,
                    category='genome',
                    message=f'Duplicate sequence IDs found',
                    details={'duplicate_ids': list(set(duplicates))}
                )
                if not self.settings.allow_duplicate_ids:
                    raise GenomeValidationError(f"Duplicate sequence IDs not allowed: {list(set(duplicates))}")

        self.logger.debug(f"✓ Edits applied, {len(self.sequences)} sequence(s) remaining")
    
    def _collect_statistics(self):
        """Collect statistics about the sequences and validate against thresholds."""
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

        # Warn about number of sequences
        if len(self.sequences) > self.settings.warn_n_sequences:
            self.logger.add_validation_issue(
                level='WARNING',
                category='genome',
                message=f"High number of sequences: {len(self.sequences)}",
                details={
                    'num_sequences': len(self.sequences),
                    'threshold': self.settings.warn_n_sequences
                }
            )

        self.logger.debug(f"Statistics: {self.statistics}")
    
    def _convert_to_fasta(self):
        """Convert sequences to FASTA format (if not already)."""
        format_str = str(self.genome_config.detected_format).lower()

        if 'fasta' in format_str:
            self.logger.debug("Already in FASTA format")
            return

        self.logger.debug(f"Converting from {format_str} to FASTA...")

        # BioPython SeqRecord objects can be written as FASTA directly
        # No conversion needed, just change the output format

        self.logger.debug("✓ Ready for FASTA output")
    
    def _save_output(self) -> Path:
        """
        Save processed genome to output directory using settings.

        Note: Compression handling will be moved to utils/file_handler.py in future.

        Returns:
            Path to output file
        """
        self.logger.debug("Saving output file...")

        # Determine output directory (with optional subdirectory)
        output_dir = self.output_dir
        if self.settings.output_subdir_name:
            output_dir = output_dir / self.settings.output_subdir_name

        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate output filename from original filename (without compression extension)
        # Get base name without any extensions
        base_name = self.genome_config.filename
        # Remove all suffixes (.fasta.gz -> remove .gz then .fasta)
        for suffix in self.input_path.suffixes:
            base_name = base_name.replace(suffix, '')

        # Add suffix if specified
        if self.settings.output_filename_suffix:
            output_filename = f"{base_name}_{self.settings.output_filename_suffix}.fasta"
        else:
            output_filename = f"{base_name}.fasta"

        # Add compression extension if requested (from settings, not input)
        if self.settings.coding_type == 'gz':
            output_filename += '.gz'
        elif self.settings.coding_type == 'bz2':
            output_filename += '.bz2'

        output_path = output_dir / output_filename

        # Write output with appropriate compression
        self.logger.debug(f"Writing output to: {output_path}")

        if self.settings.coding_type == 'gz':
            with gzip.open(output_path, 'wt') as handle:
                SeqIO.write(self.sequences, handle, 'fasta')
        elif self.settings.coding_type == 'bz2':
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

