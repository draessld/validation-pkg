"""
Genome file validator and processor.

Handles FASTA and GenBank formats with compression support.
"""

from pathlib import Path
from typing import Optional, Dict, Any, List, IO
from dataclasses import dataclass
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import bz2
import shutil

from validation_pkg.logger import get_logger
from validation_pkg.utils.settings import BaseSettings
from validation_pkg.utils.formats import CodingType as CT
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
    2. Parse and validate using BioPython
    3. Apply editing specifications
    4. Convert to FASTA format
    5. Compress if requested
    6. Save to output directory
    """

    @dataclass
    class Settings(BaseSettings):
        """
        Settings for genome validation and processing.
        Attributes:
            allow_empty_sequences: Allow SeqRecord empty seequence (default: False) 
            allow_empty_id: Allow SeqRecord empty ID (default: False)
            warn_n_sequences: Warn if number of sequences exceeds this (default: 2)
            plasmid_split: Separate plasmid sequences into different file (default: True)
            replace_id_with: Prefix to add to sequence IDs, None = no prefix (default: None)
            min_sequence_length: Minimum sequence length to keep in bp, remove shorter (default: 100)
            coding_type: Output compression type: 'gz', 'bz2', or None (default: None)
            output_filename_suffix: Suffix to add to output filename (default: None)
            output_subdir_name: Subdirectory name for output files (default: None)
        """

        # Validation thresholds
        allow_empty_sequences: bool = False #   Raise ERROR
        allow_empty_id: bool = False    #   Raise ERROR
        warn_n_sequences: int = 2   #   Raise Warning - set plasmid_split to True

        # Editing specifications
        plasmid_split: bool = True 
        replace_id_with: Optional[str] = None 
        min_sequence_length: int = 100

        # Output format
        coding_type: Optional[str] = None
        output_filename_suffix: Optional[str] = None
        output_subdir_name: Optional[str] = None

    def __init__(self, genome_config, output_dir: Path, settings: Optional[Settings] = None) -> None:
        """
        Initialize genome validator.

        Args:
            genome_config: GenomeConfig object from ConfigManager with file info
            output_dir: Directory for output files (Path object)
            settings: Settings object with validation parameters (uses defaults if None)

        Note: In future versions, output_dir may be moved to settings.

        Example:
            >>> settings = GenomeValidator.Settings()
            >>> settings = settings.update(replace_id_with="chr1")
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
    
    def validate(self) -> None:
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

    def _open_file(self, mode: str = 'rt') -> IO:
        """
        Open file with automatic decompression based on genome_config.

        Uses centralized file opening from file_handler.py to eliminate code duplication.

        Args:
            mode: File opening mode (default: 'rt' for text read)

        Returns:
            File handle

        Raises:
            CompressionError: If file cannot be opened or decompressed
        """
        try:
            from validation_pkg.utils.file_handler import open_file_with_coding_type
            return open_file_with_coding_type(
                self.input_path,
                self.genome_config.coding_type,
                mode
            )
        except CompressionError as e:
            self.logger.add_validation_issue(
                level='ERROR',
                category='genome',
                message=str(e),
                details={'file': str(self.input_path)}
            )
            raise
    
    def _parse_file(self) -> None:
        """Parse file using BioPython and validate format from genome_config."""
        self.logger.debug(f"Parsing {self.genome_config.detected_format} file...")

        try:
            with self._open_file() as handle:
                # Parse using BioPython with clean enum conversion
                biopython_format = self.genome_config.detected_format.to_biopython()
                self.sequences = list(SeqIO.parse(handle, biopython_format))

            # Validate we got sequences
            if not self.sequences:
                error_msg = f"No sequences found in {self.genome_config.detected_format} file"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='genome',
                    message=error_msg,
                    details={'file': self.genome_config.filename, 'format': str(self.genome_config.detected_format)}
                )
                raise GenomeValidationError(error_msg)

            self.logger.debug(f"Parsed {len(self.sequences)} sequence(s)")

            # Validate sequences
            self._validate_sequences()

        except FileFormatError:
            raise
        except Exception as e:
            error_msg = f"Failed to parse {self.genome_config.detected_format} file: {e}"

            # Determine exception type based on format using enum
            from validation_pkg.utils.formats import GenomeFormat
            if self.genome_config.detected_format == GenomeFormat.FASTA:
                exception_class = FastaFormatError
            elif self.genome_config.detected_format == GenomeFormat.GENBANK:
                exception_class = GenBankFormatError
            else:
                exception_class = FileFormatError

            self.logger.add_validation_issue(
                level='ERROR',
                category='genome',
                message=error_msg,
                details={
                    'file': self.genome_config.filename,
                    'format': str(self.genome_config.detected_format),
                    'error': str(e)
                }
            )
            raise exception_class(error_msg) from e
    
    def _validate_sequences(self) -> None:
        """Validate parsed sequences."""
        self.logger.debug("Validating sequences...")
        
        # Warn about number of sequences
        if len(self.sequences) >= self.settings.warn_n_sequences:
            self.logger.add_validation_issue(
                level='WARNING',
                category='genome',
                message=f"High number of sequences: {len(self.sequences)}",
                details={
                    'num_sequences': len(self.sequences),
                    'threshold': self.settings.warn_n_sequences
                }
            )

        for idx, record in enumerate(self.sequences):
            # Check sequence ID
            if not record.id and not self.settings.allow_empty_id:
                error_msg = f"Sequence at index {idx} has no ID"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='genome',
                    message=error_msg,
                    details={'sequence_index': idx, 'sequence_id': str(record.id)}
                )
                raise GenomeValidationError(error_msg)
            
            # Check sequence length
            if len(record.seq) == 0 and not self.settings.allow_empty_sequences:
                error_msg = f"Sequence '{record.id}' has zero length"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='genome',
                    message=error_msg,
                    details={'sequence_id': record.id, 'index': idx}
                )
                raise GenomeValidationError(error_msg)
            
        
        self.logger.debug("✓ Sequence validation passed")
    
    def _apply_edits(self) -> None:
        """
        Apply editing specifications to sequences based on settings.

        Applies the following edits (if enabled in settings):
        - Remove short sequences (min_sequence_length)
        - Add sequence prefix (replace_id_with)
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
        if self.settings.replace_id_with:
            prefix = self.settings.replace_id_with
            for record in self.sequences:
                # Add prefix to sequence ID
                record.description = f"{record.id}"
                record.id = f"{prefix}"
            self.logger.debug(f"Added prefix '{prefix}' to all sequence IDs")

        # 3. Split plasmid sequences
        if self.settings.plasmid_split and len(self.sequences) > 1:
            self._split_plasmids()

        self.logger.debug(f"✓ Edits applied, {len(self.sequences)} sequence(s) remaining")

    def _split_plasmids(self) -> None:
        """
        Split plasmid sequences into a separate file.

        When more than 2 sequences are present:
        - Keep the longest sequence (assumed to be the chromosome) in main output
        - Save remaining sequences (plasmids) to a separate file with "_plasmid" suffix
        - Uses same output settings (compression, subdirectory) as main file
        """
        self.logger.debug(f"Splitting plasmids from {len(self.sequences)} sequences...")

        # Sort sequences by length (longest first)
        sorted_sequences = sorted(self.sequences, key=lambda x: len(x.seq), reverse=True)

        # Keep longest as main chromosome
        main_sequence = [sorted_sequences[0]]
        plasmid_sequences = sorted_sequences[1:]

        self.logger.info(
            f"Splitting genome: keeping longest sequence ({sorted_sequences[0].id}, "
            f"{len(sorted_sequences[0].seq)} bp) as main, "
            f"saving {len(plasmid_sequences)} sequence(s) as plasmids"
        )

        for i,plasmid in enumerate(plasmid_sequences):
            # Save plasmid sequences to separate file
            self._save_plasmid_file(plasmid_sequences,i)

        # Update main sequences to only include chromosome
        self.sequences = main_sequence

    def _save_plasmid_file(self, plasmid_sequences: SeqRecord, index: int) -> None:
        """
        Save plasmid sequences to a separate file.

        Args:
            plasmid_sequences: List of SeqRecord objects for plasmids
        """
        # Determine output directory (with optional subdirectory)
        output_dir = self.output_dir
        if self.settings.output_subdir_name:
            output_dir = output_dir / self.settings.output_subdir_name

        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate plasmid filename
        base_name = self.genome_config.filename
        # Remove all suffixes
        for suffix in self.input_path.suffixes:
            base_name = base_name.replace(suffix, '')

        # Add plasmid suffix and optional user suffix
        if self.settings.output_filename_suffix:
            plasmid_filename = f"{base_name}_{self.settings.output_filename_suffix}_plasmid{index}.fasta"
        else:
            plasmid_filename = f"{base_name}_plasmid{index}.fasta"

        # Add compression extension if requested
        coding = self.settings.coding_type

        if coding in ('gz', 'gzip', CT.GZIP):
            plasmid_filename += '.gz'
        elif coding in ('bz2', 'bzip2', CT.BZIP2):
            plasmid_filename += '.bz2'
            
        plasmid_path = output_dir / plasmid_filename

        # Write plasmid sequences with appropriate compression
        self.logger.debug(f"Writing plasmid sequences to: {plasmid_path}")


        if coding in ('gz', 'gzip', CT.GZIP):
            with gzip.open(plasmid_path, 'wt') as handle:
                SeqIO.write(plasmid_sequences, handle, 'fasta')
        elif coding in ('bz2', 'bzip2', CT.BZIP2):
            with bz2.open(plasmid_path, 'wt') as handle:
                SeqIO.write(plasmid_sequences, handle, 'fasta')
        else:
            with open(plasmid_path, 'w') as handle:
                SeqIO.write(plasmid_sequences, handle, 'fasta')

        self.logger.info(f"Plasmid sequences saved: {plasmid_path}")

        # Log details about each plasmid
        for seq in plasmid_sequences:
            self.logger.debug(f"  Plasmid: {seq.id} ({len(seq.seq)} bp)")

    def _collect_statistics(self) -> None:
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

        self.logger.debug(f"Statistics: {self.statistics}")
    
    def _convert_to_fasta(self) -> None:
        """Convert sequences to FASTA format (if not already)."""
        # Check if already FASTA using enum comparison
        from validation_pkg.utils.formats import GenomeFormat

        if self.genome_config.detected_format == GenomeFormat.FASTA:
            self.logger.debug("Already in FASTA format")
            return

        self.logger.debug(f"Converting from {self.genome_config.detected_format} to FASTA...")

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
        # Support both string ('gz', 'bz2') and enum (CodingType.GZIP, CodingType.BZIP2)
        coding = self.settings.coding_type

        if coding in ('gz', 'gzip', CT.GZIP):
            output_filename += '.gz'
        elif coding in ('bz2', 'bzip2', CT.BZIP2):
            output_filename += '.bz2'

        output_path = output_dir / output_filename

        # Write output with appropriate compression
        self.logger.debug(f"Writing output to: {output_path}")

        if coding in ('gz', 'gzip', CT.GZIP):
            with gzip.open(output_path, 'wt') as handle:
                SeqIO.write(self.sequences, handle, 'fasta')
        elif coding in ('bz2', 'bzip2', CT.BZIP2):
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
