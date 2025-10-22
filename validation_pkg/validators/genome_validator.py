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
from validation_pkg.utils.formats import CodingType as CT, GenomeFormat
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
            is_plasmid: Treat all sequences as plasmids (no main chromosome) (default: False)
            plasmid_split: Separate plasmid sequences into different files (default: False)
            plasmids_to_one: Merge all plasmid sequences into one file (default: False)
            main_longest: Select longest sequence as main chromosome (default: True)
            main_first: Select first sequence as main chromosome (default: False)
            replace_id_with: Prefix to add to sequence IDs, None = no prefix (default: None)
            min_sequence_length: Minimum sequence length to keep in bp, remove shorter (default: 100)
            coding_type: Output compression type: 'gz', 'bz2', or None (default: None)
            output_filename_suffix: Suffix to add to output filename (default: None)
            output_subdir_name: Subdirectory name for output files (default: None)

        Note:
            - plasmid_split and plasmids_to_one are mutually exclusive - only one can be True.
            - main_longest and main_first are mutually exclusive - only one can be True.
            - When is_plasmid=True, all sequences are treated as plasmids:
              - If plasmid_split=True: each sequence saved to separate file
              - If plasmids_to_one=True: all sequences saved to one merged file
              - If both False: all sequences saved to main output file
            - Main selection (main_longest/main_first) only applies when is_plasmid=False
        """
        validation_level: str = 'strict'  # 'strict', 'trust', or 'minimal'

        # Validation thresholds
        allow_empty_sequences: bool = False #   Raise ERROR
        allow_empty_id: bool = False    #   Raise ERROR
        warn_n_sequences: int = 2   #   Raise Warning - set plasmid_split to True

        # Editing specifications
        #   TODO: plasmid handeling better
        is_plasmid: bool = False
        plasmid_split: bool = False
        plasmids_to_one: bool = False
        main_longest: bool = True
        main_first: bool = False

        replace_id_with: Optional[str] = None
        min_sequence_length: int = 100

        # Output format
        coding_type: Optional[str] = None 
        output_filename_suffix: Optional[str] = None
        output_subdir_name: Optional[str] = None

        def __post_init__(self):
            """Validate settings after initialization."""
            if self.plasmid_split and self.plasmids_to_one:
                raise ValueError(
                    "plasmid_split and plasmids_to_one cannot both be True. "
                    "Choose one: plasmid_split creates separate files for each plasmid, "
                    "plasmids_to_one merges all plasmids into a single file."
                )

            if self.main_longest and self.main_first:
                raise ValueError(
                    "main_longest and main_first cannot both be True. "
                    "Choose one: main_longest selects the longest sequence as main, "
                    "main_first selects the first sequence as main."
                )

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
        self.output_dir = output_dir
        self.settings = settings if settings is not None else self.Settings()

        # Validate validation_level setting
        valid_levels = {'strict', 'trust', 'minimal'}
        if self.settings.validation_level not in valid_levels:
            raise ValueError(
                f"Invalid validation_level '{self.settings.validation_level}'. "
                f"Must be one of: {', '.join(valid_levels)}"
            )

        # Log settings being used
        self.logger.debug(f"Initializing GenomeValidator with settings:\n{self.settings}")

        # Resolved paths
        self.input_path = genome_config.filepath

        # Parsed data
        self.sequences = []  # List of SeqRecord objects

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
            self._parse_file()

            self._validate_sequences()

            self._convert_to_fasta()

            self._apply_edits() # include plasmid handle

            self._save_output()

            self.logger.info(f"✓ Genome validation completed")

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
        """
        Parse file using BioPython and validate format from genome_config.

        Behavior depends on validation_level:
        - 'strict': Parse all sequences, validate all
        - 'trust': Parse all sequences, validate only first one
        - 'minimal': Skip parsing entirely
        """
        self.logger.debug(f"Parsing {self.genome_config.detected_format} file (validation_level={self.settings.validation_level})...")

        # Minimal mode - skip parsing
        if self.settings.validation_level == 'minimal':
            self.logger.info("Minimal validation mode - skipping file parsing")
            self.sequences = []
            return

        # Trust and Strict modes - parse all sequences
        try:
            with self._open_file() as handle:
                # Parse using BioPython with clean enum conversion
                biopython_format = self.genome_config.detected_format.to_biopython()
                self.sequences = list(SeqIO.parse(handle, biopython_format))

            # Validate sequences
            if not self.sequences:
                error_msg = f"No sequences found in {self.genome_config.detected_format} file"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='genome',
                    message=error_msg,
                    details={'file': self.genome_config.filename, 'format': str(self.genome_config.detected_format)}
                )
                raise GenomeValidationError(error_msg)

            if self.settings.validation_level == 'trust':
                self.logger.info(f"Trust mode - parsed {len(self.sequences)} sequence(s), will validate first only")
            else:
                self.logger.debug(f"Parsed {len(self.sequences)} sequence(s)")

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
        """
        Validate parsed sequences.

        Behavior depends on validation_level:
        - 'strict': Validate all sequences
        - 'trust': Validate only first sequence
        - 'minimal': Skip validation (no parsing done)
        """
        self.logger.debug("Validating sequences...")

        # Minimal mode - no sequences to validate
        if self.settings.validation_level == 'minimal':
            self.logger.debug("Minimal mode - skipping sequence validation")
            return

        # Trust mode - validate only first sequence
        if self.settings.validation_level == 'trust':
            self.logger.debug("Trust mode - validating first sequence only")
            if len(self.sequences) > 0:
                record = self.sequences[0]
                # Check sequence ID
                if not record.id and not self.settings.allow_empty_id:
                    error_msg = f"First sequence has no ID"
                    self.logger.add_validation_issue(
                        level='ERROR',
                        category='genome',
                        message=error_msg,
                        details={'sequence_index': 0, 'sequence_id': str(record.id)}
                    )
                    raise GenomeValidationError(error_msg)

                # Check sequence length
                if len(record.seq) == 0 and not self.settings.allow_empty_sequences:
                    error_msg = f"First sequence '{record.id}' has zero length"
                    self.logger.add_validation_issue(
                        level='ERROR',
                        category='genome',
                        message=error_msg,
                        details={'sequence_id': record.id, 'index': 0}
                    )
                    raise GenomeValidationError(error_msg)

                self.logger.debug(f"✓ First sequence validated: {record.id} ({len(record.seq)} bp)")
            return

        # Strict mode - full validation
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

        Behavior depends on validation_level:
        - 'strict': Apply all edits
        - 'trust': Apply all edits (all sequences loaded)
        - 'minimal': Skip edits (file will be copied as-is)
        """
        self.logger.debug("Applying editing specifications from settings...")

        # Minimal mode - skip edits, file will be copied as-is
        if self.settings.validation_level == 'minimal':
            self.logger.debug("Minimal mode - skipping edits")
            return

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

            # If all sequences were filtered out, return early
            if len(self.sequences) == 0:
                self.logger.debug("All sequences filtered out")
                return

        # 2. Add sequence prefix
        #   TODO: before or after plasmid handling? - what is plasmid file correspond to ref genome, where the id is changes??
        #   TODO: plasmidy nepotrebuji stejny ID
        #   TODO: _ref_plasmid mergnout s _ref_plasmid ze zvlastniho souboru.
        if self.settings.replace_id_with:
            prefix = self.settings.replace_id_with
            for record in self.sequences:
                # Add prefix to sequence ID
                record.description = f"{record.id}"
                record.id = f"{prefix}"
            self.logger.debug(f"Added prefix '{prefix}' to all sequence IDs")

        # 3. Handle plasmid sequences
        # select main sequence
        if self.settings.is_plasmid:
            # Treat all sequences as plasmids (no main chromosome)
            self._handle_plasmids(self.sequences)
            self.sequences = []
        elif self.settings.plasmid_split or self.settings.plasmids_to_one:
            # Only select main sequence if we're actually doing plasmid handling
            self.main_sequence,plasmid_sequences = self._select_main_sequence(self.sequences)
            self._handle_plasmids(plasmid_sequences)
            self.sequences = [self.main_sequence]
        # else: keep all sequences in main output file (default behavior)

        self.logger.debug(f"✓ Edits applied, {len(self.sequences)} sequence(s) remaining")

    def _select_main_sequence(self, sequences: List[SeqRecord]) -> tuple[SeqRecord, List[SeqRecord]]:
        """
        Select main chromosome from sequences based on settings.

        Args:
            sequences: List of SeqRecord objects to choose from

        Returns:
            Tuple of (main_sequence, plasmid_sequences)

        Raises:
            ValueError: If neither main_longest nor main_first is True
        """
        if self.settings.main_longest:
            # Sort by length (longest first)
            sorted_sequences = sorted(sequences, key=lambda x: len(x.seq), reverse=True)
            main_sequence = sorted_sequences[0]
            plasmid_sequences = sorted_sequences[1:]
            self.logger.debug(f"Selected longest sequence as main: {main_sequence.id} ({len(main_sequence.seq)} bp)")
        elif self.settings.main_first:
            # Keep original order, first is main
            main_sequence = sequences[0]
            plasmid_sequences = sequences[1:]
            self.logger.debug(f"Selected first sequence as main: {main_sequence.id} ({len(main_sequence.seq)} bp)")
        else:
            # This should not happen due to __post_init__ validation
            raise ValueError("Either main_longest or main_first must be True")

        return main_sequence, plasmid_sequences

    def _handle_plasmids(self, plasmid_sequences:List[SeqRecord]):
        if len(plasmid_sequences) == 0:
            return

        self.logger.debug(f"Handling plasmid file with {len(plasmid_sequences)} sequence(s)...")

        if self.settings.plasmid_split:
            # Save each plasmid to separate file
            self.logger.info(
                f"Processing plasmid file: saving {len(plasmid_sequences)} plasmid(s) "
                f"to separate files"
            )
            for i, plasmid in enumerate(plasmid_sequences):
                # Save each plasmid individually (pass as list to reuse existing method)
                self._save_plasmid_file([plasmid], i)

        elif self.settings.plasmids_to_one:
            # Save all plasmids to one merged file
            self.logger.info(
                f"Processing plasmid file: merging {len(plasmid_sequences)} plasmid(s) "
                f"into one file"
            )
            self._save_plasmid_file(plasmid_sequences,"")

        else:
            # Keep all sequences in main output file
            self.logger.info(
                f"Processing plasmid file: keeping all {len(plasmid_sequences)} "
                f"sequence(s) in main output file"
            )

    def _save_plasmid_file(self, plasmid_sequences: List[SeqRecord], index: str) -> None:
        """
        Save plasmid sequences to a separate file.

        Args:
            plasmid_sequences: List of SeqRecord objects for plasmids (usually one plasmid)
            index: Index number to append to filename
        """
        # Determine output directory (with optional subdirectory)
        output_dir = self.output_dir
        if self.settings.output_subdir_name and self.settings.output_subdir_name != "plasmid":
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

        Behavior depends on validation_level:
        - 'strict': Write sequences using BioPython (may convert to FASTA)
        - 'trust': Write sequences using BioPython (all sequences with edits applied)
        - 'minimal': Copy file as-is (preserve original format)

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

        # Minimal mode - copy file as-is without parsing
        # Required: FASTA format + NO compression
        if self.settings.validation_level == 'minimal':
            self.logger.debug("Minimal mode - validating format and coding requirements")

            # Check format - must be FASTA
            if self.genome_config.detected_format != GenomeFormat.FASTA:
                error_msg = f'Minimal mode requires FASTA format, got {self.genome_config.detected_format}. Use validation_level "trust" or "strict" to convert.'
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='genome',
                    message=error_msg,
                    details={'file': self.genome_config.filename, 'detected_format': str(self.genome_config.detected_format)}
                )
                raise GenomeValidationError(error_msg)

            # Check coding - must be uncompressed (NONE)
            if self.genome_config.coding_type and self.genome_config.coding_type != CT.NONE:
                error_msg = f'Minimal mode requires uncompressed FASTA, got {self.genome_config.coding_type}. Use validation_level "trust" or "strict" to change compression.'
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='genome',
                    message=error_msg,
                    details={'file': self.genome_config.filename, 'coding_type': str(self.genome_config.coding_type)}
                )
                raise GenomeValidationError(error_msg)

            # Normalize output filename - always .fasta extension
            if self.settings.output_filename_suffix:
                output_filename_minimal = f"{base_name}_{self.settings.output_filename_suffix}.fasta"
            else:
                output_filename_minimal = f"{base_name}.fasta"

            output_path_minimal = output_dir / output_filename_minimal

            # Copy file
            self.logger.debug(f"Copying {self.input_path} to {output_path_minimal}")
            shutil.copy2(self.input_path, output_path_minimal)

            self.logger.info(f"Output saved: {output_path_minimal}")
            return output_path_minimal

        # Strict and Trust modes - write sequences using BioPython
        # Check if we have sequences to write
        if self.sequences == []:
            self.logger.warning("No sequences to write")
            return None

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

        return
