"""
Read file validator and processor.

Handles FASTQ and BAM formats with compression support.
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
from validation_pkg.exceptions import (
    ReadValidationError,
    FileFormatError,
    FastqFormatError,
    BamFormatError,
    CompressionError
)
from validation_pkg.utils.formats import ReadFormat
from validation_pkg.utils.formats import CodingType as CT
from validation_pkg.utils.file_handler import (
    bz2_to_gz,
    none_to_gz,
    gz_to_bz2,
    gz_to_none,
    bz2_to_none,
    none_to_bz2
)

class ReadValidator:
    """
    Validates and processes read files (FASTQ and BAM formats).

    Workflow:
    1. Detect and handle compression
    2. Detect file format
    3. Parse and validate using BioPython
    4. Apply editing specifications
    5. Convert to FASTQ format
    6. Compress if requested
    7. Save to output directory
    """

    @dataclass
    class Settings(BaseSettings):
        """
        Settings for read validation and processing.

        Attributes:
            validation_level: Validation thoroughness level:
                - 'strict': Full decompression, parsing, and validation of all sequences (default)
                - 'trust': Fast validation - check line count and first record only
                - 'minimal': Only verify file exists and has correct extension
            check_invalid_chars: If True, validates that sequences contain only valid nucleotides (ATCGN)
            allow_empty_id: If True, allows sequences without IDs
            allow_duplicate_ids: If True, allows duplicate sequence IDs across the file
            keep_bam: If True, keeps a copy of the original BAM file when converting to FASTQ
            coding_type: Output compression type ('gz', 'bz2', or None for uncompressed)
            output_filename_suffix: Optional suffix to add to output filename (e.g., 'filtered')
            output_subdir_name: Optional subdirectory name within output_dir for saving files
            outdir_by_ngs_type: If True, automatically set output_subdir_name to ngs_type (default: False)

        Note:
            When outdir_by_ngs_type=True, output_subdir_name will be automatically set to the
            read's ngs_type (illumina/ont/pacbio), overriding any manually set value.

        Trust Mode Behavior:
            - validation_level='strict': Full validation (current behavior)
            - validation_level='trust': FASTQ - check line count % 4 == 0 and first record format
                                        BAM - verify header and first record only
            - validation_level='minimal': Only file existence and extension check

        Example:
            >>> settings = ReadValidator.Settings()
            >>> print(settings)  # View all settings
            >>> ref_settings = settings.update(sequence_prefix="chr1", output_filename_suffix="ref")
            >>> validator = ReadValidator(read_config, output_dir, ref_settings)
        """
        # Validation level
        validation_level: str = 'strict'  # 'strict', 'trust', or 'minimal'

        # Validation options
        check_invalid_chars: bool = False 
        allow_empty_id: bool = False
        allow_duplicate_ids: bool = True

        # Editing specifications
        keep_bam: bool = True 
        ignore_bam: bool = True

        # Output format
        coding_type: Optional[str] = None
        output_filename_suffix: Optional[str] = None
        output_subdir_name: Optional[str] = None
        outdir_by_ngs_type: bool = False

    def __init__(self, read_config, output_dir: Path, settings: Optional[Settings] = None) -> None:
        """
        Initialize read validator.

        Args:
            read_config: ReadConfig object from ConfigManager with file info
            output_dir: Directory for output files (Path object)
            settings: Settings object with validation parameters (uses defaults if None)

        Note: In future versions, output_dir may be moved to settings.

        Example:
            >>> settings = ReadValidator.Settings()
            >>> settings = settings.update(sequence_prefix="chr1")
            >>> validator = ReadValidator(read_config, output_dir, settings)
        """
        self.logger = get_logger()
        self.read_config = read_config
        self.output_dir = Path(output_dir)
        self.settings = settings if settings is not None else self.Settings()

        # Validate validation_level setting
        valid_levels = {'strict', 'trust', 'minimal'}
        if self.settings.validation_level not in valid_levels:
            raise ValueError(
                f"Invalid validation_level '{self.settings.validation_level}'. "
                f"Must be one of: {', '.join(valid_levels)}"
            )

        # Apply outdir_by_ngs_type if enabled
        if self.settings.outdir_by_ngs_type:
            # Override output_subdir_name with ngs_type
            self.settings = self.settings.update(output_subdir_name=read_config.ngs_type)
            self.logger.debug(f"Applied outdir_by_ngs_type: output_subdir_name set to '{read_config.ngs_type}'")

        # Log settings being used
        self.logger.debug(f"Initializing ReadValidator with settings:\n{self.settings}")

        # Resolved paths
        self.input_path = read_config.filepath

        # Parsed data
        self.sequences = []  # List of SeqRecord objects
    
    def validate(self) -> None:
        """
        Main validation and processing workflow.

        Uses read_config data (format, compression) provided by ConfigManager.

        Workflow differs based on input format:
        - FASTQ: Parse → Validate → Edit → Statistics → Save
        - BAM: Copy original → Convert to FASTQ → Parse → Validate → Edit → Statistics → Save

        Raises:
            ReadValidationError: If validation fails
        """
        self.logger.info(f"Processing read file: {self.read_config.filename}")
        self.logger.debug(f"Format: {self.read_config.detected_format}, Compression: {self.read_config.coding_type}")

        try:
            # Special workflow for BAM files
            if self.read_config.detected_format == ReadFormat.BAM:
                if self.settings.ignore_bam:
                    self.logger.warning(f"Cannot proccess BAM files, the file will be ignored")
                    return

                # Step 1: Copy original BAM to output (if keep_bam is enabled)
                if self.settings.keep_bam:
                    self._copy_bam_to_output()

                # Step 2: Convert BAM to FASTQ (populates self.sequences)
                self._convert_bam_to_fastq()

                # Step 3: Validate sequences
                self._validate_sequences()

                # Step 4: Apply editing specifications
                self._apply_edits()

                # Step 6: Save FASTQ output
                output_path = self._save_output()

            else:
                # Standard FASTQ workflow
                # Step 1: Parse and validate
                self._parse_file()

                # Validate sequences
                self._validate_sequences()
                
                # Step 2: Apply editing specifications
                self._apply_edits()

                # Step 4: Save to output directory (already FASTQ)
                self._save_output()

            self.logger.info(f"✓ Read validation completed")

        except Exception as e:
            self.logger.error(f"Read validation failed: {e}")
            raise

    def _open_file(self, mode: str = 'rt') -> IO:
        """
        Open file with automatic decompression based on read_config.

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
                self.read_config.coding_type,
                mode
            )
        except CompressionError as e:
            self.logger.add_validation_issue(
                level='ERROR',
                category='read',
                message=str(e),
                details={'file': str(self.input_path)}
            )
            raise

    def _count_lines_fast(self) -> int:
        """
        Fast line counting using zcat/bzcat for compressed files.

        Uses shell commands for efficient line counting without loading file into memory:
        - gzip (.gz): zcat file.gz | wc -l
        - bzip2 (.bz2): bzcat file.bz2 | wc -l
        - uncompressed: wc -l file

        Note: wc -l counts newline characters, not lines. If a file doesn't end with
        a newline, we need to add 1 to the count. We check if the last byte is a newline.

        Returns:
            Number of lines in the file

        Raises:
            ReadValidationError: If line counting fails
        """
        import subprocess

        try:
            # Count newline characters
            if self.read_config.coding_type == CT.GZIP:
                # Use zcat for gzip files
                result = subprocess.run(
                    ['sh', '-c', f'zcat "{self.input_path}" | wc -l'],
                    capture_output=True,
                    text=True,
                    timeout=300  # 5 minutes timeout
                )
            elif self.read_config.coding_type == CT.BZIP2:
                # Use bzcat for bzip2 files
                result = subprocess.run(
                    ['sh', '-c', f'bzcat "{self.input_path}" | wc -l'],
                    capture_output=True,
                    text=True,
                    timeout=300  # 5 minutes timeout
                )
            else:
                # Use wc -l for uncompressed files
                result = subprocess.run(
                    ['wc', '-l', str(self.input_path)],
                    capture_output=True,
                    text=True,
                    timeout=60
                )

            if result.returncode != 0:
                raise ReadValidationError(f"Line counting failed: {result.stderr}")

            # Parse output - wc -l returns "count filename" or just "count"
            line_count = int(result.stdout.strip().split()[0])

            # Check if file ends with newline - wc -l counts newline chars, not lines
            # If file doesn't end with newline, we need to add 1
            if self.read_config.coding_type == CT.GZIP:
                # Check last byte of decompressed content
                last_byte_result = subprocess.run(
                    ['sh', '-c', f'zcat "{self.input_path}" | tail -c 1'],
                    capture_output=True,
                    timeout=30
                )
                if last_byte_result.returncode == 0 and last_byte_result.stdout and last_byte_result.stdout[-1:] != b'\n':
                    line_count += 1
            elif self.read_config.coding_type == CT.BZIP2:
                # Check last byte of decompressed content
                last_byte_result = subprocess.run(
                    ['sh', '-c', f'bzcat "{self.input_path}" | tail -c 1'],
                    capture_output=True,
                    timeout=30
                )
                if last_byte_result.returncode == 0 and last_byte_result.stdout and last_byte_result.stdout[-1:] != b'\n':
                    line_count += 1
            else:
                # Check last byte of uncompressed file
                if self.input_path.stat().st_size > 0:
                    with open(self.input_path, 'rb') as f:
                        f.seek(-1, 2)  # Seek to last byte
                        last_byte = f.read(1)
                        if last_byte != b'\n':
                            line_count += 1

            return line_count

        except subprocess.TimeoutExpired:
            raise ReadValidationError("Line counting timed out - file too large")
        except (ValueError, IndexError) as e:
            raise ReadValidationError(f"Failed to parse line count: {e}")
        except Exception as e:
            raise ReadValidationError(f"Line counting failed: {e}")

    def _validate_first_fastq_record(self) -> Optional[SeqRecord]:
        """
        Validate only the first FASTQ record for format correctness.

        Checks:
        1. Line 1 starts with '@' (sequence ID)
        2. Line 2 contains sequence
        3. Line 3 starts with '+' (separator)
        4. Line 4 contains quality scores

        Returns:
            First SeqRecord if valid, None if validation fails

        Raises:
            FastqFormatError: If first record is not in valid FASTQ format
        """
        self.logger.debug("Validating first FASTQ record...")

        try:
            handle = self._open_file()
            try:
                # Read first 4 lines
                lines = []
                for i in range(4):
                    line = handle.readline()
                    if not line:
                        raise FastqFormatError(
                            f"File has fewer than 4 lines - not a valid FASTQ file"
                        )
                    lines.append(line.rstrip('\n\r'))

                # Validate FASTQ format
                if not lines[0].startswith('@'):
                    raise FastqFormatError(
                        f"First line must start with '@', found: {lines[0][:50]}"
                    )

                if not lines[2].startswith('+'):
                    raise FastqFormatError(
                        f"Third line must start with '+', found: {lines[2][:50]}"
                    )

                # Check that sequence and quality have same length
                seq_len = len(lines[1])
                qual_len = len(lines[3])
                if seq_len != qual_len:
                    raise FastqFormatError(
                        f"Sequence length ({seq_len}) doesn't match quality length ({qual_len})"
                    )

                # Parse first record using BioPython to ensure it's valid
                from io import StringIO
                first_record_text = '\n'.join(lines)
                record = next(SeqIO.parse(StringIO(first_record_text), 'fastq'))

                self.logger.debug(f"✓ First record validated: {record.id}, length={len(record.seq)}")
                return record

            finally:
                handle.close()

        except FastqFormatError:
            raise
        except Exception as e:
            error_msg = f"Failed to validate first FASTQ record: {e}"
            self.logger.add_validation_issue(
                level='ERROR',
                category='read',
                message=error_msg,
                details={'file': self.read_config.filename, 'error': str(e)}
            )
            raise FastqFormatError(error_msg) from e
    
    def _parse_file(self) -> None:
        """
        Parse FASTQ file using BioPython and validate format from read_config.

        Behavior depends on validation_level:
        - 'strict': Full parsing and validation of all sequences
        - 'trust': Only validate first record and check line count
        - 'minimal': No parsing or validation

        Note: BAM files are handled separately in the validate() method.
        This method only processes FASTQ files.
        """
        self.logger.debug(f"Parsing {self.read_config.detected_format} file (validation_level={self.settings.validation_level})...")

        # Minimal mode - skip all parsing and validation except CodingType
        if self.settings.validation_level == 'minimal':
            self.logger.info("Minimal validation mode - skipping file parsing")
            # Set empty sequences list - file will be copied as-is in _save_output
            self.sequences = []
            return

        # Trust mode: Parse only first 10 sequences for validation
        # Strict mode: Parse all sequences
        if self.settings.validation_level == 'trust':
            self.logger.info("Trust mode - parsing first 10 sequences only for validation")
            parse_limit = 10
        else:
            # Strict mode - parse all
            parse_limit = None

        try:
            # Get total line count for progress reporting (strict mode only)
            if self.settings.validation_level == 'strict':
                try:
                    line_count = self._count_lines_fast()
                    estimated_sequences = line_count // 4  # FASTQ has 4 lines per sequence
                    self.logger.info(f"Processing {estimated_sequences:,} reads from FASTQ file...")
                    progress_interval = max(1, estimated_sequences // 20)  # Report every 5%
                except Exception as e:
                    # If counting fails, just proceed without progress reporting
                    self.logger.debug(f"Could not count lines for progress reporting: {e}")
                    estimated_sequences = None
                    progress_interval = 100000  # Report every 100k reads
            else:
                estimated_sequences = None
                progress_interval = None

            # Open file with automatic decompression
            handle = self._open_file()
            try:
                # Parse using BioPython with clean enum conversion
                biopython_format = self.read_config.detected_format.to_biopython()

                # Parse sequences
                self.sequences = []
                processed = 0
                for record in SeqIO.parse(handle, biopython_format):
                    self.sequences.append(record)
                    processed += 1

                    # Trust mode - stop after parsing limit
                    if parse_limit and processed >= parse_limit:
                        break

                    # Progress reporting (strict mode only)
                    if progress_interval and processed % progress_interval == 0:
                        if estimated_sequences:
                            percent = (processed / estimated_sequences) * 100
                            self.logger.info(f"Progress: {processed:,}/{estimated_sequences:,} reads ({percent:.1f}%)")
                        else:
                            self.logger.info(f"Progress: {processed:,} reads processed...")
            finally:
                handle.close()

            # Validate we got sequences
            if not self.sequences:
                error_msg = f"No sequences found in {self.read_config.detected_format} file"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='read',
                    message=error_msg,
                    details={'file': self.read_config.filename, 'format': self.read_config.detected_format}
                )
                raise ReadValidationError(error_msg)

            self.logger.info(f"✓ Parsed {len(self.sequences):,} sequence(s)")

        except FileFormatError:
            raise
        except Exception as e:
            error_msg = f"Failed to parse {self.read_config.detected_format} file: {e}"

            if self.read_config.detected_format == ReadFormat.FASTQ:
                exception_class = FastqFormatError
            elif self.read_config.detected_format == ReadFormat.BAM :
                exception_class = BamFormatError
            else:
                exception_class = FileFormatError

            self.logger.add_validation_issue(
                level='ERROR',
                category='read',
                message=error_msg,
                details={
                    'file': self.read_config.filename,
                    'format': self.read_config.detected_format,
                    'error': str(e)
                }
            )
            raise exception_class(error_msg) from e
    
    def _validate_sequences(self) -> None:
        """
        Validate parsed sequences.

        Behavior depends on validation_level:
        - 'strict': Validate all sequences
        - 'trust': Validate only first 10 sequences
        - 'minimal': Skip (no sequences parsed)
        """
        # Determine how many sequences to validate
        if self.settings.validation_level == 'trust':
            validate_count = min(10, len(self.sequences))
            self.logger.info(f"Trust mode - validating first {validate_count} of {len(self.sequences):,} sequences")
        else:
            validate_count = len(self.sequences)
            self.logger.debug(f"Validating {validate_count:,} sequences...")

        # Validate sequences
        for idx in range(validate_count):
            record = self.sequences[idx]
            # Check sequence ID
            if not record.id and not self.settings.allow_empty_id:
                error_msg = f"Sequence at index {idx} has no ID"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='read',
                    message=error_msg,
                    details={'sequence_index': idx, 'sequence_id': str(record.id)}
                )
                raise ReadValidationError(error_msg)
            
            
            # Check for valid nucleotides
            if self.settings.check_invalid_chars:
                valid_chars = set('ATCGNatcgn')
                seq_str = str(record.seq)
                invalid_chars = set(seq_str) - valid_chars

                if invalid_chars:
                    char_count = sum(1 for c in seq_str if c in invalid_chars)
                    error_msg = f"Sequence '{record.id}' contains {char_count} invalid character(s): {', '.join(sorted(invalid_chars))}"
                    self.logger.add_validation_issue(
                        level='ERROR',
                        category='read',
                        message=error_msg,
                        details={
                            'sequence_id': record.id,
                            'invalid_chars': list(invalid_chars),
                            'count': char_count
                        }
                    )
                    raise ReadValidationError(error_msg)

        # 3. Check for duplicate IDs
        if not self.settings.allow_duplicate_ids:
            seq_ids = [record.id for record in self.sequences]
            if len(seq_ids) != len(set(seq_ids)):
                duplicates = [sid for sid in set(seq_ids) if seq_ids.count(sid) > 1]
                error_msg = f"Duplicate sequence IDs not allowed: {duplicates}"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='read',
                    message=error_msg,
                    details={'duplicate_ids': duplicates}
                )
                raise ReadValidationError(error_msg)
        
        self.logger.debug("✓ Sequence validation passed")
    
    def _apply_edits(self) -> None:
        """
        Apply editing specifications to sequences based on settings.

        Future edits may include:
        - Remove short sequences (min_read_length)
        - Quality filtering
        - Adapter trimming
        - Read ID prefix/suffix modification
        """
        self.logger.debug("Applying editing specifications from settings...")

        # Currently no edits are implemented
        # This is a placeholder for future functionality

        self.logger.debug(f"✓ Edits applied, {len(self.sequences)} sequence(s) remaining")
    
    def _convert_bam_to_fastq(self) -> None:
        """
        Convert BAM file to FASTQ format using pysam or samtools.

        Behavior depends on validation_level:
        - 'strict': Full conversion of all reads
        - 'trust': Only validate BAM header and first few reads
        - 'minimal': Skip conversion entirely

        This method tries to use pysam first, and falls back to samtools
        command-line tool if pysam is not available.

        Raises:
            ReadValidationError: If conversion fails or neither tool is available
        """
        import subprocess

        # Minimal mode - skip conversion
        if self.settings.validation_level == 'minimal':
            self.logger.info("Minimal validation mode - skipping BAM conversion")
            self.sequences = []
            return

        # Trust mode - validate header and first few reads only
        if self.settings.validation_level == 'trust':
            self.logger.info("Trust mode - validating BAM header and first records only")
            try:
                import pysam  # type: ignore

                with pysam.AlignmentFile(str(self.input_path), "rb") as bam_file:
                    # Validate BAM header exists
                    if not bam_file.header:
                        raise BamFormatError("BAM file has no header")

                    self.logger.debug(f"✓ BAM header validated: {len(bam_file.header.get('SQ', []))} sequences in reference")

                    # Read and validate first few records (up to 10)
                    first_records = []
                    for i, read in enumerate(bam_file):
                        if i >= 10:
                            break

                        # Skip unmapped or secondary/supplementary alignments
                        if read.is_unmapped or read.is_secondary or read.is_supplementary:
                            continue

                        # Create SeqRecord from BAM read
                        seq = Seq(read.query_sequence if read.query_sequence else "")
                        qualities = read.query_qualities if read.query_qualities else []

                        record = SeqRecord(
                            seq,
                            id=read.query_name,
                            description="",
                            letter_annotations={"phred_quality": list(qualities)} if qualities else {}
                        )
                        first_records.append(record)

                    if not first_records:
                        raise BamFormatError("No valid reads found in BAM file")

                    self.logger.info(f"✓ Trust mode BAM validation passed - validated {len(first_records)} records")
                    self.sequences = first_records
                    return

            except ImportError:
                # If pysam not available, try samtools for header validation
                try:
                    result = subprocess.run(
                        ["samtools", "view", "-H", str(self.input_path)],
                        capture_output=True,
                        text=True,
                        timeout=30
                    )

                    if result.returncode != 0:
                        raise BamFormatError(f"Failed to read BAM header: {result.stderr}")

                    if not result.stdout.strip():
                        raise BamFormatError("BAM file has no header")

                    self.logger.info("✓ Trust mode BAM validation passed (header only, pysam not available)")
                    self.sequences = []
                    return

                except (FileNotFoundError, subprocess.TimeoutExpired) as e:
                    raise ReadValidationError(
                        "Trust mode BAM validation requires either 'pysam' or 'samtools'. "
                        "Please install one of them."
                    ) from e

        # Strict mode - full conversion (original behavior)
        # Try using pysam first (optional dependency)
        try:
            import pysam  # type: ignore

            self.logger.debug("Using pysam for BAM to FASTQ conversion")
            self.sequences = []

            with pysam.AlignmentFile(str(self.input_path), "rb") as bam_file:
                # Try to get total read count for progress reporting
                try:
                    total_reads = bam_file.count(until_eof=True)
                    bam_file.reset()  # Reset file pointer after counting
                    self.logger.info(f"Processing {total_reads:,} reads from BAM file...")
                    progress_interval = max(1, total_reads // 20)  # Report every 5%
                except:
                    # If counting fails, just proceed without progress reporting
                    total_reads = None
                    progress_interval = 100000  # Report every 100k reads

                processed = 0
                for read in bam_file:
                    # Skip unmapped or secondary/supplementary alignments
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue

                    # Create SeqRecord from BAM read
                    seq = Seq(read.query_sequence if read.query_sequence else "")
                    qualities = read.query_qualities if read.query_qualities else []

                    record = SeqRecord(
                        seq,
                        id=read.query_name,
                        description="",
                        letter_annotations={"phred_quality": list(qualities)} if qualities else {}
                    )
                    self.sequences.append(record)

                    processed += 1
                    # Progress reporting
                    if processed % progress_interval == 0:
                        if total_reads:
                            percent = (processed / total_reads) * 100
                            self.logger.info(f"Progress: {processed:,}/{total_reads:,} reads ({percent:.1f}%)")
                        else:
                            self.logger.info(f"Progress: {processed:,} reads processed...")

            self.logger.info(f"✓ Converted {len(self.sequences):,} reads from BAM to FASTQ")
            return

        except ImportError:
            self.logger.debug("pysam not available, trying samtools...")

        # Fall back to samtools
        try:
            # Check if samtools is available
            result = subprocess.run(
                ["samtools", "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )

            if result.returncode != 0:
                raise FileNotFoundError("samtools not found")

            self.logger.debug("Using samtools for BAM to FASTQ conversion")

            # Use samtools fastq to convert BAM to FASTQ
            # samtools fastq writes to stdout by default
            result = subprocess.run(
                ["samtools", "fastq", str(self.input_path)],
                capture_output=True,
                text=True,
                timeout=300  # 5 minutes timeout
            )

            if result.returncode != 0:
                error_msg = f"samtools conversion failed: {result.stderr}"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='read',
                    message=error_msg,
                    details={'file': self.read_config.filename, 'error': result.stderr}
                )
                raise ReadValidationError(error_msg)

            # Parse FASTQ output from samtools
            from io import StringIO
            fastq_data = StringIO(result.stdout)
            self.sequences = list(SeqIO.parse(fastq_data, "fastq"))

            self.logger.info(f"Converted {len(self.sequences)} reads from BAM to FASTQ using samtools")

        except (FileNotFoundError, subprocess.TimeoutExpired, subprocess.SubprocessError) as e:
            error_msg = (
                "BAM to FASTQ conversion requires either 'pysam' Python package or 'samtools' command-line tool. "
                "Please install one of them:\n"
                "  - pip install pysam\n"
                "  - or install samtools (https://www.htslib.org/)"
            )
            self.logger.add_validation_issue(
                level='ERROR',
                category='read',
                message=error_msg,
                details={'error': str(e)}
            )
            raise ReadValidationError(error_msg) from e
    
    def _copy_bam_to_output(self) -> Path:
        """
        Copy BAM file to output directory without processing.

        Returns:
            Path to output file
        """
        self.logger.debug("Copying BAM file to output directory...")

        # Determine output directory (with optional subdirectory)
        output_dir = self.output_dir
        if self.settings.output_subdir_name:
            output_dir = output_dir / self.settings.output_subdir_name

        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate output filename - keep original name and compression
        output_filename = self.read_config.filename
        if self.settings.output_filename_suffix:
            # Insert suffix before extensions
            base_name = self.read_config.filename
            for suffix in self.input_path.suffixes:
                base_name = base_name.replace(suffix, '')
            # Reconstruct with suffix and original extensions
            extensions = ''.join(self.input_path.suffixes)
            output_filename = f"{base_name}_{self.settings.output_filename_suffix}{extensions}"

        output_path = output_dir / output_filename

        # Copy file
        self.logger.debug(f"Copying {self.input_path} to {output_path}")
        shutil.copy2(self.input_path, output_path)

        self.logger.info(f"BAM file copied: {output_path}")
        return output_path

    def _save_output(self) -> Path:
        """
        Save processed read to output directory using settings.

        Behavior depends on validation_level:
        - 'strict': Write sequences using BioPython
        - 'trust' or 'minimal': Copy file as-is (preserve original format)

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
        base_name = self.read_config.filename

        # Remove all suffixes (.fastq.gz -> remove .gz then .fastq)
        for suffix in self.input_path.suffixes:
            base_name = base_name.replace(suffix, '')
        # Add suffix if specified
        if self.settings.output_filename_suffix:
            output_filename = f"{base_name}_{self.settings.output_filename_suffix}.fastq"
        else:
            output_filename = f"{base_name}.fastq"


        # Add compression extension if requested (from settings)
        # Normalize settings.coding_type to CodingType enum
        coding = CT(self.settings.coding_type) if self.settings.coding_type else CT.NONE

        if coding == CT.GZIP:
            output_filename += '.gz'
        elif coding == CT.BZIP2:
            output_filename += '.bz2'

        output_path = output_dir / output_filename

        # Minimal mode - copy file as-is without parsing
        # Required: FASTQ format + coding must match settings.coding_type
        if self.settings.validation_level == 'minimal':
            self.logger.debug("Minimal mode - validating format and coding requirements")

            # Check format - must be FASTQ
            if self.read_config.detected_format != ReadFormat.FASTQ:
                error_msg = f'Minimal mode requires FASTQ format, got {self.read_config.detected_format}. Use validation_level "trust" or "strict" to convert.'
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='read',
                    message=error_msg,
                    details={'file': self.read_config.filename, 'detected_format': str(self.read_config.detected_format)}
                )
                raise ReadValidationError(error_msg)

            # Check coding - must match settings.coding_type
            # read_config.coding_type is CodingType enum, coding is normalized above
            input_coding = self.read_config.coding_type
            required_coding = coding

            if input_coding != required_coding:
                error_msg = f'Minimal mode requires input coding to match output coding. Input: {input_coding}, Required: {required_coding}. Use validation_level "trust" or "strict" to change compression.'
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='read',
                    message=error_msg,
                    details={'file': self.read_config.filename, 'input_coding': str(input_coding), 'required_coding': str(required_coding)}
                )
                raise ReadValidationError(error_msg)

            # Check line count divisible by 4 (fast validation)
            try:
                line_count = self._count_lines_fast()
                if line_count % 4 != 0:
                    error_msg = f"Invalid FASTQ: {line_count:,} lines not divisible by 4"
                    self.logger.add_validation_issue(
                        level='ERROR',
                        category='read',
                        message=error_msg,
                        details={'file': self.read_config.filename, 'line_count': line_count}
                    )
                    raise FastqFormatError(error_msg)
                num_sequences = line_count // 4
                self.logger.info(f"Minimal mode - verified {num_sequences:,} sequences (line count check)")
            except FastqFormatError:
                raise
            except Exception as e:
                error_msg = f"Minimal mode validation failed: {e}"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='read',
                    message=error_msg,
                    details={'file': self.read_config.filename, 'error': str(e)}
                )
                raise ReadValidationError(error_msg) from e

            # Use the output_filename already constructed (includes correct extension)
            output_path_minimal = output_dir / output_filename

            # Copy file
            self.logger.debug(f"Copying {self.input_path} to {output_path_minimal}")
            shutil.copy2(self.input_path, output_path_minimal)

            self.logger.info(f"Output saved: {output_path_minimal}")
            return output_path_minimal

        # Trust mode - copy original file with coding conversion
        # Trust mode parsed only first 10 sequences for validation, so we copy the original file
        if self.settings.validation_level == 'trust':
            self.logger.debug("Trust mode - copying original file with coding conversion")

            # Both are CodingType enum: read_config.coding_type and coding (normalized above)
            input_coding = self.read_config.coding_type
            output_coding = coding

            # If input and output coding match, simple copy
            if input_coding == output_coding:
                self.logger.debug(f"Copying {self.input_path} to {output_path} (same coding)")
                shutil.copy2(self.input_path, output_path)
            else:
                # Convert coding using file_handler utilities
                self.logger.debug(f"Converting {input_coding} -> {output_coding}")

                # Map (input, output) pairs to conversion functions
                conversion_map = {
                    (CT.BZIP2, CT.GZIP): bz2_to_gz,
                    (CT.NONE, CT.GZIP): none_to_gz,
                    (CT.GZIP, CT.NONE): gz_to_none,
                    (CT.BZIP2, CT.NONE): bz2_to_none,
                    (CT.NONE, CT.BZIP2): none_to_bz2,
                    (CT.GZIP, CT.BZIP2): gz_to_bz2,
                }

                conversion_func = conversion_map.get((input_coding, output_coding))
                if conversion_func:
                    conversion_func(self.input_path, output_path)
                else:
                    error_msg = f"Unsupported coding conversion: {input_coding} -> {output_coding}"
                    raise ReadValidationError(error_msg)

            self.logger.info(f"Output saved: {output_path}")
            return output_path

        # Strict mode - write sequences using BioPython (original behavior)
        # Write output with appropriate compression
        self.logger.debug(f"Writing output to: {output_path}")

        if coding == CT.GZIP:
            with gzip.open(output_path, 'wt') as handle:
                SeqIO.write(self.sequences, handle, 'fastq')
        elif coding == CT.BZIP2:
            with bz2.open(output_path, 'wt') as handle:
                SeqIO.write(self.sequences, handle, 'fastq')
        else:
            with open(output_path, 'w') as handle:
                SeqIO.write(self.sequences, handle, 'fastq')

        self.logger.info(f"Output saved: {output_path}")

        return output_path
    