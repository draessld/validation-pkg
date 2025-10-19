"""
Feature file validator and processor.

Handles GFF, GTF, and BED formats with compression support.
"""

from pathlib import Path
from typing import Optional, Dict, Any, List, IO
from dataclasses import dataclass
import gzip
import bz2

from validation_pkg.logger import get_logger
from validation_pkg.utils.settings import BaseSettings
from validation_pkg.exceptions import (
    FeatureValidationError,
    FileFormatError,
    BedFormatError,
    GffFormatError,
    CompressionError,
)
from validation_pkg.utils.formats import FeatureFormat

@dataclass
class Feature:
    """Represents a genomic feature (generic container for GFF/BED data)."""
    seqname: str
    start: int
    end: int
    feature_type: str = ""
    score: str = "."
    strand: str = "+"
    source: str = "."
    frame: str = "."
    attributes: str = ""
    original_format: str = "gff"

    @property
    def length(self) -> int:
        """Calculate feature length."""
        return self.end - self.start


class FeatureValidator:
    """
    Validates and processes feature annotation files (GFF, GTF, BED formats).

    Workflow:
    1. Detect and handle compression
    2. Parse features from file
    3. Validate feature coordinates and structure
    4. Apply editing specifications (sort, replace IDs)
    5. Collect statistics
    6. Convert to GFF format (if needed)
    7. Save to output directory
    """

    @dataclass
    class Settings(BaseSettings):
        """
        Settings for feature validation and processing.

        Attributes:
            sort_by_position: Sort features by position (default: True)
            check_coordinates: Validate that start < end (default: True)
            allow_zero_length: Allow zero-length features (start == end) (default: False)
            allow_negative_coords: Allow negative coordinates (default: False)
            replace_id_with: Replace the seqname/chromosome field (column 1) for all features
                           with this value. Original seqname stored in Note attribute.
                           (e.g., 'chr1' replaces all seqnames with 'chr1') (default: None)
            coding_type: Output compression type: 'gz', 'bz2', or None (default: None)
            output_filename_suffix: Suffix to add to output filename (default: None)
            output_subdir_name: Subdirectory name for output files (default: None)

        Example:
            >>> settings = FeatureValidator.Settings()
            >>> settings = settings.update(sort_by_position=True, coding_type='gz', replace_id_with='chr1')
            >>> validator = FeatureValidator(feature_config, output_dir, settings)
        """
        # Validation options
        sort_by_position: bool = True
        check_coordinates: bool = True
        allow_zero_length: bool = False
        allow_negative_coords: bool = False

        # Editing options
        replace_id_with: Optional[str] = None

        # Output format
        coding_type: Optional[str] = None
        output_filename_suffix: Optional[str] = None
        output_subdir_name: Optional[str] = None


    def __init__(self, feature_config, output_dir: Path, settings: Optional[Settings] = None) -> None:
        """
        Initialize feature validator.

        Args:
            feature_config: FeatureConfig object from ConfigManager with file info
            output_dir: Directory for output files (Path object)
            settings: Settings object with validation parameters (uses defaults if None)

        Example:
            >>> settings = FeatureValidator.Settings()
            >>> settings = settings.update(sort_by_position=True)
            >>> validator = FeatureValidator(feature_config, output_dir, settings)
        """
        self.logger = get_logger()
        self.feature_config = feature_config
        self.output_dir = Path(output_dir)
        self.settings = settings if settings is not None else self.Settings()

        # Log settings being used
        self.logger.debug(f"Initializing FeatureValidator with settings:\n{self.settings}")

        # Resolved paths
        self.input_path = feature_config.filepath

        # Parsed data
        self.features: List[Feature] = []

        # Statistics
        self.statistics = {
            'num_features': 0,
            'num_by_type': {},
            'num_by_strand': {'+': 0, '-': 0, '.': 0},
            'total_span': 0,
            'min_length': 0,
            'max_length': 0,
            'avg_length': 0.0
        }
    
    def validate(self) -> None:
        """
        Main validation and processing workflow.

        Uses feature_config data (format, compression) provided by ConfigManager.

        Raises:
            FeatureValidationError: If validation fails
        """
        self.logger.info(f"Processing feature file: {self.feature_config.filename}")
        self.logger.debug(f"Format: {self.feature_config.detected_format}, Compression: {self.feature_config.coding_type}")

        try:
            # Step 1: Parse and validate
            self._parse_file()

            # Step 2: Apply editing specifications
            self._apply_edits()

            # Step 3: Collect statistics
            self._collect_statistics()

            # Step 4: Convert to GFF (if needed)
            self._convert_to_gff()

            # Step 5: Save to output directory
            output_path = self._save_output()

            self.logger.info(f"✓ Feature validation completed: {output_path.name}")
            self.logger.debug(f"Statistics: {self.statistics}")

        except Exception as e:
            self.logger.error(f"Feature validation failed: {e}")
            raise

    def _open_file(self, mode: str = 'rt') -> IO:
        """
        Open file with automatic decompression based on feature_config.

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
                self.feature_config.coding_type,
                mode
            )
        except CompressionError as e:
            self.logger.add_validation_issue(
                level='ERROR',
                category='feature',
                message=str(e),
                details={'file': str(self.input_path)}
            )
            raise

    def _parse_file(self) -> None:
        """Parse feature file based on detected format from feature_config."""
        format_str = str(self.feature_config.detected_format)
        self.logger.debug(f"Parsing {format_str} file...")

        try:
            with self._open_file() as handle:
                if 'gff' in format_str.lower():
                    self.features = self._parse_gff(handle)
                elif 'bed' in format_str.lower():
                    self.features = self._parse_bed(handle)
                else:
                    raise FileFormatError(f"Unsupported feature format: {format_str}")

            # Validate we got features
            if not self.features:
                error_msg = f"No features found in {format_str} file"
                self.logger.add_validation_issue(
                    level='WARNING',
                    category='feature',
                    message=error_msg,
                    details={'file': self.feature_config.filename, 'format': format_str}
                )
                # Empty feature files are allowed, just warn

            self.logger.debug(f"Parsed {len(self.features)} feature(s)")

            # Validate features
            if self.features:
                self._validate_features()

        except (FileFormatError, FeatureValidationError):
            raise
        except Exception as e:
            error_msg = f"Failed to parse {format_str} file: {e}"

            if 'gff' in format_str.lower():
                exception_class = GffFormatError
            elif 'bed' in format_str.lower():
                exception_class = BedFormatError
            else:
                exception_class = FileFormatError

            self.logger.add_validation_issue(
                level='ERROR',
                category='feature',
                message=error_msg,
                details={
                    'file': self.feature_config.filename,
                    'format': format_str,
                    'error': str(e)
                }
            )
            raise exception_class(error_msg) from e

    def _normalize_strand(self, strand: str) -> str:
        """Normalize strand value to +, -, or ."""
        return strand if strand in ['+', '-', '.'] else '.'

    def _skip_line(self, line: str) -> bool:
        """Check if line should be skipped (empty or comment)."""
        return not line or line.startswith('#')

    def _parse_gff(self, handle) -> List[Feature]:
        """
        Parse GFF/GTF format file.

        GFF format (tab-separated):
        seqname source feature start end score strand frame attributes

        Returns:
            List of Feature objects
        """
        features = []
        for line_num, line in enumerate(handle, start=1):
            line = line.strip()
            if self._skip_line(line):
                continue

            try:
                parts = line.split('\t')
                if len(parts) < 8:
                    self._log_parse_warning(
                        f"Invalid GFF line (expected 8-9 columns, got {len(parts)})",
                        line_num, line
                    )
                    continue

                feature = Feature(
                    seqname=parts[0],
                    source=parts[1],
                    feature_type=parts[2],
                    start=int(parts[3]),
                    end=int(parts[4]),
                    score=parts[5],
                    strand=self._normalize_strand(parts[6]),
                    frame=parts[7],
                    attributes=parts[8] if len(parts) > 8 else "",
                    original_format="gff"
                )
                features.append(feature)

            except (ValueError, IndexError) as e:
                self._log_parse_warning(f"Failed to parse GFF line: {e}", line_num)
                continue

        return features

    def _log_parse_warning(self, message: str, line_num: int, content: str = "") -> None:
        """Log a parsing warning."""
        details = {'line': line_num}
        if content:
            details['content'] = content[:100]
        self.logger.add_validation_issue(
            level='WARNING',
            category='feature',
            message=message,
            details=details
        )

    def _parse_bed(self, handle) -> List[Feature]:
        """
        Parse BED format file.

        BED format (tab-separated, minimum 3 columns):
        chrom start end [name] [score] [strand] [...]

        Returns:
            List of Feature objects
        """
        features = []
        for line_num, line in enumerate(handle, start=1):
            line = line.strip()
            if self._skip_line(line):
                continue

            try:
                parts = line.split('\t')
                if len(parts) < 3:
                    self._log_parse_warning(
                        f"Invalid BED line (expected at least 3 columns, got {len(parts)})",
                        line_num, line
                    )
                    continue

                feature = Feature(
                    seqname=parts[0],
                    start=int(parts[1]),
                    end=int(parts[2]),
                    attributes=parts[3] if len(parts) > 3 else "",
                    score=parts[4] if len(parts) > 4 else ".",
                    strand=self._normalize_strand(parts[5]) if len(parts) > 5 else '.',
                    source=".",
                    feature_type="region",
                    frame=".",
                    original_format="bed"
                )
                features.append(feature)

            except (ValueError, IndexError) as e:
                self._log_parse_warning(f"Failed to parse BED line: {e}", line_num)
                continue

        return features

    def _validate_features(self) -> None:
        """Validate parsed features."""
        self.logger.debug("Validating features...")

        for idx, feature in enumerate(self.features):
            self._validate_coordinates(feature, idx)
            self._validate_negative_coords(feature, idx)

        self.logger.debug("✓ Feature validation passed")

    def _validate_coordinates(self, feature: Feature, idx: int) -> None:
        """Validate feature coordinates."""
        if not self.settings.check_coordinates:
            return

        if feature.start > feature.end:
            self._raise_validation_error(
                f"Feature has start > end: {feature.seqname}:{feature.start}-{feature.end}",
                {'feature_index': idx, 'seqname': feature.seqname,
                 'start': feature.start, 'end': feature.end}
            )

        if feature.start == feature.end and not self.settings.allow_zero_length:
            self._raise_validation_error(
                f"Zero-length feature not allowed: {feature.seqname}:{feature.start}-{feature.end}",
                {'feature_index': idx, 'seqname': feature.seqname, 'position': feature.start}
            )

    def _validate_negative_coords(self, feature: Feature, idx: int) -> None:
        """Validate feature doesn't have negative coordinates."""
        if not self.settings.allow_negative_coords and (feature.start < 0 or feature.end < 0):
            self._raise_validation_error(
                f"Negative coordinates not allowed: {feature.seqname}:{feature.start}-{feature.end}",
                {'feature_index': idx, 'seqname': feature.seqname}
            )

    def _raise_validation_error(self, message: str, details: Dict[str, Any]) -> None:
        """Log and raise a validation error."""
        self.logger.add_validation_issue(
            level='ERROR',
            category='feature',
            message=message,
            details=details
        )
        raise FeatureValidationError(message)

    def _parse_gff_attributes(self, attr_string: str) -> Dict[str, str]:
        """
        Parse GFF attributes string into dictionary.

        GFF attributes are semicolon-separated key=value pairs:
        "ID=gene1;Name=testGene;Note=description"

        Args:
            attr_string: GFF attributes string

        Returns:
            Dictionary of attribute key-value pairs

        Example:
            >>> attrs = self._parse_gff_attributes("ID=gene1;Name=test")
            >>> attrs
            {'ID': 'gene1', 'Name': 'test'}
        """
        if not attr_string or attr_string == '.':
            return {}

        attrs = {}
        # Split by semicolon
        for pair in attr_string.split(';'):
            pair = pair.strip()
            if not pair:
                continue

            # Split by first '=' only (values might contain '=')
            if '=' in pair:
                key, value = pair.split('=', 1)
                attrs[key.strip()] = value.strip()
            else:
                # Malformed attribute, store as-is
                self.logger.debug(f"Malformed GFF attribute (no '='): {pair}")

        return attrs

    def _build_gff_attributes(self, attr_dict: Dict[str, str]) -> str:
        """
        Build GFF attributes string from dictionary.

        Creates semicolon-separated key=value pairs, with ID first if present.

        Args:
            attr_dict: Dictionary of attribute key-value pairs

        Returns:
            GFF attributes string

        Example:
            >>> attr_string = self._build_gff_attributes({'ID': 'gene1', 'Name': 'test'})
            >>> attr_string
            'ID=gene1;Name=test'
        """
        if not attr_dict:
            return '.'

        # Build list of key=value pairs, with ID first
        parts = []

        # ID should come first (GFF3 convention)
        if 'ID' in attr_dict:
            parts.append(f"ID={attr_dict['ID']}")

        # Add all other attributes in sorted order (for consistency)
        for key in sorted(attr_dict.keys()):
            if key != 'ID':  # Skip ID, already added
                parts.append(f"{key}={attr_dict[key]}")

        return ';'.join(parts)

    def _apply_edits(self) -> None:
        """
        Apply editing specifications to features based on settings.

        Applies the following edits (if enabled in settings):
        - Sort features by position (sort_by_position)
        - Replace feature seqname/chromosome field (replace_id_with)
        """
        self.logger.debug("Applying editing specifications from settings...")

        # Sort by position (seqname, then start, then end)
        if self.settings.sort_by_position:
            original_count = len(self.features)
            self.features.sort(key=lambda f: (f.seqname, f.start, f.end))
            self.logger.debug(f"Sorted {original_count} features by position")

        # Replace feature seqnames (chromosome/sequence identifiers)
        if self.settings.replace_id_with:
            new_seqname = self.settings.replace_id_with
            self.logger.debug(f"Replacing feature seqnames with '{new_seqname}'...")

            for feature in self.features:
                # Store original seqname in attributes
                old_seqname = feature.seqname

                if feature.original_format == 'gff':
                    # Parse GFF attributes
                    attrs = self._parse_gff_attributes(feature.attributes)

                    # Store original seqname in Note field (GFF3 standard)
                    note_text = f"Original_seqname:{old_seqname}"
                    if 'Note' in attrs and attrs['Note']:
                        # Append to existing Note
                        attrs['Note'] = f"{attrs['Note']};{note_text}"
                    else:
                        # Create new Note
                        attrs['Note'] = note_text

                    # Rebuild attributes string
                    feature.attributes = self._build_gff_attributes(attrs)

                # Replace seqname for both GFF and BED
                feature.seqname = new_seqname

            self.logger.debug(f"Replaced {len(self.features)} feature seqnames with '{new_seqname}'")

        self.logger.debug(f"✓ Edits applied, {len(self.features)} feature(s) remaining")

    def _collect_statistics(self) -> None:
        """Collect statistics about the features."""
        self.logger.debug("Collecting statistics...")

        if not self.features:
            return

        type_counts = {}
        strand_counts = {'+': 0, '-': 0, '.': 0}
        lengths = [f.length for f in self.features]

        for feature in self.features:
            type_counts[feature.feature_type] = type_counts.get(feature.feature_type, 0) + 1
            if feature.strand in strand_counts:
                strand_counts[feature.strand] += 1

        self.statistics = {
            'num_features': len(self.features),
            'num_by_type': type_counts,
            'num_by_strand': strand_counts,
            'total_span': sum(lengths),
            'min_length': min(lengths) if lengths else 0,
            'max_length': max(lengths) if lengths else 0,
            'avg_length': sum(lengths) / len(lengths) if lengths else 0.0
        }

        self.logger.debug(f"Statistics: {self.statistics}")

    def _convert_to_gff(self) -> None:
        """Convert features to GFF format (if needed and requested)."""
        if self.feature_config.detected_format == FeatureFormat.GFF:
            self.logger.debug("Keeping GFF format")
            return

        # Check if we need conversion
        bed_features = [f for f in self.features if f.original_format == 'bed']

        if bed_features:
            self.logger.debug(f"Converting {len(bed_features)} BED features to GFF format...")

            # BED is 0-based half-open, GFF is 1-based closed
            # Convert: BED [start, end) -> GFF [start+1, end]
            for feature in bed_features:
                if feature.original_format == 'bed':
                    feature.start += 1  # Convert from 0-based to 1-based
                    feature.original_format = 'gff'  # Mark as converted

            self.logger.debug("✓ BED to GFF conversion complete")
        else:
            self.logger.debug("All features already in GFF format")

    def _save_output(self) -> Path:
        """
        Save processed features to output directory using settings.

        Returns:
            Path to output file
        """
        self.logger.debug("Saving output file...")

        output_path = self._get_output_path()
        self._write_output_file(output_path)

        self.logger.info(f"Output saved: {output_path}")
        return output_path

    def _get_output_path(self) -> Path:
        """Determine the output file path."""
        output_dir = self.output_dir
        if self.settings.output_subdir_name:
            output_dir = output_dir / self.settings.output_subdir_name
        output_dir.mkdir(parents=True, exist_ok=True)

        base_name = self.feature_config.filename
        for suffix in self.input_path.suffixes:
            base_name = base_name.replace(suffix, '')

        format_ext = '.gff' if any(f.original_format == 'gff' for f in self.features) else '.bed'

        if self.settings.output_filename_suffix:
            output_filename = f"{base_name}_{self.settings.output_filename_suffix}{format_ext}"
        else:
            output_filename = f"{base_name}{format_ext}"

        if self.settings.coding_type == 'gz':
            output_filename += '.gz'
        elif self.settings.coding_type == 'bz2':
            output_filename += '.bz2'

        return output_dir / output_filename

    def _write_output_file(self, output_path: Path) -> None:
        """Write features to output file with appropriate compression."""
        self.logger.debug(f"Writing output to: {output_path}")

        open_func = {
            'gz': lambda p: gzip.open(p, 'wt'),
            'bz2': lambda p: bz2.open(p, 'wt'),
            None: lambda p: open(p, 'w')
        }.get(self.settings.coding_type, lambda p: open(p, 'w'))

        with open_func(output_path) as handle:
            self._write_features(handle)

    def _write_features(self, handle: IO) -> None:
        """
        Write features to file handle in appropriate format.

        Args:
            handle: File handle to write to
        """
        # Determine output format
        use_gff = any(f.original_format == 'gff' for f in self.features)

        if use_gff:
            # Write GFF header
            handle.write("##gff-version 3\n")

            # Write features
            for feature in self.features:
                line = '\t'.join([
                    feature.seqname,
                    feature.source,
                    feature.feature_type,
                    str(feature.start),
                    str(feature.end),
                    feature.score,
                    feature.strand,
                    feature.frame,
                    feature.attributes
                ])
                handle.write(line + '\n')
        else:
            # Write BED format (convert back if needed)
            for feature in self.features:
                # If feature was converted to GFF, convert back to BED coordinates
                start = feature.start - 1 if feature.original_format == 'gff' else feature.start

                # BED minimum 3 columns, write up to 6
                parts = [
                    feature.seqname,
                    str(start),
                    str(feature.end)
                ]

                if feature.attributes:
                    parts.append(feature.attributes)
                if feature.score and feature.score != '.':
                    while len(parts) < 4:
                        parts.append('.')
                    parts.append(feature.score)
                if feature.strand and feature.strand != '.':
                    while len(parts) < 5:
                        parts.append('.')
                    parts.append(feature.strand)

                line = '\t'.join(parts)
                handle.write(line + '\n')

    def get_statistics(self) -> Dict[str, Any]:
        """
        Get statistics about the processed features.

        Returns:
            Dictionary of statistics
        """
        return self.statistics.copy()