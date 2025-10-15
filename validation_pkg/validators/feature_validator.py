"""
Feature file validator and processor.

Handles GFF, GTF, and BED formats with compression support.
"""

from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
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
    FileNotFoundError as ValidationFileNotFoundError
)
from validation_pkg.utils.formats import FeatureFormat, CodingType


@dataclass
class Feature:
    """Represents a genomic feature (generic container for GFF/BED data)."""
    seqname: str           # Sequence/chromosome name
    start: int             # Start position (1-based for GFF, 0-based for BED)
    end: int               # End position
    feature_type: str = "" # Feature type (gene, CDS, etc.)
    score: str = "."       # Score
    strand: str = "+"      # Strand (+, -, .)
    source: str = "."      # Source
    frame: str = "."       # Frame (for CDS)
    attributes: str = ""   # Attributes column (GFF) or name (BED)

    # Store original format to preserve information
    original_format: str = "gff"  # 'gff' or 'bed'


class FeatureValidator:
    """
    Validates and processes feature annotation files (GFF, GTF, BED formats).

    Workflow:
    1. Detect and handle compression
    2. Parse features from file
    3. Validate feature coordinates and structure
    4. Apply editing specifications
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
            coding_type: Output compression type: 'gz', 'bz2', or None (default: None)
            output_filename_suffix: Suffix to add to output filename (default: None)
            output_subdir_name: Subdirectory name for output files (default: None)

        Example:
            >>> settings = FeatureValidator.Settings()
            >>> settings = settings.update(sort_by_position=True, coding_type='gz')
            >>> validator = FeatureValidator(feature_config, output_dir, settings)
        """
        # Validation options
        sort_by_position: bool = True
        check_coordinates: bool = True
        allow_zero_length: bool = False
        allow_negative_coords: bool = True

        # Output format
        coding_type: Optional[str] = None
        output_filename_suffix: Optional[str] = None
        output_subdir_name: Optional[str] = None


    def __init__(self, feature_config, output_dir, settings: Optional[Settings] = None):
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
    
    def validate(self):
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

    def _open_file(self, mode='rt'):
        """
        Open file with automatic decompression based on feature_config.

        Args:
            mode: File opening mode (default: 'rt' for text read)

        Returns:
            File handle

        Raises:
            CompressionError: If file cannot be opened or decompressed
        """
        try:
            coding_type_str = str(self.feature_config.coding_type).lower()

            if 'gzip' in coding_type_str or coding_type_str == 'gz':
                return gzip.open(self.input_path, mode)
            elif 'bzip2' in coding_type_str or coding_type_str == 'bz2':
                return bz2.open(self.input_path, mode)
            else:
                # No compression or unrecognized
                return open(self.input_path, mode)
        except Exception as e:
            error_msg = f"Failed to open file: {e}"
            self.logger.add_validation_issue(
                level='ERROR',
                category='feature',
                message=error_msg,
                details={'file': str(self.input_path), 'error': str(e)}
            )
            raise CompressionError(error_msg) from e

    def _parse_file(self):
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

    def _parse_gff(self, handle) -> List[Feature]:
        """
        Parse GFF/GTF format file.

        GFF format (tab-separated):
        seqname source feature start end score strand frame attributes

        Returns:
            List of Feature objects
        """
        features = []
        line_num = 0

        for line in handle:
            line_num += 1
            line = line.strip()

            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue

            # Parse GFF line
            try:
                parts = line.split('\t')
                if len(parts) < 8:
                    self.logger.add_validation_issue(
                        level='WARNING',
                        category='feature',
                        message=f"Invalid GFF line (expected 8-9 columns, got {len(parts)})",
                        details={'line': line_num, 'content': line[:100]}
                    )
                    continue

                feature = Feature(
                    seqname=parts[0],
                    source=parts[1],
                    feature_type=parts[2],
                    start=int(parts[3]),  # GFF is 1-based
                    end=int(parts[4]),
                    score=parts[5],
                    strand=parts[6] if parts[6] in ['+', '-', '.'] else '.',
                    frame=parts[7],
                    attributes=parts[8] if len(parts) > 8 else "",
                    original_format="gff"
                )
                features.append(feature)

            except (ValueError, IndexError) as e:
                self.logger.add_validation_issue(
                    level='WARNING',
                    category='feature',
                    message=f"Failed to parse GFF line: {e}",
                    details={'line': line_num, 'error': str(e)}
                )
                continue

        return features

    def _parse_bed(self, handle) -> List[Feature]:
        """
        Parse BED format file.

        BED format (tab-separated, minimum 3 columns):
        chrom start end [name] [score] [strand] [...]

        Returns:
            List of Feature objects
        """
        features = []
        line_num = 0

        for line in handle:
            line_num += 1
            line = line.strip()

            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue

            # Parse BED line
            try:
                parts = line.split('\t')
                if len(parts) < 3:
                    self.logger.add_validation_issue(
                        level='WARNING',
                        category='feature',
                        message=f"Invalid BED line (expected at least 3 columns, got {len(parts)})",
                        details={'line': line_num, 'content': line[:100]}
                    )
                    continue

                feature = Feature(
                    seqname=parts[0],
                    start=int(parts[1]),  # BED is 0-based
                    end=int(parts[2]),
                    attributes=parts[3] if len(parts) > 3 else "",  # name column
                    score=parts[4] if len(parts) > 4 else ".",
                    strand=parts[5] if len(parts) > 5 and parts[5] in ['+', '-', '.'] else '.',
                    source=".",
                    feature_type="region",  # BED doesn't have type, use generic
                    frame=".",
                    original_format="bed"
                )
                features.append(feature)

            except (ValueError, IndexError) as e:
                self.logger.add_validation_issue(
                    level='WARNING',
                    category='feature',
                    message=f"Failed to parse BED line: {e}",
                    details={'line': line_num, 'error': str(e)}
                )
                continue

        return features

    def _validate_features(self):
        """Validate parsed features."""
        self.logger.debug("Validating features...")

        for idx, feature in enumerate(self.features):
            # Check coordinates
            if self.settings.check_coordinates:
                if feature.start > feature.end:
                    error_msg = f"Feature has start > end: {feature.seqname}:{feature.start}-{feature.end}"
                    self.logger.add_validation_issue(
                        level='ERROR',
                        category='feature',
                        message=error_msg,
                        details={
                            'feature_index': idx,
                            'seqname': feature.seqname,
                            'start': feature.start,
                            'end': feature.end
                        }
                    )
                    raise FeatureValidationError(error_msg)

                if feature.start == feature.end and not self.settings.allow_zero_length:
                    error_msg = f"Zero-length feature not allowed: {feature.seqname}:{feature.start}-{feature.end}"
                    self.logger.add_validation_issue(
                        level='ERROR',
                        category='feature',
                        message=error_msg,
                        details={
                            'feature_index': idx,
                            'seqname': feature.seqname,
                            'position': feature.start
                        }
                    )
                    raise FeatureValidationError(error_msg)

            # Check for negative coordinates
            if not self.settings.allow_negative_coords and (feature.start < 0 or feature.end < 0):
                error_msg = f"Negative coordinates not allowed: {feature.seqname}:{feature.start}-{feature.end}"
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='feature',
                    message=error_msg,
                    details={'feature_index': idx, 'seqname': feature.seqname}
                )
                raise FeatureValidationError(error_msg)

        self.logger.debug("✓ Feature validation passed")

    def _apply_edits(self):
        """
        Apply editing specifications to features based on settings.

        Applies the following edits (if enabled in settings):
        - Sort features by position (sort_by_position)
        """
        self.logger.debug("Applying editing specifications from settings...")

        # Sort by position (seqname, then start, then end)
        if self.settings.sort_by_position:
            original_count = len(self.features)
            self.features.sort(key=lambda f: (f.seqname, f.start, f.end))
            self.logger.debug(f"Sorted {original_count} features by position")

        self.logger.debug(f"✓ Edits applied, {len(self.features)} feature(s) remaining")

    def _collect_statistics(self):
        """Collect statistics about the features."""
        self.logger.debug("Collecting statistics...")

        if not self.features:
            return

        # Count by type
        type_counts = {}
        strand_counts = {'+': 0, '-': 0, '.': 0}
        lengths = []

        for feature in self.features:
            # Count by type
            ftype = feature.feature_type
            type_counts[ftype] = type_counts.get(ftype, 0) + 1

            # Count by strand
            if feature.strand in strand_counts:
                strand_counts[feature.strand] += 1

            # Calculate length
            length = feature.end - feature.start
            lengths.append(length)

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

    def _convert_to_gff(self):
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

        # Determine output directory (with optional subdirectory)
        output_dir = self.output_dir
        if self.settings.output_subdir_name:
            output_dir = output_dir / self.settings.output_subdir_name

        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate output filename from original filename (without compression extension)
        base_name = self.feature_config.filename
        # Remove all suffixes
        for suffix in self.input_path.suffixes:
            base_name = base_name.replace(suffix, '')

        # Determine output format extension
        if any(f.original_format == 'gff' for f in self.features):
            format_ext = '.gff'
        else:
            format_ext = '.bed'

        # Add suffix if specified
        if self.settings.output_filename_suffix:
            output_filename = f"{base_name}_{self.settings.output_filename_suffix}{format_ext}"
        else:
            output_filename = f"{base_name}{format_ext}"

        # Add compression extension if requested
        if self.settings.coding_type == 'gz':
            output_filename += '.gz'
        elif self.settings.coding_type == 'bz2':
            output_filename += '.bz2'

        output_path = output_dir / output_filename

        # Write output with appropriate compression
        self.logger.debug(f"Writing output to: {output_path}")

        if self.settings.coding_type == 'gz':
            with gzip.open(output_path, 'wt') as handle:
                self._write_features(handle)
        elif self.settings.coding_type == 'bz2':
            with bz2.open(output_path, 'wt') as handle:
                self._write_features(handle)
        else:
            with open(output_path, 'w') as handle:
                self._write_features(handle)

        self.logger.info(f"Output saved: {output_path}")

        return output_path

    def _write_features(self, handle):
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