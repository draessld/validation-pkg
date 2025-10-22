"""
Feature file validator and processor.

Handles GFF, GTF, and BED formats with compression support.
"""

from pathlib import Path
from typing import Optional, Dict, Any, List, IO
from dataclasses import dataclass
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
)
from validation_pkg.utils.formats import CodingType as CT
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
    6. Convert to GFF format (if needed)
    7. Save to output directory
    """

    @dataclass
    class Settings(BaseSettings):
        """
        Settings for feature validation and processing.

        Attributes:
            validation_level: Validation thoroughness level:
                - 'strict': Full parsing and validation of all features (default)
                - 'trust': Fast validation - check basic format and first few features
                - 'minimal': Only verify file exists and has correct extension
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

        Validation Level Behavior:
            - validation_level='strict': Parse all features, validate all, apply edits to all
            - validation_level='trust': Parse all features, validate first 10 only, apply edits to all
            - validation_level='minimal': Skip parsing, copy file as-is

        Example:
            >>> settings = FeatureValidator.Settings()
            >>> settings = settings.update(sort_by_position=True, coding_type='gz', replace_id_with='chr1')
            >>> validator = FeatureValidator(feature_config, output_dir, settings)
        """
        # Validation level
        validation_level: str = 'strict'  # 'strict', 'trust', or 'minimal'

        # Validation options
        sort_by_position: bool = True
        check_coordinates: bool = True

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

        # Validate validation_level setting
        valid_levels = {'strict', 'trust', 'minimal'}
        if self.settings.validation_level not in valid_levels:
            raise ValueError(
                f"Invalid validation_level '{self.settings.validation_level}'. "
                f"Must be one of: {', '.join(valid_levels)}"
            )

        # Log settings being used
        self.logger.debug(f"Initializing FeatureValidator with settings:\n{self.settings}")

        # Resolved paths
        self.input_path = feature_config.filepath

        # Parsed data
        self.features: List[Feature] = []

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
            self._parse_file()

            self._validate_features()

            self._convert_to_gff()

            self._apply_edits()

            self._save_output()

            self.logger.info(f"✓ Feature validation completed")

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
        """
        Parse feature file based on detected format from feature_config.

        Behavior depends on validation_level:
        - 'strict': Full parsing and validation of all features
        - 'trust': Parse all features, but validate only first 10
        - 'minimal': Skip parsing entirely
        """
        format_str = str(self.feature_config.detected_format)
        self.logger.debug(f"Parsing {format_str} file (validation_level={self.settings.validation_level})...")

        # Minimal mode - skip parsing
        if self.settings.validation_level == 'minimal':
            self.logger.info("Minimal validation mode - skipping file parsing")
            self.features = []
            return

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

            if self.settings.validation_level == 'trust':
                self.logger.info(f"Trust mode - parsed {len(self.features)} feature(s), will validate first 10 only")
            else:
                self.logger.debug(f"Parsed {len(self.features)} feature(s)")

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

        Parses all features regardless of validation_level.

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
                )
                features.append(feature)

            except (ValueError, IndexError) as e:
                self._log_parse_warning(f"Failed to parse GFF line: {e}", line_num)
                continue

        return features

    def _parse_bed(self, handle) -> List[Feature]:
        """
        Parse BED format file.

        BED format (tab-separated, minimum 3 columns):
        chrom start end [name] [score] [strand] [...]

        Parses all features regardless of validation_level.

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
                )
                features.append(feature)

            except (ValueError, IndexError) as e:
                self._log_parse_warning(f"Failed to parse BED line: {e}", line_num)
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

    def _validate_features(self) -> None:
        """
        Validate parsed features.

        Behavior depends on validation_level:
        - 'strict': Validate all features
        - 'trust': Validate only first 10 features
        - 'minimal': Skip validation (no parsing done)
        """
        self.logger.debug("Validating features...")

        # Minimal mode - no features to validate
        if self.settings.validation_level == 'minimal':
            self.logger.debug("Minimal mode - skipping feature validation")
            return

        # Determine how many features to validate
        if self.settings.validation_level == 'trust':
            # Trust mode - validate only first 10 features
            features_to_validate = min(10, len(self.features))
            self.logger.debug(f"Trust mode - validating first {features_to_validate} of {len(self.features)} features")
        else:
            # Strict mode - validate all features
            features_to_validate = len(self.features)

        # Validate features
        if self.settings.check_coordinates:
            for idx in range(features_to_validate):
                self._validate_coordinates(self.features[idx], idx)

        self.logger.debug("✓ Feature validation passed")

    def _validate_coordinates(self, feature: Feature, idx: int) -> None:
        """Validate feature coordinates."""
        if feature.start > feature.end:
            error_message = f"Feature has start > end: {feature.seqname}:{feature.start}-{feature.end}"
            self.logger.add_validation_issue(
                level='ERROR',
                category='feature',
                message=error_message,
                details={'feature_index': idx, 'seqname': feature.seqname, 'start': feature.start, 'end': feature.end}
            )
            raise FeatureValidationError(error_message)

        if feature.start < 0 or feature.end < 0:
            error_message = f"Negative coordinates not allowed: {feature.seqname}:{feature.start}-{feature.end}"
            self.logger.add_validation_issue(
                level='ERROR',
                category='feature',
                message=error_message,
                details={'feature_index': idx, 'seqname': feature.seqname}
            )
            raise FeatureValidationError(error_message)

        if feature.start == feature.end:
            error_message = f"Zero-length feature not allowed: {feature.seqname}:{feature.start}-{feature.end}"
            self.logger.add_validation_issue(
                level='ERROR',
                category='feature',
                message=error_message,
                details={'feature_index': idx, 'seqname': feature.seqname, 'position': feature.start}
            )
            raise FeatureValidationError(error_message)

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

        Behavior depends on validation_level:
        - 'strict': Apply all edits to all features
        - 'trust': Apply all edits to parsed features (first 10)
        - 'minimal': Skip edits (file will be copied as-is)
        """
        self.logger.debug("Applying editing specifications from settings...")

        # Minimal mode - skip edits, file will be copied as-is
        if self.settings.validation_level == 'minimal':
            self.logger.debug("Minimal mode - skipping edits")
            return

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

    def _convert_to_gff(self) -> None:
        """Convert features to GFF format (if needed and requested)."""
        if self.feature_config.detected_format == FeatureFormat.GFF:
            self.logger.debug("Keeping GFF format")
            return

        self.logger.debug(f"Converting {len(self.features)} BED features to GFF format...")

        # BED is 0-based half-open, GFF is 1-based closed
        # Convert: BED [start, end) -> GFF [start+1, end]
        for feature in self.features:
            feature.start += 1  # Convert from 0-based to 1-based

        self.logger.debug("✓ BED to GFF conversion complete")

    def _save_output(self) -> Path:
        """
        Save processed features to output directory using settings.

        Behavior depends on validation_level:
        - 'strict': Write features using parser (may convert format)
        - 'trust': Write parsed features (first 10 only) - NOT RECOMMENDED, use strict instead
        - 'minimal': Copy file as-is (preserve original format)

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
        base_name = self.feature_config.filename

        # Remove all suffixes (.gff.gz -> remove .gz then .gff)
        for suffix in self.input_path.suffixes:
            base_name = base_name.replace(suffix, '')
        # Add suffix if specified
        if self.settings.output_filename_suffix:
            output_filename = f"{base_name}_{self.settings.output_filename_suffix}.gff"
        else:
            output_filename = f"{base_name}.gff"

        # Add compression extension if requested (from settings, not input)
        # Support both string ('gz', 'bz2') and enum (CodingType.GZIP, CodingType.BZIP2)
        coding = self.settings.coding_type

        if coding in ('gz', 'gzip', CT.GZIP):
            output_filename += '.gz'
        elif coding in ('bz2', 'bzip2', CT.BZIP2):
            output_filename += '.bz2'

        # Minimal mode - copy file as-is without parsing
        # Required: GFF format + coding must match settings.coding_type
        if self.settings.validation_level == 'minimal':
            self.logger.debug("Minimal mode - validating format and coding requirements")

            # Check format - must be GFF (reject BED files)
            if self.feature_config.detected_format != FeatureFormat.GFF:
                error_msg = f'Minimal mode requires GFF format, got {self.feature_config.detected_format}. Use validation_level "trust" or "strict" to convert.'
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='feature',
                    message=error_msg,
                    details={'file': self.feature_config.filename, 'detected_format': str(self.feature_config.detected_format)}
                )
                raise FeatureValidationError(error_msg)

            # Check coding - must match settings.coding_type
            # Normalize both for comparison
            input_coding = self.feature_config.coding_type
            required_coding = coding  # This is already normalized from settings.coding_type

            # Map None to CT.NONE for comparison
            if input_coding is None:
                input_coding = CT.NONE
            if required_coding is None:
                required_coding = CT.NONE

            # Convert string to CT if needed
            if isinstance(required_coding, str):
                if required_coding in ('gz', 'gzip'):
                    required_coding = CT.GZIP
                elif required_coding in ('bz2', 'bzip2'):
                    required_coding = CT.BZIP2

            if input_coding != required_coding:
                error_msg = f'Minimal mode requires input coding to match output coding. Input: {input_coding}, Required: {required_coding}. Use validation_level "trust" or "strict" to change compression.'
                self.logger.add_validation_issue(
                    level='ERROR',
                    category='feature',
                    message=error_msg,
                    details={'file': self.feature_config.filename, 'input_coding': str(input_coding), 'required_coding': str(required_coding)}
                )
                raise FeatureValidationError(error_msg)

            # Use the output_filename already constructed (includes correct extension)
            output_path_minimal = output_dir / output_filename

            # Copy file
            self.logger.debug(f"Copying {self.input_path} to {output_path_minimal}")
            shutil.copy2(self.input_path, output_path_minimal)

            self.logger.info(f"Output saved: {output_path_minimal}")
            return output_path_minimal

        # Strict and Trust modes - write edited features
        output_path = output_dir / output_filename
        self.logger.debug(f"Writing output to: {output_path}")

        open_func = {
            'gz': lambda p: gzip.open(p, 'wt'),
            'bz2': lambda p: bz2.open(p, 'wt'),
            None: lambda p: open(p, 'w')
        }.get(self.settings.coding_type, lambda p: open(p, 'w'))

        with open_func(output_path) as handle:
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

        self.logger.info(f"Output saved: {output_path}")
        return output_path