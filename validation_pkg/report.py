"""
Validation report generation module.

Provides comprehensive reporting for validation results including
input settings, output metadata, and inter-file validation results.
"""

from dataclasses import dataclass, field
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional, Any, Union
import json


@dataclass
class FileValidationRecord:
    """Complete record of a single file's validation."""
    # Input information
    input_file: str
    validator_type: str  # "genome", "read", "feature"
    input_settings: Optional[dict] = None  # Serialized Settings object

    # Output information
    output_file: str = None
    output_metadata: dict = field(default_factory=dict)

    # Timing
    elapsed_time: Optional[float] = None
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())


@dataclass
class InterFileValidationRecord:
    """Record of inter-file validation check."""
    validation_type: str  # "genomexgenome", "readxread", etc.
    status: str  # "PASSED", "FAILED"
    passed: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())


class ValidationReport:
    """
    Comprehensive validation report builder.

    Collects results from validators and inter-file validation,
    tracks input settings and output metadata, then generates a
    detailed report when flush() is called.

    Features:
    - Input settings tracking (what parameters were used)
    - Output metadata tracking (what was produced)
    - Inter-file validation results
    - Validation issues (errors/warnings per file)
    - Timing information
    - Text and JSON export formats

    Usage:
        report = ValidationReport(Path("report.txt"))

        # Add validator results
        genome_res = validate_genome(config.ref_genome, settings)
        report.write(genome_res, file_type="genome",
                    input_file=config.ref_genome.filename,
                    settings=settings)

        # Add inter-file validation results
        interfile_res = genomexgenome_validation(ref_res, mod_res, settings)
        report.write(interfile_res, file_type="genomexgenome")

        # Generate final report
        report.flush()  # Writes text report
        report.flush(format="json")  # Writes JSON report
    """

    def __init__(self, report_path: Path):
        """
        Initialize validation report.

        Args:
            report_path: Path where the report will be written
        """
        self.report_path = Path(report_path)
        self.report_path.parent.mkdir(parents=True, exist_ok=True)

        # Storage for file validation records
        self.file_records: List[FileValidationRecord] = []

        # Storage for inter-file validation records
        self.interfile_records: List[InterFileValidationRecord] = []

        # Track start time
        self.start_time = datetime.now()

    def _metadata_to_dict(self, metadata_obj: Any) -> Dict[str, Any]:
        """
        Convert OutputMetadata object to dictionary.

        Args:
            metadata_obj: OutputMetadata object from a validator

        Returns:
            Dictionary representation of the metadata
        """
        if hasattr(metadata_obj, '__dict__'):
            # Convert dataclass to dict, filtering out None values for cleaner output
            return {k: v for k, v in metadata_obj.__dict__.items() if v is not None}
        return {}

    def write(
        self,
        result: Union[dict, List[dict]],
        file_type: str,
        input_file: Optional[str] = None,
        settings: Optional[Any] = None,
        elapsed_time: Optional[float] = None
    ) -> None:
        """
        Add a validation result to the report.

        Args:
            result: Result dict from validator or inter-file validation
            file_type: Type of result:
                - "genome": GenomeValidator result
                - "read": ReadValidator result (can be single or list)
                - "feature": FeatureValidator result
                - "genomexgenome": Genome inter-file validation
                - "readxread": Read inter-file validation
                - "featurexgenome": Feature-genome inter-file validation
            input_file: Original input filename (for file validators)
            settings: Settings object used for validation (for file validators)
            elapsed_time: Time taken for validation (optional)
        """
        # Handle file validation results
        if file_type in ("genome", "read", "feature"):
            # Handle both single result and list of results
            results_list = result if isinstance(result, list) else [result]

            for single_result in results_list:
                # Handle both OutputMetadata objects and dicts (backward compatibility)
                if hasattr(single_result, 'output_file'):
                    # New OutputMetadata object
                    output_file = single_result.output_file
                    output_metadata_dict = self._metadata_to_dict(single_result)
                else:
                    # Old dict format
                    output_file = single_result.get('output_file', 'unknown')
                    output_metadata_dict = single_result.get('metadata') or single_result.get('output_metadata', {})

                # Extract input filename if not provided
                if input_file is None:
                    input_file = Path(output_file).name if output_file else 'unknown'

                # Serialize settings if provided
                settings_dict = None
                if settings is not None:
                    settings_dict = settings.to_dict() if hasattr(settings, 'to_dict') else {}

                # Create file record
                record = FileValidationRecord(
                    input_file=input_file,
                    validator_type=file_type,
                    input_settings=settings_dict,
                    output_file=output_file,
                    output_metadata=output_metadata_dict,
                    elapsed_time=elapsed_time
                )

                self.file_records.append(record)

        # Handle inter-file validation results
        elif file_type in ("genomexgenome", "readxread", "featurexgenome"):
            passed = result.get('passed', False)
            status = "PASSED" if passed else "FAILED"

            record = InterFileValidationRecord(
                validation_type=file_type,
                status=status,
                passed=passed,
                errors=result.get('errors', []),
                warnings=result.get('warnings', []),
                metadata=result.get('metadata', {})
            )

            self.interfile_records.append(record)
        else:
            raise ValueError(f"Unknown file_type: {file_type}")

    def flush(self, format: str = "text") -> None:
        """
        Generate and write the final report to file.

        Args:
            format: Report format - "text" (default) or "json"

        Creates a comprehensive report with:
        - Summary statistics
        - Per-file validation results with settings
        - Inter-file validation results
        - Issue listings
        """
        if format == "text":
            self._flush_text()
        elif format == "json":
            self._flush_json()
        else:
            raise ValueError(f"Unknown format: {format}. Use 'text' or 'json'")

    def _flush_text(self) -> None:
        """Generate and write text format report."""
        end_time = datetime.now()
        duration = (end_time - self.start_time).total_seconds()

        lines = []

        # Header
        lines.append("=" * 100)
        lines.append("BIOINFORMATICS VALIDATION REPORT".center(100))
        lines.append("=" * 100)
        lines.append(f"Generated: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append(f"Total Duration: {duration:.2f}s")
        lines.append("")

        # Summary
        lines.extend(self._generate_summary_section())

        # File-specific results
        if self.file_records:
            lines.extend(self._generate_file_results_section())

        # Inter-file validation results
        if self.interfile_records:
            lines.extend(self._generate_interfile_section())

        # Write to file
        report_text = '\n'.join(lines)
        self.report_path.write_text(report_text, encoding='utf-8')

        print(f"\n✓ Report written to: {self.report_path}")

    def _flush_json(self) -> None:
        """Generate and write JSON format report."""
        end_time = datetime.now()
        duration = (end_time - self.start_time).total_seconds()

        report_data = {
            "report_metadata": {
                "generated": end_time.isoformat(),
                "duration_seconds": duration,
                "version": "1.0.0"
            },
            "summary": self._get_summary_data(),
            "file_validations": [
                {
                    "input_file": rec.input_file,
                    "validator_type": rec.validator_type,
                    "input_settings": rec.input_settings,
                    "output_file": rec.output_file,
                    "output_metadata": rec.output_metadata,
                    "elapsed_time": rec.elapsed_time,
                    "timestamp": rec.timestamp
                }
                for rec in self.file_records
            ],
            "inter_file_validations": [
                {
                    "validation_type": rec.validation_type,
                    "status": rec.status,
                    "passed": rec.passed,
                    "errors": rec.errors,
                    "warnings": rec.warnings,
                    "metadata": rec.metadata,
                    "timestamp": rec.timestamp
                }
                for rec in self.interfile_records
            ]
        }

        json_path = self.report_path.with_suffix('.json')
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(report_data, f, indent=2, ensure_ascii=False)

        print(f"\n✓ JSON report written to: {json_path}")

    def _generate_summary_section(self) -> List[str]:
        """Generate summary section of report."""
        lines = []

        # Count files by type
        genome_count = sum(1 for r in self.file_records if r.validator_type == "genome")
        read_count = sum(1 for r in self.file_records if r.validator_type == "read")
        feature_count = sum(1 for r in self.file_records if r.validator_type == "feature")
        total_files = len(self.file_records)

        # Count inter-file validations
        interfile_passed = sum(1 for r in self.interfile_records if r.passed)
        interfile_failed = len(self.interfile_records) - interfile_passed

        # Overall status
        overall_status = "✓ PASSED" if interfile_failed == 0 else "✗ FAILED"

        lines.append("SUMMARY")
        lines.append("=" * 100)
        lines.append(f"Overall Status: {overall_status}")
        lines.append("")
        lines.append(f"Files Processed: {total_files}")
        lines.append(f"  Genomes:  {genome_count}")
        lines.append(f"  Reads:    {read_count}")
        lines.append(f"  Features: {feature_count}")
        lines.append("")
        lines.append(f"Inter-file Validations: {len(self.interfile_records)}")
        lines.append(f"  Passed:  {interfile_passed}")
        lines.append(f"  Failed:  {interfile_failed}")
        lines.append("")

        return lines

    def _get_summary_data(self) -> Dict[str, Any]:
        """Get summary data for JSON export."""
        genome_count = sum(1 for r in self.file_records if r.validator_type == "genome")
        read_count = sum(1 for r in self.file_records if r.validator_type == "read")
        feature_count = sum(1 for r in self.file_records if r.validator_type == "feature")

        interfile_passed = sum(1 for r in self.interfile_records if r.passed)
        interfile_failed = len(self.interfile_records) - interfile_passed

        return {
            "total_files": len(self.file_records),
            "genome_files": genome_count,
            "read_files": read_count,
            "feature_files": feature_count,
            "interfile_validations": len(self.interfile_records),
            "interfile_passed": interfile_passed,
            "interfile_failed": interfile_failed,
            "overall_status": "PASSED" if interfile_failed == 0 else "FAILED"
        }

    def _generate_file_results_section(self) -> List[str]:
        """Generate file-specific results section."""
        lines = []
        lines.append("=" * 100)
        lines.append("FILE VALIDATION RESULTS".center(100))
        lines.append("=" * 100)
        lines.append("")

        for idx, record in enumerate(self.file_records, 1):
            lines.extend(self._format_file_record(idx, record))

        return lines

    def _format_file_record(self, idx: int, record: FileValidationRecord) -> List[str]:
        """Format a single file validation record."""
        lines = []

        # Header with file type and name
        file_type_upper = record.validator_type.upper()
        output_name = Path(record.output_file).name if record.output_file else record.input_file
        lines.append(f"[{idx}] {file_type_upper}: {output_name}")
        lines.append("-" * 100)

        # Input information
        lines.append(f"Input File:  {record.input_file}")

        # Input settings
        if record.input_settings:
            lines.append("Input Settings:")
            for key, value in sorted(record.input_settings.items()):
                lines.append(f"  {key}: {value}")

        # Output information
        if record.output_file:
            lines.append(f"Output File: {record.output_file}")

        # Output metadata (validator-specific)
        if record.output_metadata:
            lines.append("Output Metadata:")
            lines.extend(self._format_metadata(record.output_metadata, record.validator_type))

        # Timing
        if record.elapsed_time is not None:
            lines.append(f"Processing Time: {record.elapsed_time:.2f}s")

        lines.append("")
        return lines

    def _format_metadata(self, metadata: Dict[str, Any], validator_type: str) -> List[str]:
        """Format metadata based on validator type."""
        lines = []

        if validator_type == "genome":
            # Genome-specific metadata
            if 'num_sequences' in metadata:
                lines.append(f"  Sequences: {metadata['num_sequences']}")

            if 'sequence_ids' in metadata:
                seq_ids = metadata['sequence_ids']
                if len(seq_ids) <= 5:
                    lines.append(f"    IDs: {', '.join(seq_ids)}")
                else:
                    lines.append(f"    IDs: {', '.join(seq_ids[:5])}, ... ({len(seq_ids)} total)")

            if 'sequence_lengths' in metadata:
                # Handle both list and dict formats
                seq_lengths = metadata['sequence_lengths']
                if isinstance(seq_lengths, dict):
                    total_len = sum(seq_lengths.values())
                else:  # Assume list
                    total_len = sum(seq_lengths)
                lines.append(f"  Total Length: {total_len:,} bp")

        elif validator_type == "read":
            # Read-specific metadata
            if 'num_reads' in metadata:
                lines.append(f"  Reads: {metadata['num_reads']:,}")

            if 'ngs_type_detected' in metadata:
                lines.append(f"  NGS Type: {metadata['ngs_type_detected']}")

            if 'base_name' in metadata and 'read_number' in metadata:
                if metadata['base_name'] and metadata['read_number']:
                    lines.append(f"  Paired-End: R{metadata['read_number']} (base: {metadata['base_name']})")

        elif validator_type == "feature":
            # Feature-specific metadata
            if 'num_features' in metadata:
                lines.append(f"  Features: {metadata['num_features']}")

            if 'feature_types' in metadata:
                types = metadata['feature_types']
                lines.append(f"    Types: {', '.join(types)}")

            if 'sequence_ids' in metadata:
                seq_ids = metadata['sequence_ids']
                if len(seq_ids) <= 5:
                    lines.append(f"    Sequences: {', '.join(seq_ids)}")
                else:
                    lines.append(f"    Sequences: {', '.join(seq_ids[:5])}, ... ({len(seq_ids)} total)")

        return lines

    def _generate_interfile_section(self) -> List[str]:
        """Generate inter-file validation section."""
        lines = []
        lines.append("=" * 100)
        lines.append("INTER-FILE VALIDATION RESULTS".center(100))
        lines.append("=" * 100)
        lines.append("")

        for idx, record in enumerate(self.interfile_records, 1):
            lines.extend(self._format_interfile_record(idx, record))

        return lines

    def _format_interfile_record(self, idx: int, record: InterFileValidationRecord) -> List[str]:
        """Format a single inter-file validation record."""
        lines = []

        # Title
        if record.validation_type == "genomexgenome":
            title = "Genome-to-Genome Validation"
        elif record.validation_type == "readxread":
            title = "Read-to-Read Validation (Paired-End Completeness)"
        elif record.validation_type == "featurexgenome":
            title = "Feature-to-Genome Validation"
        else:
            title = f"Inter-file Validation ({record.validation_type})"

        lines.append(f"[{idx}] {title}")
        lines.append("-" * 100)

        # Status with visual indicator
        status_symbol = "✓" if record.passed else "✗"
        lines.append(f"Status: {status_symbol} {record.status}")

        # Errors
        if record.errors:
            lines.append(f"\nErrors ({len(record.errors)}):")
            for error in record.errors:
                lines.append(f"  • {error}")

        # Warnings
        if record.warnings:
            lines.append(f"\nWarnings ({len(record.warnings)}):")
            for warning in record.warnings:
                lines.append(f"  • {warning}")

        # Details/Metadata
        if record.metadata:
            lines.append("\nDetails:")
            for key, value in sorted(record.metadata.items()):
                # Format lists nicely
                if isinstance(value, list):
                    if len(value) == 0:
                        lines.append(f"  {key}: (none)")
                    elif len(value) <= 5:
                        lines.append(f"  {key}: {value}")
                    else:
                        lines.append(f"  {key}: {len(value)} items")
                # Format dicts
                elif isinstance(value, dict):
                    if len(str(value)) > 100:
                        lines.append(f"  {key}: {len(value)} entries")
                    else:
                        lines.append(f"  {key}: {value}")
                else:
                    lines.append(f"  {key}: {value}")

        lines.append("")
        return lines
