"""
Validation Coordinator - Orchestrates the entire validation workflow.

This module provides the main entry point for running complete validation workflows,
coordinating all validators and generating comprehensive reports.
"""

from pathlib import Path
from typing import Dict, Any, Optional, List
from dataclasses import dataclass

from validation_pkg.config_manager import ConfigManager, Config
from validation_pkg.validators.genome_validator import GenomeValidator
from validation_pkg.validators.read_validator import ReadValidator
from validation_pkg.validators.feature_validator import FeatureValidator
from validation_pkg.logger import setup_logging, get_logger
from validation_pkg.exceptions import ValidationError


@dataclass
class ValidationReport:
    """
    Comprehensive report of validation results.

    Attributes:
        passed: True if all validations passed
        errors: Number of validation errors
        warnings: Number of validation warnings
        validated_files: List of successfully validated files
        failed_files: List of files that failed validation
        statistics: Combined statistics from all validators
        output_dir: Directory containing output files
    """
    passed: bool
    errors: int
    warnings: int
    validated_files: List[str]
    failed_files: List[Dict[str, str]]  # filename -> error message
    statistics: Dict[str, Any]
    output_dir: Path

    def summary(self) -> str:
        """Generate human-readable summary of validation results."""
        lines = []
        lines.append("=" * 60)
        lines.append("VALIDATION REPORT")
        lines.append("=" * 60)
        lines.append(f"Status: {'✓ PASSED' if self.passed else '✗ FAILED'}")
        lines.append(f"Errors: {self.errors}")
        lines.append(f"Warnings: {self.warnings}")
        lines.append(f"Successfully validated: {len(self.validated_files)} file(s)")

        if self.validated_files:
            lines.append("\nValidated Files:")
            for filename in self.validated_files:
                lines.append(f"  ✓ {filename}")

        if self.failed_files:
            lines.append(f"\nFailed: {len(self.failed_files)} file(s)")
            for item in self.failed_files:
                lines.append(f"  ✗ {item['filename']}: {item['error']}")

        lines.append(f"\nOutput directory: {self.output_dir}")
        lines.append("=" * 60)

        return "\n".join(lines)


class ValidationCoordinator:
    """
    Orchestrates the entire validation workflow.

    Coordinates validation of genomes, reads, and features, providing
    a unified interface for running complete validation pipelines.

    Example:
        >>> coordinator = ValidationCoordinator("config.json")
        >>> report = coordinator.validate_all()
        >>> if report.passed:
        ...     print("All validations passed!")
        >>> else:
        ...     print(f"Found {report.errors} errors")
        ...     print(report.summary())
    """

    def __init__(self, config_path: str, output_dir: Optional[str] = None):
        """
        Initialize coordinator with configuration.

        Args:
            config_path: Path to configuration JSON file
            output_dir: Optional output directory (defaults to config's output_dir)
        """
        self.logger = get_logger()
        # Clear previous issues to prevent contamination across runs
        self.logger.clear_issues()

        self.config = ConfigManager.load(config_path)
        self.output_dir = Path(output_dir) if output_dir else Path(self.config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.validated_files: List[str] = []
        self.failed_files: List[Dict[str, str]] = []
        self.statistics: Dict[str, Any] = {}

        self.logger.info(f"ValidationCoordinator initialized with config: {config_path}")
        self.logger.info(f"Output directory: {self.output_dir}")

    def validate_all(self) -> ValidationReport:
        """
        Run all validations and return comprehensive report.

        Validates in order:
        1. Reference genome
        2. Modified genome
        3. Features (if provided)
        4. All read files

        Returns:
            ValidationReport with complete results
        """
        self.logger.info("Starting validation workflow")

        try:
            # 1. Validate reference genome
            if self.config.ref_genome:
                self._validate_genome(self.config.ref_genome, "Reference Genome")

            # 2. Validate modified genome
            if self.config.mod_genome:
                self._validate_genome(self.config.mod_genome, "Modified Genome")

            # 3. Validate features
            if self.config.feature:
                self._validate_features(self.config.feature)

            # 4. Validate all reads
            if self.config.reads:
                self._validate_reads(self.config.reads)

            # Generate report
            report = self._generate_report()

            self.logger.info(f"Validation workflow completed: {'PASSED' if report.passed else 'FAILED'}")
            return report

        except Exception as e:
            self.logger.error(f"Validation workflow failed: {e}")
            raise

    def validate_genomes(self) -> ValidationReport:
        """Validate only genome files."""
        self.logger.info("Validating genome files only")

        if self.config.ref_genome:
            self._validate_genome(self.config.ref_genome, "Reference Genome")

        if self.config.mod_genome:
            self._validate_genome(self.config.mod_genome, "Modified Genome")

        return self._generate_report()

    def validate_reads(self) -> ValidationReport:
        """Validate only read files."""
        self.logger.info("Validating read files only")

        if self.config.reads:
            self._validate_reads(self.config.reads)

        return self._generate_report()

    def validate_features(self) -> ValidationReport:
        """Validate only feature file."""
        self.logger.info("Validating feature file only")

        if self.config.feature:
            self._validate_features(self.config.feature)

        return self._generate_report()

    def _validate_genome(self, genome_config, label: str) -> None:
        """
        Validate a single genome file.

        Args:
            genome_config: GenomeConfig object
            label: Human-readable label for logging
        """
        try:
            self.logger.info(f"Validating {label}: {genome_config.filename}")

            # Create settings from config if present
            settings = self._create_genome_settings(genome_config)

            # Run validation
            validator = GenomeValidator(genome_config, self.output_dir, settings)
            validator.validate()

            # Collect results
            self.validated_files.append(genome_config.filename)
            self.statistics[genome_config.filename] = validator.get_statistics()

            self.logger.info(f"✓ {label} validated successfully")

        except ValidationError as e:
            self.failed_files.append({
                'filename': genome_config.filename,
                'error': str(e)
            })
            self.logger.error(f"✗ {label} validation failed: {e}")

    def _validate_reads(self, read_configs) -> None:
        """
        Validate all read files.

        Args:
            read_configs: List of ReadConfig objects
        """
        for idx, read_config in enumerate(read_configs, 1):
            try:
                self.logger.info(f"Validating Read file {idx}/{len(read_configs)}: {read_config.filename}")

                # Create settings from config if present
                settings = self._create_read_settings(read_config)

                # Run validation
                validator = ReadValidator(read_config, self.output_dir, settings)
                validator.validate()

                # Collect results
                self.validated_files.append(read_config.filename)
                self.statistics[read_config.filename] = validator.get_statistics()

                self.logger.info(f"✓ Read file {idx} validated successfully")

            except ValidationError as e:
                self.failed_files.append({
                    'filename': read_config.filename,
                    'error': str(e)
                })
                self.logger.error(f"✗ Read file {idx} validation failed: {e}")

    def _validate_features(self, feature_config) -> None:
        """
        Validate feature annotation file.

        Args:
            feature_config: FeatureConfig object
        """
        try:
            self.logger.info(f"Validating Feature file: {feature_config.filename}")

            # Create settings from config if present
            settings = self._create_feature_settings(feature_config)

            # Run validation
            validator = FeatureValidator(feature_config, self.output_dir, settings)
            validator.validate()

            # Collect results
            self.validated_files.append(feature_config.filename)
            self.statistics[feature_config.filename] = validator.get_statistics()

            self.logger.info(f"✓ Feature file validated successfully")

        except ValidationError as e:
            self.failed_files.append({
                'filename': feature_config.filename,
                'error': str(e)
            })
            self.logger.error(f"✗ Feature file validation failed: {e}")

    def _create_genome_settings(self, genome_config) -> Optional[GenomeValidator.Settings]:
        """Create GenomeValidator settings from config extras."""
        # Check if config has extra settings
        extras = getattr(genome_config, 'extras', {})
        if not extras:
            return None

        return GenomeValidator.Settings(**extras)

    def _create_read_settings(self, read_config) -> Optional[ReadValidator.Settings]:
        """Create ReadValidator settings from config extras."""
        # Check if config has extra settings
        extras = getattr(read_config, 'extras', {})
        if not extras:
            return None

        return ReadValidator.Settings(**extras)

    def _create_feature_settings(self, feature_config) -> Optional[FeatureValidator.Settings]:
        """Create FeatureValidator settings from config extras."""
        # Check if config has extra settings
        extras = getattr(feature_config, 'extras', {})
        if not extras:
            return None

        return FeatureValidator.Settings(**extras)

    def _generate_report(self) -> ValidationReport:
        """Generate validation report from collected results."""
        # Count issues from logger
        issues = self.logger.validation_issues
        errors = sum(1 for issue in issues if issue.get('level') == 'ERROR')
        warnings = sum(1 for issue in issues if issue.get('level') in ('WARNING', 'WARN'))

        passed = len(self.failed_files) == 0 and errors == 0

        return ValidationReport(
            passed=passed,
            errors=errors + len(self.failed_files),
            warnings=warnings,
            validated_files=self.validated_files,
            failed_files=self.failed_files,
            statistics=self.statistics,
            output_dir=self.output_dir
        )

    @staticmethod
    def load(config_path: str):
        """
        Load configuration and return coordinator instance.

        This is a convenience method that matches the ConfigManager.load() pattern.

        Args:
            config_path: Path to configuration JSON file

        Returns:
            ValidationCoordinator instance

        Example:
            >>> coordinator = ValidationCoordinator.load("config.json")
            >>> report = coordinator.validate_all()
        """
        return ValidationCoordinator(config_path)


# Backward compatibility alias
Coordinator = ValidationCoordinator
