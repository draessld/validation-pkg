"""
Logging configuration for the bioinformatics validation package.

Provides structured logging with multiple outputs:
- Console output (colored, user-friendly)
- File output (detailed, for debugging)
- Validation report (summary of all issues found)
"""

import sys
from pathlib import Path
from datetime import datetime
from typing import Optional
import structlog


# ANSI color codes for console output
class Colors:
    """ANSI color codes for terminal output."""
    RESET = '\033[0m'
    BOLD = '\033[1m'
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    GRAY = '\033[90m'


def add_log_level_colors(_, level: str, event_dict: dict) -> dict:
    """Add colors to log level in console output."""
    level_colors = {
        "debug": Colors.GRAY,
        "info": Colors.BLUE,
        "warning": Colors.YELLOW,
        "error": Colors.RED,
        "critical": Colors.RED + Colors.BOLD,
    }

    level_lower = level.lower()
    color = level_colors.get(level_lower, "")
    if color:
        event_dict["level"] = f"{color}{level.upper()}{Colors.RESET}"
    else:
        event_dict["level"] = level.upper()

    return event_dict


class ValidationLogger:
    """
    Logger for validation package.

    Features:
    - Console output (colored, user-friendly)
    - File output (detailed debug log)
    - Validation report (summary of issues)
    """

    _instance = None
    _initialized = False

    def __new__(cls):
        """Singleton pattern - only one logger instance."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        """Initialize logger (only once)."""
        if self._initialized:
            return

        # Storage for validation issues
        self.validation_issues = []

        # Will hold the structlog logger instance
        self.logger = None

        # Store report file path
        self.report_file = None

        self._initialized = True

    def setup(
        self,
        console_level: str = "INFO",
        log_file: Optional[Path] = None,
        report_file: Optional[Path] = None,
        clear_previous_issues: bool = True
    ):
        """
        Set up logging handlers.

        Args:
            console_level: Level for console output (DEBUG, INFO, WARNING, ERROR)
            log_file: Path to detailed log file (optional)
            report_file: Path to validation report file (optional)
            clear_previous_issues: Clear validation issues from previous runs (default: True)

        Note:
            Setting clear_previous_issues=True helps prevent test contamination and
            ensures each validation run starts with a clean slate.
        """
        # Clear validation issues from previous runs (prevents contamination)
        if clear_previous_issues:
            self.clear_issues()

        # Store report file path
        self.report_file = report_file
        if report_file:
            self.report_file = Path(report_file)
            self.report_file.parent.mkdir(parents=True, exist_ok=True)

        # File logging setup (if requested)
        import logging

        # Setup stdlib logging first (if file logging is requested)
        if log_file:
            log_file = Path(log_file)
            log_file.parent.mkdir(parents=True, exist_ok=True)

            # Use stdlib logging with structlog
            processors = [
                structlog.contextvars.merge_contextvars,
                structlog.stdlib.add_log_level,
                structlog.stdlib.add_logger_name,
                structlog.processors.TimeStamper(fmt="iso"),
                structlog.processors.StackInfoRenderer(),
                structlog.processors.format_exc_info,
                structlog.stdlib.ProcessorFormatter.wrap_for_formatter,
            ]

            # Configure structlog to use stdlib logging
            structlog.configure(
                processors=processors,
                wrapper_class=structlog.stdlib.BoundLogger,
                context_class=dict,
                logger_factory=structlog.stdlib.LoggerFactory(),
                cache_logger_on_first_use=False,
            )

            # Setup stdlib logger for file and console output
            stdlib_logger = logging.getLogger("bioinformatics_validator")
            stdlib_logger.handlers.clear()  # Clear any existing handlers
            stdlib_logger.setLevel(logging.DEBUG)
            stdlib_logger.propagate = False

            # Configure file handler with ProcessorFormatter for structured logs
            file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(
                structlog.stdlib.ProcessorFormatter(
                    processor=structlog.processors.JSONRenderer(),
                    foreign_pre_chain=processors,
                )
            )
            stdlib_logger.addHandler(file_handler)

            # Configure the ProcessorFormatter for console with colors
            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.setLevel(getattr(logging, console_level.upper()))
            console_handler.setFormatter(
                structlog.stdlib.ProcessorFormatter(
                    processor=structlog.dev.ConsoleRenderer(colors=True),
                    foreign_pre_chain=processors,
                )
            )
            stdlib_logger.addHandler(console_handler)

        else:
            # No file logging - use PrintLogger for simplicity
            processors = [
                structlog.contextvars.merge_contextvars,
                structlog.processors.add_log_level,
                structlog.processors.TimeStamper(fmt="iso"),
                structlog.processors.StackInfoRenderer(),
                structlog.processors.format_exc_info,
            ]

            # Console renderer with colors
            console_processors = processors + [
                add_log_level_colors,
                structlog.dev.ConsoleRenderer(colors=True)
            ]

            # Configure structlog
            structlog.configure(
                processors=console_processors,
                wrapper_class=structlog.make_filtering_bound_logger(
                    getattr(structlog.stdlib.logging, console_level.upper(), structlog.stdlib.logging.INFO)
                ),
                context_class=dict,
                logger_factory=structlog.PrintLoggerFactory(file=sys.stdout),
                cache_logger_on_first_use=False,
            )

        # Get logger instance
        self.logger = structlog.get_logger("bioinformatics_validator")

        if log_file:
            self.info(f"Detailed log file: {log_file}")

        if report_file:
            self.info(f"Validation report will be saved to: {report_file}")

    def debug(self, message: str, **kwargs):
        """Log debug message with optional structured context."""
        if self.logger:
            self.logger.debug(message, **kwargs)

    def info(self, message: str, **kwargs):
        """Log info message with optional structured context."""
        if self.logger:
            self.logger.info(message, **kwargs)

    def warning(self, message: str, **kwargs):
        """Log warning message with optional structured context."""
        if self.logger:
            self.logger.warning(message, **kwargs)
        self.validation_issues.append(('WARNING', message))

    def error(self, message: str, **kwargs):
        """Log error message with optional structured context."""
        if self.logger:
            self.logger.error(message, **kwargs)
        self.validation_issues.append(('ERROR', message))

    def critical(self, message: str, **kwargs):
        """Log critical message with optional structured context."""
        if self.logger:
            self.logger.critical(message, **kwargs)
        self.validation_issues.append(('CRITICAL', message))

    def add_validation_issue(self, level: str, category: str, message: str, details: dict = None):
        """
        Add a structured validation issue.

        Args:
            level: ERROR, WARNING, INFO
            category: genome, feature, read, inter-file, etc.
            message: Human-readable description
            details: Additional context (file, line, position, etc.)
        """
        issue = {
            'timestamp': datetime.now().isoformat(),
            'level': level,
            'category': category,
            'message': message,
            'details': details or {}
        }
        self.validation_issues.append(issue)

        # Also log to console (but don't add to validation_issues again!)
        log_message = f"[{category}] {message}"

        # Log with structured context
        log_kwargs = {'category': category}
        if details:
            log_kwargs.update(details)

        # Use logger directly to avoid duplicate appending
        if self.logger:
            if level == 'ERROR':
                self.logger.error(log_message, **log_kwargs)
            elif level == 'WARNING':
                self.logger.warning(log_message, **log_kwargs)
            else:
                self.logger.info(log_message, **log_kwargs)

    def generate_report(self):
        """
        Generate validation report and save to file.

        Returns:
            str: Report content
        """
        if not self.report_file:
            return None

        # Count issues by level
        errors = sum(1 for issue in self.validation_issues
                    if isinstance(issue, dict) and issue.get('level') == 'ERROR')
        warnings = sum(1 for issue in self.validation_issues
                      if isinstance(issue, dict) and issue.get('level') == 'WARNING')

        # Generate report
        report_lines = [
            "=" * 80,
            "BIOINFORMATICS VALIDATION REPORT",
            "=" * 80,
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "SUMMARY",
            "-" * 80,
            f"Total Issues: {len(self.validation_issues)}",
            f"  Errors:   {errors}",
            f"  Warnings: {warnings}",
            "",
            "DETAILS",
            "-" * 80,
        ]

        if not self.validation_issues:
            report_lines.append("✓ No issues found - all validations passed!")
        else:
            for idx, issue in enumerate(self.validation_issues, 1):
                if isinstance(issue, dict):
                    report_lines.append(f"\n{idx}. [{issue['level']}] {issue['category']}")
                    report_lines.append(f"   {issue['message']}")
                    if issue.get('details'):
                        for key, value in issue['details'].items():
                            report_lines.append(f"   - {key}: {value}")
                else:
                    # Old format (tuple)
                    level, message = issue
                    report_lines.append(f"\n{idx}. [{level}] {message}")

        report_lines.append("\n" + "=" * 80)

        report_content = "\n".join(report_lines)

        # Write to file
        self.report_file.write_text(report_content, encoding='utf-8')
        self.info(f"Validation report saved to: {self.report_file}")

        return report_content

    def get_summary(self) -> dict:
        """
        Get summary of validation results.

        Returns:
            dict: Summary statistics
        """
        errors = sum(1 for issue in self.validation_issues
                    if isinstance(issue, dict) and issue.get('level') == 'ERROR')
        warnings = sum(1 for issue in self.validation_issues
                      if isinstance(issue, dict) and issue.get('level') == 'WARNING')

        return {
            'total_issues': len(self.validation_issues),
            'errors': errors,
            'warnings': warnings,
            'passed': errors == 0
        }

    def clear_issues(self):
        """
        Clear all validation issues.

        This is useful for:
        - Starting a new validation run with clean state
        - Test isolation (preventing contamination between tests)
        - Resetting after generating a report

        Example:
            >>> logger = get_logger()
            >>> logger.clear_issues()
            >>> # Now validation_issues list is empty
        """
        self.validation_issues.clear()

    def __enter__(self):
        """
        Context manager entry - clear issues at start.

        This enables using the logger with 'with' statement for automatic isolation:

        Example:
            >>> with get_logger() as logger:
            ...     logger.info("Starting validation")
            ...     # validation_issues starts empty
            ...     # and will be available throughout this block
        """
        self.clear_issues()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Context manager exit - optionally clear issues at end.

        Args:
            exc_type: Exception type if an exception occurred
            exc_val: Exception value
            exc_tb: Exception traceback

        Returns:
            False to propagate exceptions
        """
        # Don't clear issues on exit - caller may want to inspect them
        # If you want to clear, call clear_issues() explicitly
        return False


# Convenience functions for easy import
def get_logger() -> ValidationLogger:
    """Get the singleton logger instance."""
    return ValidationLogger()


def setup_logging(
    console_level: str = "INFO",
    log_file: Optional[Path] = None,
    report_file: Optional[Path] = None
):
    """
    Set up logging for the package.

    Args:
        console_level: Console output level (DEBUG, INFO, WARNING, ERROR)
        log_file: Path to detailed log file
        report_file: Path to validation report
    """
    logger = get_logger()
    logger.setup(console_level, log_file, report_file)
    return logger


# Example usage
if __name__ == "__main__":
    # Example 1: Simple setup
    logger = setup_logging(console_level="INFO")

    logger.info("Starting validation...")
    logger.debug("This is a debug message (won't show in console)")
    logger.warning("This is a warning")
    logger.error("This is an error")

    # Example 2: With files
    logger = setup_logging(
        console_level="INFO",
        log_file=Path("logs/validation.log"),
        report_file=Path("logs/report.txt")
    )

    logger.info("Processing genome file...")
    logger.add_validation_issue(
        level='ERROR',
        category='genome',
        message='Invalid FASTA format',
        details={'file': 'ref.fasta', 'line': 42, 'reason': 'Unexpected character'}
    )

    logger.add_validation_issue(
        level='WARNING',
        category='feature',
        message='Feature outside genome bounds',
        details={'feature_id': 'gene1', 'position': 5000, 'genome_length': 4500}
    )

    # Generate report
    logger.generate_report()

    # Get summary
    summary = logger.get_summary()
    print(f"\nValidation complete: {summary}")
