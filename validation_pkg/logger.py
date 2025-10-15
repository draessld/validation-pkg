"""
Logging configuration for the bioinformatics validation package.

Provides structured logging with multiple outputs:
- Console output (colored, user-friendly)
- File output (detailed, for debugging)
- Validation report (summary of all issues found)
"""

import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional


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


class ColoredFormatter(logging.Formatter):
    """Custom formatter with colors for console output."""
    
    FORMATS = {
        logging.DEBUG: Colors.GRAY + '%(levelname)s' + Colors.RESET + ' - %(message)s',
        logging.INFO: Colors.BLUE + '%(levelname)s' + Colors.RESET + ' - %(message)s',
        logging.WARNING: Colors.YELLOW + '%(levelname)s' + Colors.RESET + ' - %(message)s',
        logging.ERROR: Colors.RED + '%(levelname)s' + Colors.RESET + ' - %(message)s',
        logging.CRITICAL: Colors.RED + Colors.BOLD + '%(levelname)s' + Colors.RESET + ' - %(message)s',
    }
    
    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


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
        
        self.logger = logging.getLogger('bioinformatics_validator')
        self.logger.setLevel(logging.DEBUG)
        self.logger.propagate = False
        
        # Storage for validation issues
        self.validation_issues = []
        
        self._initialized = True
    
    def setup(
        self,
        console_level: str = "INFO",
        log_file: Optional[Path] = None,
        report_file: Optional[Path] = None
    ):
        """
        Set up logging handlers.
        
        Args:
            console_level: Level for console output (DEBUG, INFO, WARNING, ERROR)
            log_file: Path to detailed log file (optional)
            report_file: Path to validation report file (optional)
        """
        # Clear existing handlers
        self.logger.handlers.clear()
        
        # Console handler (colored, user-friendly)
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(getattr(logging, console_level.upper()))
        console_handler.setFormatter(ColoredFormatter())
        self.logger.addHandler(console_handler)
        
        # File handler (detailed debug log)
        if log_file:
            log_file = Path(log_file)
            log_file.parent.mkdir(parents=True, exist_ok=True)
            
            file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
            file_handler.setFormatter(file_formatter)
            self.logger.addHandler(file_handler)
            
            self.info(f"Detailed log file: {log_file}")
        
        # Store report file path
        self.report_file = report_file
        if report_file:
            self.report_file = Path(report_file)
            self.report_file.parent.mkdir(parents=True, exist_ok=True)
            self.info(f"Validation report will be saved to: {report_file}")
    
    def debug(self, message: str):
        """Log debug message."""
        self.logger.debug(message)
    
    def info(self, message: str):
        """Log info message."""
        self.logger.info(message)
    
    def warning(self, message: str):
        """Log warning message."""
        self.logger.warning(message)
        self.validation_issues.append(('WARNING', message))
    
    def error(self, message: str):
        """Log error message."""
        self.logger.error(message)
        self.validation_issues.append(('ERROR', message))
    
    def critical(self, message: str):
        """Log critical message."""
        self.logger.critical(message)
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
        
        # Also log to console
        log_message = f"[{category}] {message}"
        if details:
            log_message += f" | Details: {details}"
        
        if level == 'ERROR':
            self.error(log_message)
        elif level == 'WARNING':
            self.warning(log_message)
        else:
            self.info(log_message)
    
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
        """Clear all validation issues (useful for testing)."""
        self.validation_issues.clear()


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