"""
Basic tests for the logging functionality.
Tests core features without over-testing Python's logging module.
"""

import pytest
import tempfile
from pathlib import Path
from validation_pkg.logger import ValidationLogger, setup_logging, get_logger


class TestLogger:
    """Test suite for ValidationLogger."""
    
    @pytest.fixture(autouse=True)
    def reset_logger(self):
        """Reset logger state before each test."""
        logger = get_logger()
        logger.clear_issues()
        # Reset logger instance (structlog doesn't use handlers in the same way)
        logger.logger = None
        yield
        logger.clear_issues()
    
    def test_logger_singleton(self):
        """Test that logger follows singleton pattern."""
        logger1 = get_logger()
        logger2 = get_logger()
        assert logger1 is logger2
    
    def test_setup_logging_creates_log_file(self):
        """Test that setup_logging creates log file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = Path(tmpdir) / "test.log"
            
            logger = setup_logging(log_file=log_file)
            logger.info("Test message")
            
            assert log_file.exists()
            content = log_file.read_text()
            assert "Test message" in content
    
    def test_add_validation_issue_stores_issue(self):
        """Test that validation issues are stored correctly."""
        logger = get_logger()
        logger.setup()
        
        logger.add_validation_issue(
            level='ERROR',
            category='genome',
            message='Test error message',
            details={'file': 'test.fasta', 'line': 42}
        )
        
        assert len(logger.validation_issues) == 1
        issue = logger.validation_issues[0]
        assert issue['level'] == 'ERROR'
        assert issue['category'] == 'genome'
        assert issue['message'] == 'Test error message'
        assert issue['details']['file'] == 'test.fasta'
    
    def test_add_multiple_validation_issues(self):
        """Test adding multiple issues."""
        logger = get_logger()
        logger.setup()
        
        logger.add_validation_issue('ERROR', 'genome', 'Error 1')
        logger.add_validation_issue('WARNING', 'feature', 'Warning 1')
        logger.add_validation_issue('ERROR', 'read', 'Error 2')
        
        assert len(logger.validation_issues) == 3
    
    def test_get_summary_counts_issues(self):
        """Test that get_summary counts issues correctly."""
        logger = get_logger()
        logger.setup()
        
        logger.add_validation_issue('ERROR', 'genome', 'Error 1')
        logger.add_validation_issue('ERROR', 'genome', 'Error 2')
        logger.add_validation_issue('WARNING', 'feature', 'Warning 1')
        
        summary = logger.get_summary()
        
        assert summary['total_issues'] == 3
        assert summary['errors'] == 2
        assert summary['warnings'] == 1
        assert summary['passed'] == False
    
    def test_get_summary_passed_when_no_errors(self):
        """Test that passed=True when there are no errors."""
        logger = get_logger()
        logger.setup()
        
        logger.add_validation_issue('WARNING', 'feature', 'Warning only')
        
        summary = logger.get_summary()
        
        assert summary['errors'] == 0
        assert summary['passed'] == True
    
    def test_get_summary_empty(self):
        """Test get_summary with no issues."""
        logger = get_logger()
        logger.setup()
        
        summary = logger.get_summary()
        
        assert summary['total_issues'] == 0
        assert summary['errors'] == 0
        assert summary['warnings'] == 0
        assert summary['passed'] == True
    
    def test_generate_report_creates_file(self):
        """Test that generate_report creates report file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_file = Path(tmpdir) / "report.txt"
            
            logger = setup_logging(report_file=report_file)
            logger.add_validation_issue('ERROR', 'genome', 'Test error')
            
            report_content = logger.generate_report()
            
            assert report_file.exists()
            assert "BIOINFORMATICS VALIDATION REPORT" in report_content
            assert "Test error" in report_content
    
    def test_generate_report_contains_all_issues(self):
        """Test that report contains all logged issues."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_file = Path(tmpdir) / "report.txt"
            
            logger = setup_logging(report_file=report_file)
            logger.add_validation_issue('ERROR', 'genome', 'Error message')
            logger.add_validation_issue('WARNING', 'feature', 'Warning message')
            
            report_content = logger.generate_report()
            
            assert "Error message" in report_content
            assert "Warning message" in report_content
            assert "[ERROR]" in report_content
            assert "[WARNING]" in report_content
    
    def test_generate_report_shows_details(self):
        """Test that report includes issue details."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_file = Path(tmpdir) / "report.txt"
            
            logger = setup_logging(report_file=report_file)
            logger.add_validation_issue(
                level='ERROR',
                category='genome',
                message='Test error',
                details={'file': 'test.fasta', 'line': 42}
            )
            
            report_content = logger.generate_report()
            
            assert "test.fasta" in report_content
            assert "42" in report_content
    
    def test_generate_report_no_issues(self):
        """Test report when no issues found."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_file = Path(tmpdir) / "report.txt"
            
            logger = setup_logging(report_file=report_file)
            report_content = logger.generate_report()
            
            assert "No issues found" in report_content
            assert "all validations passed" in report_content.lower()
    
    def test_clear_issues(self):
        """Test that clear_issues removes all issues."""
        logger = get_logger()
        logger.setup()
        
        logger.add_validation_issue('ERROR', 'genome', 'Error 1')
        logger.add_validation_issue('WARNING', 'feature', 'Warning 1')
        
        assert len(logger.validation_issues) == 2
        
        logger.clear_issues()
        
        assert len(logger.validation_issues) == 0
    
    def test_log_file_directory_creation(self):
        """Test that log file parent directories are created."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = Path(tmpdir) / "nested" / "dir" / "test.log"
            
            logger = setup_logging(log_file=log_file)
            logger.info("Test")
            
            assert log_file.exists()
            assert log_file.parent.exists()
    
    def test_report_file_directory_creation(self):
        """Test that report file parent directories are created."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_file = Path(tmpdir) / "nested" / "dir" / "report.txt"
            
            logger = setup_logging(report_file=report_file)
            logger.generate_report()
            
            assert report_file.exists()
            assert report_file.parent.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])