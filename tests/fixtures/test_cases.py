#!/usr/bin/env python3
"""
Utility script for managing test fixtures.

NOTE: For running integration tests, use:
    pytest tests/test_integration.py -v

This script provides utility functions for:
- Listing available test fixtures
- Regenerating fixtures
- Showing statistics

For full documentation, see: tests/README.md

Usage:
    python test_cases.py --list             # List available test cases
    python test_cases.py --regenerate       # Regenerate all fixtures
    python test_cases.py --stats            # Show statistics

To run tests:
    pytest tests/test_integration.py        # Run all integration tests
"""

import sys
import subprocess
import argparse
from pathlib import Path
from typing import List, Dict


class TestCaseRunner:
    """Runner for integration tests using fixtures."""

    def __init__(self, fixtures_dir: Path = None):
        if fixtures_dir is None:
            self.fixtures_dir = Path(__file__).parent
        else:
            self.fixtures_dir = Path(fixtures_dir)

        self.tests_dir = self.fixtures_dir.parent
        self.project_root = self.tests_dir.parent

    def discover_fixtures(self) -> Dict[str, List[Path]]:
        """Discover all test fixtures in the fixtures directory.

        Returns:
            Dictionary with categories as keys and list of fixture paths as values
        """
        fixtures = {
            'config': [],
            'genome': []
        }

        # Discover config test fixtures
        for fixture_dir in sorted(self.fixtures_dir.glob("config_test_*")):
            if fixture_dir.is_dir():
                fixtures['config'].append(fixture_dir)

        # Discover genome test fixtures
        for fixture_dir in sorted(self.fixtures_dir.glob("test_case_*")):
            if fixture_dir.is_dir():
                fixtures['genome'].append(fixture_dir)

        return fixtures

    def list_test_cases(self) -> None:
        """List all available test cases."""
        fixtures = self.discover_fixtures()

        print("\n" + "=" * 70)
        print("AVAILABLE TEST FIXTURES")
        print("=" * 70)

        print(f"\nConfig Manager Test Cases ({len(fixtures['config'])} total):")
        print("-" * 70)
        for fixture in fixtures['config']:
            desc_file = fixture / "description.txt"
            if desc_file.exists():
                desc = desc_file.read_text().split('\n')[0]
                print(f"  • {fixture.name}")
                print(f"    {desc}")
            else:
                print(f"  • {fixture.name}")

        print(f"\nGenome Validator Test Cases ({len(fixtures['genome'])} total):")
        print("-" * 70)
        for fixture in fixtures['genome']:
            desc_file = fixture / "description.txt"
            if desc_file.exists():
                desc = desc_file.read_text().split('\n')[0]
                print(f"  • {fixture.name}")
                print(f"    {desc}")
            else:
                print(f"  • {fixture.name}")

        print("\n" + "=" * 70)
        print(f"Total: {len(fixtures['config']) + len(fixtures['genome'])} test fixtures")
        print("=" * 70 + "\n")

    def regenerate_fixtures(self) -> bool:
        """Regenerate all test fixtures.

        Returns:
            True if successful, False otherwise
        """
        print("\n" + "=" * 70)
        print("REGENERATING TEST FIXTURES")
        print("=" * 70 + "\n")

        scripts = [
            self.fixtures_dir / "generate_config_fixtures.sh",
            self.fixtures_dir / "generate_fixtures.sh"
        ]

        for script in scripts:
            if not script.exists():
                print(f"✗ Script not found: {script}")
                continue

            print(f"Running {script.name}...")
            try:
                result = subprocess.run(
                    [str(script)],
                    cwd=self.fixtures_dir,
                    capture_output=True,
                    text=True,
                    check=True
                )
                print(result.stdout)
                if result.stderr:
                    print("Warnings:", result.stderr)
            except subprocess.CalledProcessError as e:
                print(f"✗ Failed to run {script.name}")
                print(f"Error: {e.stderr}")
                return False

        print("\n✓ All fixtures regenerated successfully!\n")
        return True

    def run_tests(self, config_only: bool = False, genome_only: bool = False,
                  verbose: bool = False, pytest_args: List[str] = None) -> int:
        """Run integration tests.

        Args:
            config_only: Run only config manager tests
            genome_only: Run only genome validator tests
            verbose: Verbose pytest output
            pytest_args: Additional pytest arguments

        Returns:
            Exit code from pytest
        """
        if pytest_args is None:
            pytest_args = []

        # Build pytest command
        cmd = ["python", "-m", "pytest"]

        # Determine which tests to run
        if config_only:
            test_files = [str(self.tests_dir / "test_config_manager_integration.py")]
            print("\n" + "=" * 70)
            print("RUNNING CONFIG MANAGER INTEGRATION TESTS")
            print("=" * 70 + "\n")
        elif genome_only:
            test_files = [str(self.tests_dir / "test_genome_validator_integration.py")]
            print("\n" + "=" * 70)
            print("RUNNING GENOME VALIDATOR INTEGRATION TESTS")
            print("=" * 70 + "\n")
        else:
            test_files = [
                str(self.tests_dir / "test_config_manager_integration.py"),
                str(self.tests_dir / "test_genome_validator_integration.py")
            ]
            print("\n" + "=" * 70)
            print("RUNNING ALL INTEGRATION TESTS")
            print("=" * 70 + "\n")

        cmd.extend(test_files)

        # Add verbosity
        if verbose:
            cmd.append("-vv")
        else:
            cmd.append("-v")

        # Add any additional pytest args
        cmd.extend(pytest_args)

        # Run tests
        try:
            result = subprocess.run(cmd, cwd=self.project_root)
            return result.returncode
        except Exception as e:
            print(f"\n✗ Error running tests: {e}")
            return 1

    def get_statistics(self) -> Dict[str, int]:
        """Get statistics about test fixtures.

        Returns:
            Dictionary with statistics
        """
        fixtures = self.discover_fixtures()

        stats = {
            'config_tests': len(fixtures['config']),
            'genome_tests': len(fixtures['genome']),
            'total_fixtures': len(fixtures['config']) + len(fixtures['genome'])
        }

        return stats


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Run integration tests using test fixtures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python test_cases.py                # Run all integration tests
  python test_cases.py --config       # Run only config manager tests
  python test_cases.py --genome       # Run only genome validator tests
  python test_cases.py --list         # List all test cases
  python test_cases.py --regenerate   # Regenerate all fixtures
  python test_cases.py --verbose      # Verbose output
  python test_cases.py -k test_case_01  # Run specific test pattern
        """
    )

    parser.add_argument(
        '--config',
        action='store_true',
        help='Run only config manager integration tests'
    )

    parser.add_argument(
        '--genome',
        action='store_true',
        help='Run only genome validator integration tests'
    )

    parser.add_argument(
        '--list',
        action='store_true',
        help='List all available test cases'
    )

    parser.add_argument(
        '--regenerate',
        action='store_true',
        help='Regenerate all test fixtures'
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output'
    )

    parser.add_argument(
        '--stats',
        action='store_true',
        help='Show test fixture statistics'
    )

    # Allow passing through pytest args
    parser.add_argument(
        'pytest_args',
        nargs='*',
        help='Additional arguments to pass to pytest (e.g., -k pattern, -x)'
    )

    args = parser.parse_args()

    # Initialize runner
    runner = TestCaseRunner()

    # Handle list command
    if args.list:
        runner.list_test_cases()
        return 0

    # Handle regenerate command
    if args.regenerate:
        success = runner.regenerate_fixtures()
        return 0 if success else 1

    # Handle stats command
    if args.stats:
        stats = runner.get_statistics()
        print("\n" + "=" * 70)
        print("TEST FIXTURE STATISTICS")
        print("=" * 70)
        print(f"  Config Manager Tests:    {stats['config_tests']}")
        print(f"  Genome Validator Tests:  {stats['genome_tests']}")
        print(f"  Total Fixtures:          {stats['total_fixtures']}")
        print("=" * 70 + "\n")
        return 0

    # Run tests
    exit_code = runner.run_tests(
        config_only=args.config,
        genome_only=args.genome,
        verbose=args.verbose,
        pytest_args=args.pytest_args
    )

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
