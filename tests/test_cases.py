"""
Generic fixture-based integration testing.
Tests all cases in fixtures/ directory by running the package with config.json.

This test file is designed to test the ENTIRE package workflow, not just configuration.
It runs your actual validation/processing pipeline on each test case.
"""

import pytest
from pathlib import Path
import sys

# Import your main package entry point
# TODO: Replace with your actual package function
# from your_package import run_validation
# from your_package.coordinator import Coordinator


def run_package_with_config(config_file: Path):
    """
    Run your package with the given config file.
    
    This is where you call your main package function.
    Replace this with your actual implementation.
    
    Args:
        config_file: Path to config.json
        
    Raises:
        Exception: If validation/processing fails
    """
    # Example implementation:
    # from your_package.coordinator import Coordinator
    # from your_package.main import run_validation
    # 
    # config = Coordinator.load(str(config_file))
    # run_validation(config)
    
    # Temporary implementation for testing coordinator only:
    from validation_pkg.coordinator import Coordinator
    config = Coordinator.load(str(config_file))
    # When you have validators, add them here:
    # genome_validator = GenomeValidator(config.ref_genome)
    # genome_validator.validate()
    # ... etc


def test_all_cases():
    """
    Test all cases in fixtures/ directory.
    
    Naming convention:
    - Folders starting with 'valid_' should pass without exceptions
    - Folders starting with 'invalid_' or 'error_' should raise exceptions
    """
    fixtures_dir = Path(__file__).parent / "fixtures"
    
    # Check if fixtures exist
    if not fixtures_dir.exists():
        pytest.fail(
            f"Fixtures directory not found: {fixtures_dir}\n"
            f"Run './create_fixtures.sh' to create test fixtures."
        )
    
    # Get all test case directories
    test_cases = sorted([d for d in fixtures_dir.iterdir() if d.is_dir()])
    
    if not test_cases:
        pytest.fail(
            f"No test cases found in {fixtures_dir}\n"
            f"Run './create_fixtures.sh' to create test fixtures."
        )
    
    passed = []
    failed = []
    
    print(f"\n{'='*70}")
    print(f"RUNNING FIXTURE-BASED INTEGRATION TESTS")
    print(f"{'='*70}")
    print(f"Fixtures directory: {fixtures_dir}")
    print(f"Total test cases: {len(test_cases)}")
    print(f"{'='*70}\n")
    
    for case_dir in test_cases:
        config_file = case_dir / "config.json"
        
        # Skip if no config.json
        if not config_file.exists():
            failed.append(f"{case_dir.name}: Missing config.json")
            print(f"⚠️  {case_dir.name}: SKIP (no config.json)")
            continue
        
        # Determine expected behavior based on folder name
        should_fail = (
            "invalid" in case_dir.name.lower() or 
            "error" in case_dir.name.lower() or
            "bad" in case_dir.name.lower()
        )
        
        # Run the test
        try:
            run_package_with_config(config_file)
            
            # No exception was raised
            if should_fail:
                failed.append(f"{case_dir.name}: Expected exception but none was raised")
                print(f"✗  {case_dir.name}: FAIL (should have raised exception)")
            else:
                passed.append(case_dir.name)
                print(f"✓  {case_dir.name}: PASS")
                
        except Exception as e:
            # Exception was raised
            exception_type = type(e).__name__
            exception_msg = str(e)
            
            if should_fail:
                passed.append(case_dir.name)
                print(f"✓  {case_dir.name}: PASS (raised {exception_type})")
            else:
                failed.append(f"{case_dir.name}: Unexpected exception: {exception_type}: {exception_msg}")
                print(f"✗  {case_dir.name}: FAIL ({exception_type}: {exception_msg})")
    
    # Print summary
    print(f"\n{'='*70}")
    print(f"TEST SUMMARY")
    print(f"{'='*70}")
    print(f"Total:  {len(test_cases)}")
    print(f"Passed: {len(passed)} ({len(passed)*100//len(test_cases) if test_cases else 0}%)")
    print(f"Failed: {len(failed)} ({len(failed)*100//len(test_cases) if test_cases else 0}%)")
    print(f"{'='*70}")
    
    if passed:
        print(f"\n✓ Passed cases ({len(passed)}):")
        for case in passed:
            print(f"  • {case}")
    
    if failed:
        print(f"\n✗ Failed cases ({len(failed)}):")
        for fail in failed:
            print(f"  • {fail}")
        print(f"{'='*70}\n")
    
    # Fail the test if any case failed
    assert len(failed) == 0, (
        f"\n{len(failed)} test case(s) failed!\n"
        f"See details above."
    )


def test_individual_case():
    """
    Example of how to test a single case during development.
    Useful for debugging specific test cases.
    """
    # To test a specific case, uncomment and modify:
    case_dir = Path(__file__).parent / "fixtures" / "valid_minimal"
    config_file = case_dir / "config.json"
    run_package_with_config(config_file)
    pass


if __name__ == "__main__":
    # Allow running this file directly for quick testing
    print("Running fixture-based integration tests...\n")
    
    try:
        test_all_cases()
        print("\n✅ All tests passed!")
        sys.exit(0)
    except AssertionError as e:
        print(f"\n❌ Tests failed: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Error running tests: {e}")
        sys.exit(1)