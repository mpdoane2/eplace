#!/usr/bin/env python
"""
Tests for the ePLACE CLI module.

These tests verify the CLI structure and argument parsing
without requiring external dependencies like pytaxonkit.
"""

import ast
import sys
from pathlib import Path


def test_cli_module_exists():
    """Test that the CLI module file exists."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    assert cli_path.exists(), f"CLI module not found at {cli_path}"


def test_cli_syntax():
    """Test that the CLI module has valid Python syntax."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()
    
    try:
        ast.parse(source)
    except SyntaxError as e:
        raise AssertionError(f"CLI module has syntax error: {e}")


def test_cli_has_main_function():
    """Test that the CLI module has a main() function."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()
    
    tree = ast.parse(source)
    main_found = False
    
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == 'main':
            main_found = True
            break
    
    assert main_found, "main() function not found in CLI module"


def test_cli_has_command_handlers():
    """Test that the CLI module has all required command handler functions."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()
    
    tree = ast.parse(source)
    expected_handlers = {
        'download_command',
        'blast_command', 
        'grouped_command',
        'relabel_command'
    }
    found_handlers = set()
    
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name in expected_handlers:
            found_handlers.add(node.name)
    
    missing = expected_handlers - found_handlers
    assert not missing, f"Missing command handlers: {missing}"


def test_entry_point_registered():
    """Test that the entry point is registered in pyproject.toml."""
    pyproject_path = Path(__file__).parent.parent / "pyproject.toml"
    with open(pyproject_path, 'r') as f:
        content = f.read()
    
    assert '[project.scripts]' in content, "No [project.scripts] section in pyproject.toml"
    assert 'eplace = "eplace_lib.cli:main"' in content, "Entry point not correctly defined"


def test_cli_module_has_docstring():
    """Test that the CLI module has a docstring."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()
    
    tree = ast.parse(source)
    
    # Check module docstring
    docstring = ast.get_docstring(tree)
    assert docstring is not None, "CLI module missing docstring"
    assert "ePLACE" in docstring, "Docstring should mention ePLACE"


def test_cli_imports():
    """Test that the CLI module imports are present."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()
    
    # Check for essential imports
    assert 'import argparse' in source, "Missing argparse import"
    assert 'import logging' in source, "Missing logging import"
    assert 'from pathlib import Path' in source, "Missing Path import"


if __name__ == '__main__':
    # Run tests manually without pytest
    import traceback
    
    tests = [
        test_cli_module_exists,
        test_cli_syntax,
        test_cli_has_main_function,
        test_cli_has_command_handlers,
        test_entry_point_registered,
        test_cli_module_has_docstring,
        test_cli_imports,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            print(f"✓ {test.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"✗ {test.__name__}: {e}")
            failed += 1
        except Exception as e:
            print(f"✗ {test.__name__}: {e}")
            traceback.print_exc()
            failed += 1
    
    print(f"\n{passed} passed, {failed} failed")
    sys.exit(0 if failed == 0 else 1)
