#!/usr/bin/env python
"""
Tests for the ePLACE relabel command.

These tests verify that the relabel command is properly integrated
into the CLI and has the correct structure.
"""

import ast
import sys
from pathlib import Path


def test_relabel_command_exists():
    """Test that the relabel_command function exists."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()
    
    tree = ast.parse(source)
    relabel_command_found = False
    
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == 'relabel_command':
            relabel_command_found = True
            # Check that it has the expected signature (args parameter)
            assert len(node.args.args) > 0, "relabel_command should have at least one parameter"
            break
    
    assert relabel_command_found, "relabel_command function not found in CLI module"


def test_relabel_parser_exists():
    """Test that the relabel subcommand parser is defined."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()
    
    assert 'relabel_parser' in source, "relabel_parser not found in CLI module"
    assert "subparsers.add_parser(\n        'relabel'" in source or \
           'subparsers.add_parser("relabel"' in source or \
           "subparsers.add_parser('relabel'" in source, \
           "relabel subparser not properly added"


def test_relabel_command_routing():
    """Test that the relabel command is routed correctly."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()
    
    assert "args.command == 'relabel'" in source, \
           "relabel command routing not found"
    assert "relabel_command(args)" in source, \
           "relabel_command call not found in routing"


def test_relabel_required_arguments():
    """Test that the relabel parser has the required arguments."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()
    
    # Find the relabel parser section
    relabel_section_start = source.find("# Relabel subcommand")
    relabel_section_end = source.find("# Parse arguments", relabel_section_start)
    relabel_section = source[relabel_section_start:relabel_section_end]
    
    # Check for required arguments
    assert "'blast_output'" in relabel_section, "blast_output argument not defined"
    assert "'tree_file'" in relabel_section, "tree_file argument not defined"
    assert "'output_tree'" in relabel_section, "output_tree argument not defined"
    assert "'--rank'" in relabel_section, "--rank argument not defined"


def test_relabel_rank_choices():
    """Test that the relabel --rank argument has the correct choices."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()
    
    # Find the relabel parser section
    relabel_section_start = source.find("# Relabel subcommand")
    relabel_section_end = source.find("# Parse arguments", relabel_section_start)
    relabel_section = source[relabel_section_start:relabel_section_end]
    
    # Check that rank choices include all expected taxonomic ranks
    expected_ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    for rank in expected_ranks:
        assert f"'{rank}'" in relabel_section, f"rank choice '{rank}' not found"


def test_relabel_command_imports():
    """Test that relabel_command imports the required modules."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()
    
    # Find the relabel_command function
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == 'relabel_command':
            func_source = ast.get_source_segment(source, node)
            if func_source:
                # Check for required imports within the function
                assert 'BlastRunner' in func_source, "BlastRunner not imported in relabel_command"
                assert 'TaxonomyExtractor' in func_source, "TaxonomyExtractor not imported in relabel_command"
                break


def test_relabel_species_handling():
    """Test that relabel_command has special handling for species rank."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()
    
    # Find the relabel_command function
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == 'relabel_command':
            func_source = ast.get_source_segment(source, node)
            if func_source:
                # Check for species handling logic
                assert "args.rank == 'species'" in func_source or \
                       'args.rank == "species"' in func_source, \
                       "Special handling for species rank not found"
                assert 'genus' in func_source, \
                       "genus not referenced in species handling (expected for compound name)"
                break


def test_relabel_in_main_help():
    """Test that relabel command is mentioned in the main help."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()
    
    # Check that relabel is mentioned in the main help/epilog
    assert 'relabel' in source.lower(), "relabel not mentioned in CLI help"


if __name__ == '__main__':
    # Run tests manually without pytest
    import traceback
    
    tests = [
        test_relabel_command_exists,
        test_relabel_parser_exists,
        test_relabel_command_routing,
        test_relabel_required_arguments,
        test_relabel_rank_choices,
        test_relabel_command_imports,
        test_relabel_species_handling,
        test_relabel_in_main_help,
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
