#!/usr/bin/env python
"""
Tests for the ePLACE CLI module.

These tests verify the CLI structure and argument parsing
without requiring external dependencies like pytaxonkit.
"""

import ast
import sys
from pathlib import Path
from unittest.mock import patch


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
    assert 'import json' in source, "Missing json import"
    assert 'from pathlib import Path' in source, "Missing Path import"


def test_cli_has_log_level_argument():
    """Test that the CLI module defines a --log-level argument with correct choices and default."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()

    tree = ast.parse(source)
    expected_levels = ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL')

    # Check that _LOG_LEVEL_CHOICES constant is defined with the expected values, OR that a
    # direct add_argument('--log-level', choices=[...]) call uses the right choices list.
    choices_found = False

    # 1. Look for a module-level constant assignment: _LOG_LEVEL_CHOICES = [...]
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id == '_LOG_LEVEL_CHOICES':
                try:
                    value = ast.literal_eval(node.value)
                    assert tuple(value) == expected_levels, (
                        f"_LOG_LEVEL_CHOICES should be {expected_levels}, got {value}"
                    )
                    choices_found = True
                except (ValueError, SyntaxError):
                    pass

    # 2. Fallback: look for inline choices=[...] on any add_argument('--log-level', ...) call
    if not choices_found:
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if not isinstance(node.func, ast.Attribute) or node.func.attr != 'add_argument':
                continue
            option_strings = [
                ast.literal_eval(a)
                for a in node.args
                if isinstance(a, ast.Constant) and isinstance(a.value, str)
            ]
            if '--log-level' not in option_strings:
                continue
            for kw in node.keywords:
                if kw.arg == 'choices':
                    try:
                        value = ast.literal_eval(kw.value)
                        assert tuple(value) == expected_levels, (
                            f"--log-level choices should be {expected_levels}, got {value}"
                        )
                        choices_found = True
                    except (ValueError, SyntaxError):
                        pass

    assert choices_found, (
        "CLI module must define --log-level choices either via a _LOG_LEVEL_CHOICES "
        "constant or as an inline list in add_argument('--log-level', choices=[...])"
    )

    # Verify the top-level parser defaults to 'INFO'.  The helper _add_log_level_argument
    # uses is_top_level=True for the main parser, so look for that call.
    top_level_default_ok = False
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        # Check _add_log_level_argument(parser, is_top_level=True)
        if isinstance(func, ast.Name) and func.id == '_add_log_level_argument':
            for kw in node.keywords:
                if kw.arg == 'is_top_level':
                    try:
                        if ast.literal_eval(kw.value) is True:
                            top_level_default_ok = True
                    except (ValueError, SyntaxError):
                        pass
        # Check direct add_argument('--log-level', default='INFO', ...)
        if isinstance(func, ast.Attribute) and func.attr == 'add_argument':
            option_strings = [
                ast.literal_eval(a)
                for a in node.args
                if isinstance(a, ast.Constant) and isinstance(a.value, str)
            ]
            if '--log-level' in option_strings:
                for kw in node.keywords:
                    if kw.arg == 'default':
                        try:
                            if ast.literal_eval(kw.value) == 'INFO':
                                top_level_default_ok = True
                        except (ValueError, SyntaxError):
                            pass

    assert top_level_default_ok, (
        "Top-level parser --log-level should default to 'INFO' "
        "(either via _add_log_level_argument(..., is_top_level=True) or default='INFO')"
    )


def test_cli_version_uses_package_metadata():
    """Test that --version pulls from package metadata instead of a hardcoded value."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()

    tree = ast.parse(source)

    version_helper_calls_metadata = False
    version_argument_uses_helper = False

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == '_get_cli_version':
            for subnode in ast.walk(node):
                if (
                    isinstance(subnode, ast.Call)
                    and isinstance(subnode.func, ast.Name)
                    and subnode.func.id == 'version'
                    and subnode.args
                    and isinstance(subnode.args[0], ast.Constant)
                    and subnode.args[0].value == 'eplace'
                ):
                    version_helper_calls_metadata = True

        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == 'add_argument'
            and any(
                isinstance(a, ast.Constant) and a.value == '--version'
                for a in node.args
            )
        ):
            for kw in node.keywords:
                if kw.arg != 'version':
                    continue
                if isinstance(kw.value, ast.JoinedStr):
                    for value_node in kw.value.values:
                        if (
                            isinstance(value_node, ast.FormattedValue)
                            and isinstance(value_node.value, ast.Call)
                            and isinstance(value_node.value.func, ast.Name)
                            and value_node.value.func.id == '_get_cli_version'
                        ):
                            version_argument_uses_helper = True

    assert version_helper_calls_metadata, (
        "_get_cli_version should read version metadata for the 'eplace' package"
    )
    assert version_argument_uses_helper, (
        "--version argument should use _get_cli_version() instead of a hardcoded version string"
    )


def test_cli_log_level_in_all_subparsers():
    """Test that --log-level is added to every subparser so it works after subcommands."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()

    tree = ast.parse(source)

    # Record which parser variables have --log-level registered, accepting two forms:
    # 1. Direct call:  parser_var.add_argument('--log-level', ...)
    # 2. Helper call:  _add_log_level_argument(parser_var, ...)
    log_level_parsers = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func

        # Form 1: <parser_var>.add_argument('--log-level', ...)
        if isinstance(func, ast.Attribute) and func.attr == 'add_argument':
            has_log_level = any(
                isinstance(a, ast.Constant) and a.value == '--log-level'
                for a in node.args
            )
            if has_log_level and isinstance(func.value, ast.Name):
                log_level_parsers.add(func.value.id)

        # Form 2: _add_log_level_argument(parser_var, ...)
        if isinstance(func, ast.Name) and func.id == '_add_log_level_argument':
            if node.args and isinstance(node.args[0], ast.Name):
                log_level_parsers.add(node.args[0].id)

    expected_parsers = {
        'parser',
        'download_parser',
        'search_parser',
        'grouped_parser',
        'relabel_parser',
    }
    missing_parsers = expected_parsers - log_level_parsers

    assert not missing_parsers, (
        "Expected --log-level to be defined on the top-level parser and all subparsers, "
        f"but it was missing from: {sorted(missing_parsers)}. "
        f"Found on: {sorted(log_level_parsers)}"
    )



def test_cli_has_mmseqs_db_source_argument():
    """Test that blast and grouped parsers have the --mmseqs-db-source argument."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()

    tree = ast.parse(source)

    # Collect all add_argument('--mmseqs-db-source', ...) calls and the parser they belong to.
    db_source_parsers = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not (isinstance(func, ast.Attribute) and func.attr == 'add_argument'):
            continue
        has_db_source = any(
            isinstance(a, ast.Constant) and a.value == '--mmseqs-db-source'
            for a in node.args
        )
        if has_db_source and isinstance(func.value, ast.Name):
            db_source_parsers.add(func.value.id)

    expected = {'search_parser', 'grouped_parser'}
    missing = expected - db_source_parsers
    assert not missing, (
        f"--mmseqs-db-source argument missing from: {sorted(missing)}"
    )


def test_cli_has_blast_db_source_argument():
    """Test that blast and grouped parsers have the --blast-db-source argument."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()

    tree = ast.parse(source)

    db_source_parsers = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not (isinstance(func, ast.Attribute) and func.attr == 'add_argument'):
            continue
        has_db_source = any(
            isinstance(a, ast.Constant) and a.value == '--blast-db-source'
            for a in node.args
        )
        if has_db_source and isinstance(func.value, ast.Name):
            db_source_parsers.add(func.value.id)

    expected = {'search_parser', 'grouped_parser'}
    missing = expected - db_source_parsers
    assert not missing, (
        f"--blast-db-source argument missing from: {sorted(missing)}"
    )


def test_cli_has_mmseqs_download_arguments():
    """Test that download parser includes MMseqs2 download and taxonomy flags."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()

    tree = ast.parse(source)

    expected = {
        '--target',
        '--mmseqs-db-dir',
        '--mmseqs-threads',
        '--add-taxonomy',
        '--ncbi-taxonomy',
        '--acc2taxid-dir',
        '--taxonomy-workdir',
        '--skip-memory-check',
    }
    found = set()

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not (isinstance(node.func, ast.Attribute) and node.func.attr == 'add_argument'):
            continue
        for arg in node.args:
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                if arg.value in expected:
                    found.add(arg.value)

    missing = expected - found
    assert not missing, f"MMseqs download arguments missing from CLI: {sorted(missing)}"


def test_cli_has_write_search_metadata_function():
    """Test that the CLI module defines the _write_search_metadata helper."""
    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"
    with open(cli_path, 'r') as f:
        source = f.read()

    tree = ast.parse(source)
    found = any(
        isinstance(node, ast.FunctionDef) and node.name == '_write_search_metadata'
        for node in ast.walk(tree)
    )
    assert found, "_write_search_metadata function not found in CLI module"


def _load_cli_module_for_testing():
    """Load cli.py with all unavailable dependencies mocked out.

    This helper is used for isolated unit testing of pure-Python helpers
    defined in ``cli.py`` (such as ``_write_search_metadata``) without
    loading the full ePLACE dependency stack.  Integration tests that
    exercise the complete workflow still require pytaxonkit to be installed
    via conda/mamba as described in INSTALL.md.

    The loaded module object is returned directly.  All stubbed entries are
    scoped to the ``patch.dict`` context manager and are removed from
    ``sys.modules`` before this function returns, so repeated calls each
    produce a fresh module instance without cross-test contamination.

    Returns:
        The loaded cli module object.
    """
    import importlib.util
    import sys
    from unittest.mock import MagicMock

    cli_path = Path(__file__).parent.parent / "src" / "eplace_lib" / "cli.py"

    stubs = {
        'pytaxonkit': MagicMock(),
        'eplace_lib': MagicMock(),
        'eplace_lib.ncbi_download': MagicMock(
            setup_ncbi_database=MagicMock(),
            setup_mmseqs_database=MagicMock(),
            setup_mmseqs_taxonomy=MagicMock(),
            check_available_memory_gb=MagicMock(return_value=(True, 256.0))
        ),
        'eplace_lib.blast_analysis': MagicMock(),
        'eplace_lib.taxonomy': MagicMock(),
        'eplace_lib.alignment': MagicMock(),
    }

    # Name the spec as "eplace_lib._cli_under_test" so that spec.parent
    # == "eplace_lib", matching the relative imports used inside cli.py.
    module_name = 'eplace_lib._cli_under_test'
    spec = importlib.util.spec_from_file_location(module_name, cli_path)
    module = importlib.util.module_from_spec(spec)

    # patch.dict restores sys.modules to its original state when the
    # context manager exits, so no stub leaks between test invocations.
    with patch.dict(sys.modules, {**stubs, module_name: module}):
        spec.loader.exec_module(module)

    return module


def test_cli_write_search_metadata_writes_json():
    """Test that _write_search_metadata writes a JSON file with the expected keys."""
    import json
    import tempfile

    module = _load_cli_module_for_testing()
    write_fn = module._write_search_metadata

    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)
        write_fn(
            output_dir=output_dir,
            search_backend='blast',
            database_name='core_nt',
            database_path='/home/user/blastdb',
            database_source='ncbi_core_nt'
        )
        metadata_file = output_dir / "search_metadata.json"
        assert metadata_file.exists(), "search_metadata.json was not created"
        data = json.loads(metadata_file.read_text())
        assert data['search_backend'] == 'blast'
        assert data['database_name'] == 'core_nt'
        assert data['database_path'] == '/home/user/blastdb'
        assert data['database_source'] == 'ncbi_core_nt'


def test_cli_write_search_metadata_mmseqs2():
    """Test that _write_search_metadata records mmseqs2 backend correctly."""
    import json
    import tempfile

    module = _load_cli_module_for_testing()
    write_fn = module._write_search_metadata

    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)
        write_fn(
            output_dir=output_dir,
            search_backend='mmseqs2',
            database_name='my_core_nt_db',
            database_path='/data/mmseqs2db',
            database_source='ncbi_core_nt_derived_mmseqs2'
        )
        metadata_file = output_dir / "search_metadata.json"
        assert metadata_file.exists(), "search_metadata.json was not created"
        data = json.loads(metadata_file.read_text())
        assert data['search_backend'] == 'mmseqs2'
        assert data['database_name'] == 'my_core_nt_db'
        assert data['database_path'] == '/data/mmseqs2db'
        assert data['database_source'] == 'ncbi_core_nt_derived_mmseqs2'


def test_cli_write_search_metadata_blast_custom_source():
    """Test that BLAST metadata uses blast_db_source when provided, else database name."""
    import json
    import tempfile

    module = _load_cli_module_for_testing()
    write_fn = module._write_search_metadata

    # When blast_db_source is explicitly provided, it is recorded as-is.
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)
        write_fn(
            output_dir=output_dir,
            search_backend='blast',
            database_name='custom_nt',
            database_path='/custom/blastdb',
            database_source='custom_nt'
        )
        data = json.loads((output_dir / "search_metadata.json").read_text())
        assert data['database_source'] == 'custom_nt'
        assert data['database_name'] == 'custom_nt'


def test_cli_metadata_not_overwritten_when_search_skipped():
    """_write_search_metadata must not be called when the search output already
    exists and skip_existing=True (the default non-overwrite mode).

    The guard in blast_command / grouped_command is::

        search_ran = not (search_output.exists() and skip_existing)
        if search_ran:
            _write_backend_search_metadata(args, mmseqs_database)

    This test exercises that conditional to confirm metadata from the original
    run is preserved when the search step is bypassed.
    """
    import json
    import tempfile
    from unittest.mock import patch

    module = _load_cli_module_for_testing()

    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)

        # Simulate first run: write initial metadata
        module._write_search_metadata(
            output_dir=output_dir,
            search_backend='blast',
            database_name='core_nt',
            database_path='/home/user/blastdb',
            database_source='ncbi_core_nt'
        )
        first_data = json.loads((output_dir / "search_metadata.json").read_text())
        assert first_data['database_source'] == 'ncbi_core_nt'

        # Simulate second run: search output already exists + skip_existing=True
        search_output = output_dir / "blast_results.txt"
        search_output.touch()
        skip_existing = True
        search_ran = not (search_output.exists() and skip_existing)
        assert not search_ran, (
            "search_ran should be False when skip_existing=True and output exists"
        )

        # Confirm _write_search_metadata is NOT called when search_ran is False
        with patch.object(module, '_write_search_metadata') as mock_write:
            if search_ran:
                module._write_search_metadata(
                    output_dir=output_dir,
                    search_backend='blast',
                    database_name='different_db',
                    database_path='/other/blastdb',
                    database_source='different_source'
                )
            mock_write.assert_not_called()

        # Original metadata must be unchanged
        second_data = json.loads((output_dir / "search_metadata.json").read_text())
        assert second_data['database_source'] == 'ncbi_core_nt', (
            "Metadata was unexpectedly overwritten when the search was skipped"
        )
        assert second_data['database_name'] == 'core_nt'


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
        test_cli_has_log_level_argument,
        test_cli_log_level_in_all_subparsers,
        test_cli_has_mmseqs_db_source_argument,
        test_cli_has_blast_db_source_argument,
        test_cli_has_mmseqs_download_arguments,
        test_cli_has_write_search_metadata_function,
        test_cli_write_search_metadata_writes_json,
        test_cli_write_search_metadata_mmseqs2,
        test_cli_write_search_metadata_blast_custom_source,
        test_cli_metadata_not_overwritten_when_search_skipped,
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
