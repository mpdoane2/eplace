Contributing to ePLACE
=======================

Thank you for your interest in contributing to ePLACE! This document provides guidelines for contributing to the project.

Getting Started
---------------

1. Fork the repository on GitHub
2. Clone your fork locally:

   .. code-block:: bash

      git clone https://github.com/YOUR_USERNAME/eplace.git
      cd eplace

3. Create a branch for your changes:

   .. code-block:: bash

      git checkout -b feature/your-feature-name

Development Setup
-----------------

Install development dependencies:

.. code-block:: bash

   # Install in development mode with dev dependencies
   pip install -e ".[dev]"

   # Install external dependencies
   mamba install -y bioconda::blast bioconda::pytaxonkit bioconda::iqtree bioconda::mafft

Running Tests
-------------

Run the test suite to ensure your changes don't break existing functionality:

.. code-block:: bash

   # Run all tests
   pytest tests/ -v

   # Run specific test modules
   pytest tests/test_blast_analysis.py -v
   pytest tests/test_taxonomy.py -v
   pytest tests/test_workflow.py -v

   # Run with coverage
   pytest tests/ --cov=eplace_lib --cov-report=html

Code Style
----------

* Follow PEP 8 style guidelines
* Use meaningful variable and function names
* Add type hints to function signatures
* Write docstrings for all public functions and classes

Docstring Format
~~~~~~~~~~~~~~~~

Use Google-style docstrings:

.. code-block:: python

   def function_name(param1: str, param2: int) -> bool:
       """Brief description of function.

       More detailed description if needed.

       Args:
           param1: Description of param1
           param2: Description of param2

       Returns:
           Description of return value

       Raises:
           ValueError: When invalid input is provided
       """
       pass

Writing Tests
-------------

* Add tests for new features
* Ensure tests are independent and can run in any order
* Use descriptive test names
* Include both positive and negative test cases

Example test structure:

.. code-block:: python

   def test_feature_name():
       """Test that feature works correctly."""
       # Setup
       input_data = create_test_data()
       
       # Execute
       result = function_under_test(input_data)
       
       # Assert
       assert result == expected_output

Documentation
-------------

* Update documentation for any new features
* Add examples to relevant documentation pages
* Update API documentation if adding new modules or functions
* Keep README.md up to date

Building Documentation
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   cd docs
   make html

   # View documentation
   # Open docs/build/html/index.html in your browser

Pull Request Process
--------------------

1. Ensure all tests pass
2. Update documentation as needed
3. Add a clear description of your changes
4. Reference any related issues
5. Request review from maintainers

Pull Request Guidelines
~~~~~~~~~~~~~~~~~~~~~~~

* Keep changes focused and atomic
* Write clear commit messages
* Update CHANGELOG if applicable
* Ensure CI/CD checks pass

Commit Messages
---------------

Write clear, descriptive commit messages:

.. code-block:: text

   Short summary (50 chars or less)

   More detailed explanation if necessary. Wrap at 72 characters.
   Explain the problem that this commit is solving, and why this
   approach was chosen.

   - Bullet points are okay
   - Use present tense ("Add feature" not "Added feature")
   - Reference issues and pull requests

Reporting Bugs
--------------

When reporting bugs, please include:

* Operating system and version
* Python version
* ePLACE version
* Steps to reproduce the issue
* Expected behavior
* Actual behavior
* Error messages or logs
* Minimal reproducible example if possible

Use the GitHub issue tracker: https://github.com/linsalrob/eplace/issues

Feature Requests
----------------

We welcome feature requests! Please:

* Check if the feature has already been requested
* Clearly describe the feature and its use case
* Explain why this feature would be useful
* Provide examples of how it would be used

Code Review Process
-------------------

All submissions require review. We use GitHub pull requests for this purpose.

* Maintainers will review your code
* Address any feedback or requested changes
* Once approved, your changes will be merged

Areas for Contribution
----------------------

Some areas where contributions are especially welcome:

* Bug fixes
* Documentation improvements
* Additional test coverage
* Performance improvements
* Support for additional BLAST databases
* New analysis features
* User interface improvements

Code of Conduct
---------------

* Be respectful and inclusive
* Welcome newcomers
* Focus on constructive feedback
* Assume good intentions

Contact
-------

* GitHub Issues: https://github.com/linsalrob/eplace/issues
* Email: raedwards@gmail.com

License
-------

By contributing to ePLACE, you agree that your contributions will be licensed under the MIT License.

Thank You!
----------

Your contributions are greatly appreciated and help make ePLACE better for everyone!
