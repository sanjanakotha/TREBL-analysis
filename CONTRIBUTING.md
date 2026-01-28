# Contributing to TREBL Analysis

Thank you for your interest in contributing to TREBL Analysis! This document provides guidelines and instructions for contributing to the project.

## Getting Started

### Prerequisites

- Python 3.10 or higher
- Git
- Basic understanding of Python and bioinformatics workflows

### Setting Up Your Development Environment

1. **Fork the repository** on GitHub

2. **Clone your fork locally**:
   ```bash
   git clone https://github.com/YOUR_USERNAME/TREBL-analysis.git
   cd TREBL-analysis
   ```

3. **Add the upstream repository**:
   ```bash
   git remote add upstream https://github.com/sanjanakotha/TREBL-analysis.git
   ```

4. **Install dependencies**:
   ```bash
   pip install -r docs/requirements.txt
   ```

5. **Create a branch for your changes**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## How to Contribute

### Reporting Bugs

If you find a bug, please create an issue on GitHub with:

- A clear, descriptive title
- Steps to reproduce the bug
- Expected behavior
- Actual behavior
- Your environment (OS, Python version, package versions)
- Any relevant error messages or logs

### Suggesting Enhancements

Feature requests are welcome! Please create an issue with:

- A clear description of the enhancement
- Use cases and benefits
- Any implementation ideas (optional)

### Pull Requests

1. **Keep changes focused**: Each PR should address a single concern
2. **Follow the existing code style**: See Code Style Guidelines below
3. **Update documentation**: Add or update docstrings and documentation files
4. **Test your changes**: Ensure your code works as expected
5. **Write clear commit messages**: Use descriptive messages explaining what and why

#### Pull Request Process

1. Update your fork with the latest upstream changes:
   ```bash
   git fetch upstream
   git rebase upstream/main
   ```

2. Make your changes in your feature branch

3. Commit your changes:
   ```bash
   git add .
   git commit -m "Brief description of changes"
   ```

4. Push to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

5. Open a Pull Request on GitHub with:
   - A descriptive title
   - A detailed description of changes
   - References to related issues (e.g., "Fixes #123")

## Code Style Guidelines

### Python Code

- **Follow PEP 8**: Use Python's style guide for code formatting
- **Line length**: Maximum 100 characters per line
- **Imports**: Group imports in the following order:
  1. Standard library imports
  2. Third-party imports
  3. Local application imports
- **Naming conventions**:
  - Functions and variables: `snake_case`
  - Classes: `PascalCase`
  - Constants: `UPPER_CASE`

### Documentation

#### Docstrings

Use Google-style docstrings for all functions and classes:

```python
def example_function(param1, param2):
    """
    Brief description of what the function does.
    
    More detailed explanation if needed.
    
    Args:
        param1 (type): Description of param1.
        param2 (type): Description of param2.
        
    Returns:
        type: Description of return value.
        
    Raises:
        ExceptionType: When and why this exception is raised.
        
    Example:
        >>> result = example_function("hello", 42)
        >>> print(result)
        Expected output
    """
    # Function implementation
    pass
```

#### Code Comments

- Use comments to explain **why**, not **what**
- Keep comments up-to-date with code changes
- Avoid obvious comments

### reStructuredText Documentation

For Sphinx documentation files (`.rst`):

- Use clear section headers
- Include code examples where appropriate
- Cross-reference related documentation with `:doc:` and `:class:` directives

## Documentation

### Building Documentation Locally

To preview documentation changes:

```bash
cd docs
make html
```

View the generated documentation at `docs/_build/html/index.html`.

### Adding New Documentation

1. Create new `.rst` files in the `docs/` directory
2. Add references to new files in `docs/index.rst` or other relevant files
3. Build and preview locally to ensure proper rendering

## Testing

While this project currently has limited automated testing, you should:

1. **Test your changes manually**: Run scripts with test data
2. **Check for errors**: Ensure no Python exceptions occur
3. **Verify output**: Confirm results match expectations
4. **Test edge cases**: Try unusual inputs or boundary conditions

### Manual Testing Checklist

- [ ] Code runs without errors
- [ ] Results match expected output
- [ ] Documentation builds without warnings
- [ ] Examples in documentation work correctly
- [ ] Changes don't break existing functionality

## Project Structure

Understanding the project structure will help you contribute effectively:

```
TREBL-analysis/
├── scripts/              # Core analysis modules
│   ├── __init__.py      # Package initialization
│   ├── finder.py        # Barcode definitions
│   ├── initial_map.py   # Sequence mapping
│   ├── map_refiner.py   # Data refinement
│   ├── pipelines.py     # Complete workflows
│   ├── error_correct.py # Error correction
│   ├── umi_deduplicate.py # UMI handling
│   ├── plotting.py      # Visualization
│   ├── preprocess.py    # Data preprocessing
│   └── complexity.py    # Complexity analysis
├── notebooks/           # Jupyter notebooks
├── docs/               # Documentation source
│   ├── conf.py         # Sphinx configuration
│   ├── index.rst       # Documentation index
│   ├── step1_protocol.rst
│   └── scripts.rst     # API documentation
├── savio_jobs/         # HPC job scripts
└── README.md
```

## Module Guidelines

### Adding New Modules

When adding a new module to `scripts/`:

1. Create the module file: `scripts/your_module.py`
2. Add module docstring at the top
3. Import in `scripts/__init__.py` if needed
4. Add API documentation to `docs/scripts.rst`
5. Update `docs/index.rst` if adding new major features

### Modifying Existing Modules

- Maintain backward compatibility when possible
- Document breaking changes clearly
- Update related documentation
- Consider adding deprecation warnings before removing features

## Commit Message Guidelines

Write clear, concise commit messages:

- **Format**: `<type>: <subject>`
- **Types**:
  - `feat:` New feature
  - `fix:` Bug fix
  - `docs:` Documentation changes
  - `style:` Code style changes (formatting, etc.)
  - `refactor:` Code refactoring
  - `test:` Adding or updating tests
  - `chore:` Maintenance tasks

**Examples**:
```
feat: Add error correction option to pipelines
fix: Handle empty barcode sequences correctly
docs: Update installation instructions
refactor: Simplify MapRefiner initialization
```

## Questions?

If you have questions about contributing:

1. Check existing documentation
2. Search through existing issues
3. Open a new issue with the "question" label
4. Contact the maintainers

## Code of Conduct

### Our Pledge

We are committed to providing a welcoming and inclusive environment for all contributors.

### Expected Behavior

- Be respectful and considerate
- Accept constructive criticism gracefully
- Focus on what's best for the project
- Show empathy towards other contributors

### Unacceptable Behavior

- Harassment or discriminatory language
- Personal attacks
- Publishing others' private information
- Other unprofessional conduct

### Enforcement

Violations may result in temporary or permanent removal from the project. Report issues to the project maintainers.

## Recognition

Contributors will be acknowledged in:
- The project README
- Release notes
- Documentation

Thank you for contributing to TREBL Analysis!
