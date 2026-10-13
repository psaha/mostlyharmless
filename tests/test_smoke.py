"""Stage 0 smoke tests: the harness works and the environment is adequate.

These are not certificates.  They exist so that a green run means something
before any real mathematics is in place.
"""

import importlib
import sys

import pytest

# Minimum (major, minor) for each dependency, matching pyproject.toml.
MINIMUM = {
    "numpy": (2, 2),
    "scipy": (1, 15),
    "sympy": (1, 13),
}

SUBPACKAGES = [
    "mostlyharmless",
    "mostlyharmless.visylm",
    "mostlyharmless.symbolic",
]


def _major_minor(version):
    """Parse the first two components of a version string, ignoring suffixes."""
    parts = []
    for chunk in version.split(".")[:2]:
        digits = "".join(ch for ch in chunk if ch.isdigit())
        parts.append(int(digits or 0))
    return tuple(parts)


def test_interpreter_is_at_least_the_target_version():
    assert sys.version_info[:2] >= (3, 10), sys.version


def test_pytest_is_recent_enough_for_the_pythonpath_option():
    # `pythonpath` in [tool.pytest.ini_options] needs pytest >= 7.  It is what
    # lets pytest run against the src layout without installing anything.
    assert _major_minor(pytest.__version__)[0] >= 7, pytest.__version__


@pytest.mark.parametrize("package", SUBPACKAGES)
def test_package_imports(package):
    assert importlib.import_module(package) is not None


@pytest.mark.parametrize("module", sorted(MINIMUM))
def test_scientific_stack_meets_minimum_version(module):
    mod = importlib.import_module(module)
    have = _major_minor(mod.__version__)
    assert have >= MINIMUM[module], f"{module} {mod.__version__} < {MINIMUM[module]}"
