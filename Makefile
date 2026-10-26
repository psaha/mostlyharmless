# mostlyharmless -- build targets.
#
#   make check   run the certificates (pytest)
#   make docs    build the manual; plain LaTeX, no Python unless a generated file is missing
#   make all     the default: check first, then build, stopping if the certificates fail (D6)
#
# PYTHON defaults to the Thonny interpreter that Stage 0 was developed against.
# Override on the command line if you are working elsewhere:
#
#   make check PYTHON=/path/to/python3

PYTHON   ?= /Users/psaha/Applications/Thonny.app/Contents/Frameworks/Python.framework/Versions/3.10/bin/python3.10
PDFLATEX ?= pdflatex
BIBTEX   ?= bibtex

DOC := docs/mostlyharmless.tex
PDF := docs/mostlyharmless.pdf
TEX := $(wildcard docs/sections/*.tex)
GENERATED := docs/generated/coefficients.tex

.PHONY: all check check-fast retest docs data earth gui mapgui clean distclean

all: check docs

# The full suite, logged so that a run can be inspected afterwards without re-running it.
check:
	$(PYTHON) -m pytest -rP -n auto --durations=10 2>&1 | tee tests/check.log

# Everything except the three slow modules: the affine projection, the fitted row's
# least-squares fits and the deformation GUI's sampled renders.  The map GUI and the Earth
# stay in --- its renders are 61-pixel panels and its coefficients are committed --- so
# this is the loop to use while editing.  --ignore takes one path per flag, hence the
# addprefix rather than a comma-separated list.
SLOW := tests/test_affine.py tests/test_fitted.py tests/test_gui.py
check-fast:
	$(PYTHON) -m pytest -rP -n auto $(addprefix --ignore=,$(SLOW)) 2>&1 | tee tests/check-fast.log

# Only the tests that failed last time.  pytest's own cache, so `make check` has to have
# run at least once; it is the quickest way to see whether a fix landed.
retest:
	$(PYTHON) -m pytest -rP -n auto --lf 2>&1 | tee tests/check-retest.log

docs: $(PDF)

$(PDF): $(DOC) $(TEX) $(GENERATED) docs/refs.bib
	cd docs && $(PDFLATEX) -interaction=nonstopmode -halt-on-error mostlyharmless.tex
	cd docs && $(BIBTEX) mostlyharmless
	cd docs && $(PDFLATEX) -interaction=nonstopmode -halt-on-error mostlyharmless.tex
	cd docs && $(PDFLATEX) -interaction=nonstopmode -halt-on-error mostlyharmless.tex

# The emitted LaTeX fragments and figures.  The output is committed, so a fresh
# checkout builds without this; the rules fire when a file is missing.
data:
	PYTHONPATH=src $(PYTHON) -m mostlyharmless.symbolic.emit

# The interactive front end (Stage 8, D25).  Needs a display; not part of `all`.
gui:
	PYTHONPATH=src $(PYTHON) -m mostlyharmless.symbolic.gui

mapgui:
	PYTHONPATH=src $(PYTHON) -m mostlyharmless.symbolic.map_gui

# Re-extract the Earth's colour coefficients from the equirectangular image in misc/.
# Run by hand: it needs matplotlib only, writes a committed module, and is the one
# step that reads a picture.
earth:
	PYTHONPATH=src $(PYTHON) -m mostlyharmless.symbolic.earth misc/Blue_Marble_2002.png --lmax 16

$(GENERATED):
	$(MAKE) data

clean:
	rm -f docs/*.aux docs/*.log docs/*.out docs/*.toc docs/*.bbl docs/*.blg
	rm -rf .pytest_cache tests/check*.log
	find src tests -name __pycache__ -type d -prune -exec rm -rf {} +

distclean: clean
	rm -f $(PDF)
