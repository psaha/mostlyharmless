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
FIGURES := docs/generated/map.pdf docs/generated/uv.pdf docs/generated/convergence.pdf

.PHONY: all check docs data gui clean distclean

all: check docs

check:
	$(PYTHON) -m pytest -rP

docs: $(PDF)

$(PDF): $(DOC) $(TEX) $(GENERATED) $(FIGURES) docs/refs.bib
	cd docs && $(PDFLATEX) -interaction=nonstopmode -halt-on-error mostlyharmless.tex
	cd docs && $(BIBTEX) mostlyharmless
	cd docs && $(PDFLATEX) -interaction=nonstopmode -halt-on-error mostlyharmless.tex
	cd docs && $(PDFLATEX) -interaction=nonstopmode -halt-on-error mostlyharmless.tex

# The emitted LaTeX fragments and figures.  The output is committed, so a fresh
# checkout builds without this; the rules fire when a file is missing.
data:
	PYTHONPATH=src $(PYTHON) -m mostlyharmless.symbolic.emit
	PYTHONPATH=src $(PYTHON) -m mostlyharmless.symbolic.figures

# The interactive front end (Stage 8, D25).  Needs a display; not part of `all`.
gui:
	PYTHONPATH=src $(PYTHON) -m mostlyharmless.symbolic.gui

$(GENERATED) $(FIGURES):
	$(MAKE) data

clean:
	rm -f docs/*.aux docs/*.log docs/*.out docs/*.toc docs/*.bbl docs/*.blg
	rm -rf .pytest_cache
	find src tests -name __pycache__ -type d -prune -exec rm -rf {} +

distclean: clean
	rm -f $(PDF)
