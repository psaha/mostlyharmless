# mostlyharmless -- build targets.
#
#   make check   run the certificates (pytest)
#   make docs    build the manual
#   make all     same as docs; deliberately NOT strict, see decision D6
#
# PYTHON defaults to the Thonny interpreter that Stage 0 was developed against.
# Override on the command line if you are working elsewhere:
#
#   make check PYTHON=/path/to/python3

PYTHON   ?= /Users/psaha/Applications/Thonny.app/Contents/Frameworks/Python.framework/Versions/3.10/bin/python3.10
PDFLATEX ?= pdflatex

DOC := docs/mostlyharmless.tex
PDF := docs/mostlyharmless.pdf
TEX := $(wildcard docs/sections/*.tex)

.PHONY: all check docs data clean distclean

all: docs

check:
	$(PYTHON) -m pytest

docs: $(PDF)

$(PDF): $(DOC) $(TEX)
	cd docs && $(PDFLATEX) -interaction=nonstopmode -halt-on-error mostlyharmless.tex
	cd docs && $(PDFLATEX) -interaction=nonstopmode -halt-on-error mostlyharmless.tex

# The coefficient table and the emitted LaTeX fragments arrive in Stage 4.
data:
	@echo "nothing to emit yet -- the coefficient table arrives in Stage 4"

clean:
	rm -f docs/*.aux docs/*.log docs/*.out docs/*.toc docs/*.bbl docs/*.blg
	rm -rf .pytest_cache
	find src tests -name __pycache__ -type d -prune -exec rm -rf {} +

distclean: clean
	rm -f $(PDF)
