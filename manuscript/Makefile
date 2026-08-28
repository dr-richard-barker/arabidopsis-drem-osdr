# Build the manuscript as PDF and Word from ONE LaTeX source.
#
#   make            manuscript.pdf
#   make docx       manuscript.docx (pandoc, same source, same bibliography)
#   make numbers    regenerate generated_numbers.tex from results/
#   make check      fail if a citation is undefined or a DOI is wrong
#   make clean      remove build artefacts
#
# Both targets read latex/main.tex, so the PDF and the Word file cannot drift apart.

TEX  := main
DIR  := latex
BIB  := $(DIR)/references.bib
SRC  := $(DIR)/$(TEX).tex $(wildcard $(DIR)/sections/*.tex) $(DIR)/generated_numbers.tex
PY   := python3

.PHONY: all pdf docx numbers references check clean

all: pdf

pdf: manuscript.pdf
docx: manuscript.docx

# generated_numbers.tex is a build product: every number in the prose comes from
# results/, never from the author's memory.
numbers:
	cd .. && $(PY) scripts/13_manuscript_numbers.py

references:
	cd .. && $(PY) scripts/check_references.py

manuscript.pdf: $(SRC) $(BIB)
	cd $(DIR) && latexmk -pdf -interaction=nonstopmode -halt-on-error $(TEX).tex
	cp $(DIR)/$(TEX).pdf manuscript.pdf

# --citeproc renders the bibliography into the .docx, so the Word file is
# self-contained for a co-author who has no LaTeX.
manuscript.docx: $(SRC) $(BIB)
	cd $(DIR) && pandoc $(TEX).tex \
	  --from=latex --citeproc --bibliography=references.bib \
	  --resource-path=.:sections --toc \
	  --default-image-extension=png \
	  -o ../manuscript.docx

check:
	cd .. && $(PY) scripts/check_references.py --check
	cd .. && $(PY) scripts/check_citations.py

clean:
	cd $(DIR) && latexmk -C $(TEX).tex 2>/dev/null || true
	rm -f $(DIR)/*.aux $(DIR)/*.bbl $(DIR)/*.blg $(DIR)/*.log $(DIR)/*.out \
	      $(DIR)/*.fdb_latexmk $(DIR)/*.fls
