PAPERFILENAME:=Manuscript
SIMATERIAL:=SupportingInformation
BIBFILE:=references.bib
PDF:=$(PAPERFILENAME).pdf
SIPDF:=$(SIMATERIAL).pdf
IMAGES:=figures/fig1.pdf

all: $(PDF) $(SIPDF)

clean:
	latexmk -c $(PDF)
	latexmk -c $(SIPDF)
	rm -f $(PDF) $(SIPDF) $(PAPERFILENAME).tex $(SIMATERIAL).tex

%.pdf: %.Rmd
	R --no-save --quiet -e "rmarkdown::render('$<')"

$(SIMATERIAL).tex: $(SIMATERIAL).pdf

$(PAPERFILENAME).pdf: $(PAPERFILENAME).Rmd $(BIBFILE) $(IMAGES) $(SIMATERIAL).tex $(SIPDF) preamble.tex

$(SIMATERIAL).pdf: $(SIMATERIAL).Rmd $(BIBFILE) sipreamble.tex

