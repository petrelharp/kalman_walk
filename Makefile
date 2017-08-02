.PHONY : clean

PDFLATEX=pdflatex -interaction nonstopmode

clean : 
	-rm *.aux *.log *.bbl *.blg *.out *.brf

%.html : %.Rmd
	Rscript -e "templater::render_template(\"$<\", output=\"$@\", change.rootdir=TRUE)"

%.pdf : %.tex %.bbl
	while ( $(PDFLATEX) $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.bbl : %.tex
	-$(PDFLATEX) $<
	bibtex $*.aux

handnotes.pdf : krefs.bib
