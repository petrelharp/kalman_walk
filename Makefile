all: kalman_walk_paper.pdf ridge_evolution_paper.pdf 

.PHONY : clean all

PDFLATEX=pdflatex -interaction nonstopmode

kalman_walk_paper.pdf : krefs.bib review-responses.tex

ridge_evolution_paper.pdf : krefs.bib ridge-review-responses.tex

clean : 
	-rm *.aux *.log *.bbl *.blg *.out *.brf

%.html : %.Rmd
	Rscript -e "templater::render_template(\"$<\", output=\"$@\", change.rootdir=TRUE)"

%.pdf : %.tex %.bbl
	while ( $(PDFLATEX) $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.bbl : %.tex krefs.bib
	-$(PDFLATEX) $<
	bibtex $*.aux
