all: kalman_walk_paper.pdf ridge_evolution_paper.pdf submission

.PHONY : clean all submission

PDFLATEX=pdflatex -interaction nonstopmode

submission :  kalman_walk_paper_only.pdf cover_letter.pdf review_responses.pdf diff-to-submitted.pdf

kalman_walk_paper_only.pdf : kalman_walk_paper.pdf
	pdfjam --outfile $@ $< 1-25

cover_letter.pdf : kalman_walk_paper.pdf
	pdfjam --outfile $@ $< 26

review_responses.pdf : kalman_walk_paper.pdf
	pdfjam --outfile $@ $< 27-

kalman_walk_paper.pdf : krefs.bib review-responses2.tex

ridge_evolution_paper.pdf : krefs.bib ridge-review-responses.tex

diff-to-submitted.tex : kalman_walk_paper.tex
	latexdiff-git --force -r 9c53627c2f07b94da4f18e8bd6f14dce6e57a47d $<
	mv kalman_walk_paper-diff9c53627c2f07b94da4f18e8bd6f14dce6e57a47d.tex $@

clean : 
	-rm *.aux *.log *.bbl *.blg *.out *.brf

%.html : %.Rmd
	Rscript -e "templater::render_template(\"$<\", output=\"$@\", change.rootdir=TRUE)"

%.pdf : %.tex %.bbl
	while ( $(PDFLATEX) $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.bbl : %.tex krefs.bib
	-$(PDFLATEX) $<
	bibtex $*.aux
