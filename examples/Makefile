.PHONY : htmls

# don't delete intermediate files
.SECONDARY :

# MATHJAX = /usr/share/javascript/mathjax/MathJax.js
MATHJAX = https://cdn.mathjax.org/mathjax/latest/MathJax.js
PANDOC_OPTS =  --to html --from markdown-implicit_figures --self-contained --standalone --section-divs --template /usr/local/lib/R/site-library/rmarkdown/rmd/h/default.html --variable 'theme:bootstrap' --include-in-header resources/header-scripts.html --mathjax --variable 'mathjax-url:$(MATHJAX)?config=TeX-AMS-MML_HTMLorMML' --no-highlight --variable highlightjs=/usr/local/lib/R/site-library/rmarkdown/rmd/h/highlight 

mds : $(patsubst %.Rmd,%.md,$(wildcard *.Rmd))
htmls : $(patsubst %.md,%.html,$(wildcard *.md) $(mds))

%.html : %.Rmd
	cd $$(dirname $<); Rscript -e 'templater::render_template(basename("$<"),output=basename("$@"),clean=FALSE)'

%.md : %.Rmd
	# cd $$(dirname $<); Rscript -e 'knitr::knit(basename("$<"),output=basename("$@"))'
	cd $$(dirname $<); Rscript -e 'templater::render_template(basename("$<"),output=basename("$@"),clean=FALSE)'

%.html : %.md
	pandoc $< $(PANDOC_OPTS) --output $@

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@
