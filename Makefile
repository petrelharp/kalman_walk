.PHONY : clean

clean : 
	-rm *.{aux,log,bbl,blg,out}

%.html : %.Rmd
	Rscript -e "templater::render_template(\"$<\", output=\"$@\", change.rootdir=TRUE)"
