work-rbras2018: work-rbras2018.Rnw rbras2018.bib codes/*

	Rscript -e 'knitr::knit("work-rbras2018.Rnw")'
	pdflatex work-rbras2018.tex
	-bibtex work-rbras2018.aux
	pdflatex work-rbras2018.tex
	pdflatex work-rbras2018.tex

simple:
	Rscript -e 'knitr::knit("work-rbras2018.Rnw")'
	pdflatex work-rbras2018.tex

clean:
	rm -f *.aux *.bbl *.blg *.brf *.idx *.ilg *.ind *.lof *.log \
	.*lol *.lot *.out *.toc *.synctex.gz
