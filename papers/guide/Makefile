file=user_guide

all:	
	make dvi
	make bib
	make dvi
	make dvi

dvi:	$(file).tex
	latex $(file)

bib:	bibliography.bib
	bibtex $(file)

ps:	dvi
	dvips $(file) -o

pdf:	dvi bib
	pdflatex $(file)

edit:	all
	xdvi $(file) &
	'emacs' -geometry 120x70 $(file).tex &

clean:
	rm -f *.bbl *.aux *.toc *.log *~ *.aux *.blg *.dvi
