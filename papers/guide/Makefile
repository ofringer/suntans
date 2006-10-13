file=user_guide

all:	pdf

dvi:	
	latex $(file)
	bibtex $(file)
	latex $(file)
	latex $(file)

ps:	dvi
	dvips -Ppdf $(file) -o

pdf:	ps
	ps2pdf $(file).ps $(file).pdf

edit:	all
	xdvi $(file) &
	'emacs' -geometry 120x70 $(file).tex &

clean:
	rm -f *.bbl *.aux *.toc *.log *~ *.aux *.blg *.dvi
