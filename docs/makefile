all: pdf html

dvi: camfr.texi
	texi2dvi camfr.texi

ps: dvi
	dvips -Ppdf -ta4 camfr

pdf: ps
	ps2pdf camfr.ps

html:
	texi2html -init_file texi2html.init -split section camfr.texi

clean:
	rm -f *~ *.aux *.cp *.cps *.dvi *.fn *.ky *.log *.pdf *.pg *.toc \
	*.tp *.vr *.html *.ps

