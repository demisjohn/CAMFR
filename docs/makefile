all: pdf html

dvi: camfr.texi
	texi2dvi camfr.texi

pdf: dvi
	dvipdf -sPAPERSIZE=a4 camfr

html:
	texi2html -init_file texi2html.init -split section camfr.texi

clean:
	rm -f *~ *.aux *.cp *.cps *.dvi *.fn *.ky *.log *.pdf *.pg *.toc \
	*.tp *.vr *.html
