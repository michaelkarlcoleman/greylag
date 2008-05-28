# Makefile for greylag doc


.PHONY: all html pdf install publish clean
.DELETE_ON_ERROR:

DEST = /usr/local/lib/greylag/doc/

DOCS = index.txt release-notes.txt install.txt user-guide.txt \
       theory-of-operation.txt

IMAGES = favicon.png
HTML = $(addsuffix .html,$(basename $(DOCS)))
PDF = $(addsuffix .pdf,$(basename $(DOCS)))
ALL = $(DOCS) $(IMAGES) $(HTML) $(PDF)

RST2HTML = rst2html --exit-status=warning --generator --date --source-link \
	--initial-header-level=2 --cloak-email-addresses

RST2LATEX = rst2latex --exit-status=warning --documentoptions=letterpaper

all :: $(ALL)

html :: $(HTML)
pdf :: $(PDF)


%.html : %.txt
	$(RST2HTML) $^ $@.tmp
	sed -e 's|</title>|</title><link rel="shortcut icon" href="favicon.png" />|' $@.tmp > $@

%.pdf : %.tex
	pdflatex $^

%.tex : %.txt
	$(RST2LATEX) $^ $@


install :: $(ALL)
	[ -d $(DEST) ] || install -d $(DEST)
	for f in $(ALL); do \
	    install -p --mode=444 $$f $(DEST); \
	done



# push docs to sourceforge web site
publish :: $(ALL)
	rsync -v -t $(ALL) shell.sf.net:greylag/htdocs/

clean ::
	-rm -f $(HTML) *.tmp $(PDF) *.tex *.log *.dvi *.out *.aux *~

