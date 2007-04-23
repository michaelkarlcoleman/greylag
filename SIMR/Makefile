# Makefile for SIMR-local examples

.PHONY: install
.DELETE_ON_ERROR:


DEST = /n/site/inst/Linux-i686/bioinfo/greylag/

install::
	[ -d $(DEST) ] || install -d $(DEST)
	install -p greylag-process greylag-sge-slice $(DEST)/
