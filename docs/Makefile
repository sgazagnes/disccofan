# Makefile for Sphinx documentation
# See sphinx-quickstart for more information

# The default theme is "default". You can change it if you like.
SPHINXOPTS    = 
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = $(shell basename $(PWD))
SOURCEDIR     = source
BUILDDIR      = build

.PHONY: help Makefile

help:
	@echo "Please use \`make <target>\` where <target> is one of"
	@echo "  html    to make HTML files (default)"
	@echo "  help    to show this help message"

html:
	@$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR) $(SPHINXOPTS) $(O)
