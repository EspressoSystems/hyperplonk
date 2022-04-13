PROJECT_NAME = main
IMG_OUTPUT=images/output

.PHONY: all
all:
	@echo "\n\n Compiling document \n\n"
	pdflatex $(PROJECT_NAME).tex
	bibtex $(PROJECT_NAME)
	pdflatex $(PROJECT_NAME).tex
	pdflatex -halt-on-error $(PROJECT_NAME).tex

.PHONY: simple
simple:
	pdflatex -halt-on-error $(PROJECT_NAME).tex

.PHONY: images
images: img/*.dot
	dot img/algorithms.dot -Tpng > $(IMG_OUTPUT)/algorithms.png

.PHONY: watch
watch:
	latexmk -quiet -bibtex -pvc -f -pdf -pdflatex="pdflatex -synctex=1 -interaction=nonstopmode" ${PROJECT_NAME}.tex

.PHONY: clean
clean:
	@echo "\n\n\n*****Cleaning... *******"
	rm -f *.aux* *.log *.ps *.toc *.lof *.lot *~ *.dvi *.bbl *.blg *.out *.glo *.ist *.acn *.acr *.gls *.alg *.glg *.bak *bcf *xml *.brf *.fls *.fdb_latexmk
	# rm -f $(PROJECT_NAME).pdf
	rm -f $(IMG_OUTPUT)/*.png
