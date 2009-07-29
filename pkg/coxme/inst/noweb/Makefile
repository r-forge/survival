PARTS = main.Rnw \
	basic.Rnw \
	build.Rnw \
	formula.Rnw \
	varfun.Rnw \
	fit.Rnw 

SFUN =  coxme.R \
	coxme.fit.R\
	expand.nested.R \
	formula.R \
	coxvarFull.R \
	coxvarMlist.R

#PARTS = main.Rnw tail

all.pdf: all.tex
	pdflatex all

all.nw: $(PARTS)
	cat $(PARTS) > all.nw
	echo "\\end{document}" >> all.nw

all.tex: all.nw
	noweave -delay all.nw > all.tex

fun:  $(SFUN)
	touch fun

$(SFUN) : all.nw

%.R: 
	echo "# Automatically generated from all.nw using noweb\n" > $@
	notangle -R$* all.nw >> $@

