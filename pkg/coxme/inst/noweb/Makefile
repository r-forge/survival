PARTS = main.Rnw \
	basic.Rnw \
	build.Rnw \
	formula.Rnw \
	varfun.Rnw  varfun2.Rnw \
	fit.Rnw \
	reconcile.Rnw

SFUN =  coxme.R \
	coxme.fit.R\
	expand.nested.R \
	formula.R \
	coxvarFull.R \
	coxvarMlist.R 
#	matrix.reconcile.R

#PARTS = main.Rnw tail

all.pdf: all.tex
	pdflatex all

all.nw: $(PARTS)
	cat $(PARTS) > all.nw
	echo "\\end{document}" >> all.nw

all.tex: all.nw
	noweave -delay all.nw > all.tex

fun:  $(SFUN)

$(SFUN) : all.nw

%.R: 
	echo "# Automatically generated from all.nw using noweb" > $@
	notangle -R$* all.nw >> $@

