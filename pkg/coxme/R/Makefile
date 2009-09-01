SFUN =  coxme.R \
        coxme.fit.R\
        expand.nested.R \
        formula.R \
        coxvarFull.R \
        coxvarMlist.R


fun:  $(SFUN)

$(SFUN) : ../inst/noweb/all.nw

%.R: 
	echo "# Automatically generated from all.nw using noweb" > $@
	notangle -R$* ../inst/noweb/all.nw >> $@

clean:
	rm $(SFUN)
