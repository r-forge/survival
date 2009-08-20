SFUN =  coxme.R \
        coxme.fit.R\
        expand.nested.R \
        formula.R \
        coxvarFull.R \
        coxvarMlist.R

NWsrc = ../inst/noweb   

fun:  $(SFUN)

$(SFUN) : $(NWsrc)/all.nw

$(NWsrc)/all.nw:
	$(MAKE) -C $(NWsrc) all.nw

%.R: 
	echo "# Automatically generated from all.nw using noweb" > $@
	notangle -R$* ../inst/noweb/all.nw >> $@

clean:
	rm $(SFUN)
	rm ../inst/noweb/all.nw
