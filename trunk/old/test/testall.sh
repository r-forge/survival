#
# Do the full shooting match of tests
#
 cat setup.s > temp
 cat byhand.s >> temp
 cat mktest1.s dotest1.s mktest2.s dotest2.s dotest2b.s >> temp
 cat mkaml.s doaml.s >> temp
 cat stratatest.s >> temp
 cat infcox.s >> temp
 cat mktest1b.s dotest3.s >> temp
 cat doweight1.s doweight2.s >> temp
 cat doweight3.s doweight4.s doweight5.s doweight6.s >> temp
 cat mkturnbull.s doturnbull.s >> temp
 cat read.ovarian.s doovarian.s >> temp
 cat read.jasa.s jasa1.s kalb1.s kalb2.s jasa2.s  >> temp
# cat resid.s >> temp    #takes forever
 cat read.cancer.s docancer.s >> temp
 cat read.bladder.s fit.bladder.s >> temp
 cat doexpect.s  doexpect2.s expect3.s >> temp
 cat testnull.s >> temp
 cat testreg.s dsurvreg.s >> temp
 cat singtest.s >> temp
 cat difftest.s >> temp
 cat coxsurv.s >> temp
 cat pyear1.s  pyear2.s >> temp
 cat read.smoke.s >> temp
 cat survtest.s turntest.s turntest2.s >> temp
 cat survtest2.s >> temp
 cat concordance.s >> temp
 cat aareg1.s aareg2.s aareg3.s >> temp
 cat testci.s gray.s >> temp
 echo 'q()' >> temp
 rm -rf .Data
 Splus CHAPTER
 Splus BATCH temp testall.out &
