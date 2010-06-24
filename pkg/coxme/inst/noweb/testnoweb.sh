# The Makefile adds a line "Automatically generated from all.nw using noweb"
# That should be the only difference found between the version created
#  using my R functions (via make) and those based on the notangle program.

notangle -Rcoxme all.nw > xxx
make coxme.R
diff xxx coxme.R
rm xxx coxme.R

notangle -RcoxmeMlist all.nw > xxx
make coxmeMlist.R
diff xxx coxmeMlist.R
rm xxx coxmeMlist.R

notangle -Rcoxme.fit all.nw > xxx
make coxme.fit.R
diff xxx coxme.fit.R
rm xxx coxme.fit.R

notangle -Rexpand.nested all.nw > xxx
make expand.nested.R
diff xxx expand.nested.R
rm xxx expand.nested.R

notangle -Rformula all.nw > xxx
make formula.R
diff xxx formula.R
rm xxx formula.R

notangle -RcoxmeFull all.nw > xxx
make coxmeFull.R
diff xxx coxmeFull.R
rm xxx coxmeFull.R

