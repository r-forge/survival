# To add a new set of functions --
#   Add any ".o" files that it uses to the OBJ list
#   Add a line for each function to the FUN list
#   Add a line for each documented function (should be all of them!) to DOC
# If the group is large, you may want to separate it out just to help keep
#   this file more readable, as is done with "SURV" below.
#
#   Use  "Splus make" or "Splus make all" to make a local copy of S.  You may want to
#     run S in this directory just to make sure all of your changes work.
#   DO NOT run "Splus CHAPTER" as this will overwrite this Makefile with one
#     auto-generated by Splus.  This Makefile uses "Splus CHAPTER -m" which does
#     not overwrite the Makefile.
#   make fun - functions only
#   make doc - documentation only
#

# Additional compilation flags
#  The ones below are perhaps temporary, they avoid certain symbol conflicts
#    between my internal routine names and Splus.  With S5's dyn.open, I
#    can't override routines that are called by other C routines.  Any
#    that get changed go here

OBJ =   coxfit2.o \
	agexact.o   \
	agfit3.o   \
	agfit5.o   \
	agmart.o \
	agmart2.o \
	agscore.o \
	agsurv1.o  \
	agsurv2.o  \
	agsurv3.o   \
	cholesky2.o \
	cholesky3.o \
	chsolve2.o  \
	chsolve3.o  \
	chinv2.o   \
	chinv3.o   \
	cox_callback.o \
	coxdetail.o \
	coxfit5.o \
	coxmart.o \
	coxph_wtest.o \
	coxscore.o \
	coxscho.o  \
	dmatrix.o  \
	doloop.o  \
	pyears1.o \
	pyears2.o \
	pyears3.o \
	pystep.o  \
	surv_callback.o \
	survConcordance.o \
	survdiff2.o \
	survfit4.o   \
	survreg2.o  \
	survreg3.o  \
	survreg4.o  \
	survreg5.o

FUN =	aareg                  \
	aareg.taper		\
	agexact.fit            \
	agreg.fit              \
	cluster                 \
	coxph                   \
	coxph.fit              \
	coxph.control		\
	coxph.detail           \
	coxph.getdata          \
	cox.zph                \
	coxpenal.fit	\
	coxpenal.df	\
	coxph.wtest	\
	dsurvreg	\
	format.Surv            \
	frailty         \
	frailty.brent	\
	frailty.controlaic	\
	frailty.controldf	\
	frailty.controlgam	\
	frailty.controlgauss	\
	frailty.gamma		\
	frailty.gammacon	\
	frailty.gaussian	\
	frailty.t		\
	is.na.coxph.penalty	\
	is.ratetable           \
	plot.aareg	\
	plot.cox.zph            \
	print.aareg            \
	print.coxph.penal	\
	print.coxph            \
	print.ratetable         \
	print.survdiff         \
	pspline			\
	pyears			\
	survexp		\
	survdiff             \
	survdiff.fit         \
	survobrien  \
	residuals.coxph      \
	residuals.coxph.penal \
	summary.coxph        \
	summary.coxph.penal   \
	Surv                   \
	survfit              \
	survfit.ci	\
	labels.survreg		\
	lines.aareg		\
	lines.survfit          \
	match.ratetable         \
	model.frame.survreg    \
	model.frame.coxph      \
	model.newframe         \
	plot.survfit           \
	points.survfit         \
	predict.coxph          \
	predict.coxph.penal	\
	predict.survreg        \
	predict.survreg.penal  \
	print.coxph.null       \
	print.cox.zph           \
	print.survexp          \
	print.survfit          \
	print.survreg          \
	print.survreg.penal    \
	print.summary.coxph	\
	print.summary.survreg  \
	print.summary.survfit  \
	pyears                  \
	ratetable               \
	residuals.coxph.null \
	residuals.survreg    \
	residuals.survreg.penal \
	strata               \
	summary.aareg	     \
	summary.pyears	     \
	summary.ratetable    \
	summary.survfit      \
	summary.survreg      \
	survexp.fit          \
	survexp.cfit         \
	survConcordance      \
	survfit.km           \
	survfit.ci	    \
	survfit.coxph.null   \
	survfit.coxph        \
	survfit.turnbull     \
	survpenal.fit	     \
	survreg              \
	survreg.control      \
	survreg.distributions\
	survreg.fit          \
	tcut                 \
	untangle.specials    \
	$(FIX)

RATES=	survexp.us           \
	survexp.usr          \
	survexp.mn	\
	survexp.mnwhite		\
	survexp.fl	   \
	survexp.flr	    \
	survexp.az	    \
	survexp.azr	    \
	survexp.wnc	


# I don't think I need to redo these in S5
#	$(MISS)              \
#	survsum              \


DOC  =  aareg	\
	coxph		\
	cluster         \
	coxph.detail    \
	coxph.object    \
	cox.zph         \
	dsurvreg	\
	is.ratetable    \
	frailty		\
	lines.survfit   \
	print.survfit   \
	print.summary.survfit	\
	print.summary.coxph	\
	pyears          \
	residuals.coxph \
	survdiff        \
	survexp         \
	survfit         \
	survfit.ci	\
	survfit.object  \
	survobrien	\
	plot.aareg      \
	plot.cox.zph    \
	plot.survfit    \
	predict.survreg \
	predict.coxph   \
	print.aareg     \
	residuals.survreg \
	strata          \
	survexp.fit     \
	survexp.us      \
	survreg         \
	survreg.object  \
        summary.aareg   \
	Surv            \
	summary.coxph   \
	summary.survfit \
	survConcordance \
	survsum         \
	untangle.specials 

#Missing value functions
MISS =  na.omit       \
	naresid.omit  \
	naprint.omit  \
	naresid       \
	naprint


# Functions that exist in S4, but I had to enhance them.
#  It is to be hoped that every one of these will go away, as the home
#  office responds to our complaints.
FIX     =  #subscript.model.matrix
FIXD    =



# Define the version 6 command on your system
SPLUS    = Splus
SHELP    = .Data/__Shelp
JHELP    = .Data/__Jhelp
HHELP    = .Data/__Hhelp

# Nothing from here on down should ever have to change

# These are the implicit rules that drive the makefile.  If you can
#   read them, consider yourself an expert!

.SUFFIXES:
.SUFFIXES: .o .c .c~ .r .r~ .f .f~ .h .h~ .s .s~ .sgm .sgm~

.sgml~.sgml:
	sccs get $<     

$(SHELP)/%.sgm: %.sgml
	cp -f $< $@

$(HHELP)/%.html: %.sgml
	$(SPLUS) sgml2html $< > $@

.Data/%: %.s
	$(SPLUS) < $<

HDOC = $(DOC:%=$(HHELP)/%.html)
SDOC = $(DOC:%=$(SHELP)/%.sgm)   #no use for this yet, but planned
INDEX= $(JHELP)/survivalMap.jhm
LFUN = $(FUN:%=.Data/%)
LRATE= $(RATES:%=.Data/%)
LDOC = $(DOC:%=%.sgml)            # this keep make from deleting sgml files

# "make all" to make a local version of S
#   test it out before you do the install
#
# Added .Data by Eric Lunde 2005-10-07 so that the build process can be done
# in one step (without having to do Splus CHAPTER -m).
all:    .Data doc fun  S.so

# only the docs, or only the functions
doc: $(LDOC) $(HDOC) $(SDOC) $(INDEX)
fun: $(LFUN) $(LRATE)

$(INDEX): $(HDOC)
	$(SHOME)/cmd/BUILD_JHELP


.Data :
	$(SPLUS) CHAPTER -m

temp:
	echo $(LDOC)

temp2: 
	echo $(LFUN)

include $(SHOME)/library/S_FLAGS

S.so:  $(OBJ) 
	$(SPLUS) LIBRARY S.so  $(OBJ) 

coxfit2.o agsurv3.o dmatrix.o pyears1.o pyears2.o pyears3.o: survS.h
coxfit2.o: survproto.h

# The year 2000 extrapolation depends on a data set, which can disappear
#   as soon as the rate tables are set up
.Data/survexp.us .Data/survexp.usr: survexp2000.sdump
.Data/survexp.mn .Data/survexp.mnwhite .Data/survexp.wnc: survexp2000.sdump
.Data/survexp.fl .Data/survexp.flr: survexp2000.sdump
.Data/survexp.az .Data/survexp.azr: survexp2000.sdump

# Commented out by Eric Lunde 2005-10-07.  sccs command is no longer needed
#survexp2000.sdump: SCCS/s.survexp2000.sdump
#	sccs get survexp2000.sdump




