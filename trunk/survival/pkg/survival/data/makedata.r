#
# This file creates various data sets from ascii sources
#   (Mostly from the Therneau and Grambsch book).
pbc <- read.table('pbc.dat', header=F, 
                  col.names=c('id', 'time', 'status', 'trt',  'age', 'sex',
                              'ascites',  'hepato',  'spiders',  'edema',
                              'bili',  'chol',  'albumin',  'copper', 
                              'alk.phos',  'ast',  'trig',  'platelet',
                              'protime',  'stage'),
                  na.strings='.')

pbc$age <- pbc$age/365.25
pbc$sex <- factor(pbc$sex, levels=0:1, labels=c("m", "f"))
save(pbc, file='pbc.rda')

pbcseq <- read.table('pbcseq.dat', header=F,
                     col.names=c('id', 'futime', 'status', 'trt', 'age',
                                 'sex', 'day', 'ascites', 'hepato',
                                 'spiders', 'edema', 'bili', 'chol',
                                 'albumin', 'alk.phos', 'ast', 'platelet',
                                 'protime', 'stage'),
                     na.strings='.')
pbcseq$age <- pbcseq$age/365.25
pbcseq$sex <- factor(pbcseq$sex, levels=0:1, labels=c("m", "f"))
save(pbcseq, file='pbc.rda')
