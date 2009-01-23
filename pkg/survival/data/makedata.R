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
#save(pbc, file='pbc.rda')

pbcseq <- read.table('pbcseq.dat', header=F,
                     col.names=c('id', 'futime', 'status', 'trt', 'age',
                                 'sex', 'day', 'ascites', 'hepato',
                                 'spiders', 'edema', 'bili', 'chol',
                                 'albumin', 'alk.phos', 'ast', 'platelet',
                                 'protime', 'stage'),
                     na.strings='.')
pbcseq$age <- pbcseq$age/365.25
pbcseq$sex <- factor(pbcseq$sex, levels=0:1, labels=c("m", "f"))
#save(pbcseq, file='pbc.rda')

lung <- read.table("lung.dat", header=T, na.strings='.')

# Monoclongal gammopathy
mgus <- read.table("mgus.dat", header=T, sep=',')
mgus$sex <- factor(mgus$sex, 1:2, c("male", "female"))

# Make the Wei,Lin,Weissfeld style data set, with risks for death,
#  multiple myeloma, and other plasma cell disease as parallel risks
temp1 <- cbind(data.frame(id=mgus$id,
                    time = pmin(mgus$pctime, mgus$futime, na.rm=T),
                    status= 1*(mgus$pcdx=='MM'),
                    event='myeloma'), 
               mgus[,c('age', 'sex', 'alb','creat', 'hgb', 'mspike')])
temp2 <- cbind(data.frame(id=mgus$id,
                    time = pmin(mgus$pctime, mgus$futime, na.rm=T),
                    status= 1*(!is.na(match(mgus$pcdx, c('LP','MA','AM')))),
                    event='other'), 
               mgus[,c('age', 'sex', 'alb','creat', 'hgb', 'mspike')])
temp3 <- cbind(data.frame(id=mgus$id,
                          time= mgus$futime,
                          status = mgus$death,
                          event='death'),
               mgus[,c('age', 'sex', 'alb','creat', 'hgb', 'mspike')])

mgus2 <- rbind(temp1, temp2, temp3)
mgus2 <- mgus2[order(mgus2$id, mgus2$event),]                    
row.names(mgus2) <- NULL
rm(temp1, temp2, temp3)

leukemia <- read.table('leukemia.dat', header=T)
aml <- leukemia

#quit(save='yes')
