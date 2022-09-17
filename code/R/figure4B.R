library(ggplot2)       
library(RColorBrewer)
library(MASS)
library(dplyr)
library(stringr)
library(mixOmics)

##
data=read.table("F:/zhulin/1datasz1/RNA/trainingdataset.csv",header=TRUE,sep=",")
aa=data$AA  
xname=data$Name                   
s=c(5,45,46,57)
name=xname[-s]
aa=aa[-s]
class2=data$class
class2=class2[-s]
acid=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
dacid=c("AA","CC","DD","EE","FF","GG","HH","II","KK","LL","MM","NN","PP","QQ","RR","SS","TT","VV","WW","YY")
amino.acid=combn(acid,2)
acid.1=paste(c("^"),acid,c("$"),sep="")
dacid.1=paste(c("^"),dacid,c("$"),sep="")
dual.s=c(paste(c("^"),amino.acid[2,],amino.acid[1,],c("$"),sep=""),paste(c("^"),amino.acid[1,],amino.acid[2,],c("$"),sep=""),dacid.1,acid.1)
dual.sing=c(paste(amino.acid[2,],amino.acid[1,],sep=""),paste(amino.acid[1,],amino.acid[2,],sep=""),dacid,acid)

###
find=function(sequence){
k=sequence
str.one=str_split(k,"")[[1]]
kk=c(paste(str.one[-length(str.one)],str.one[-1],sep=""),str.one)
single=rep(0,420)
for(i in 1:420){
single[i]=length(grep(kk,pattern=dual.s[i]))
}
single=single/sum(single)
return(single)
}
SD.matrix=matrix(0,63,420)
for(i in 1:63){
SD.matrix[i,]=find(aa[i])
}
rownames(SD.matrix)=name
colnames(SD.matrix)=dual.sing
#length(nearZeroVar(SD.matrix,freqCut=25,uniqueCut=18)$Position)
#near.position=nearZeroVar(SD.matrix,freqCut=25,uniqueCut=17.5)$Position
#dual.sing[-near.position]
#SD.matrix[,-near.position]
#SD.matrix[,1]
###筛选过滤

##直接使用文中给的123个特征
featuredata=read.table("F:/zhulin/1datasz1/RNA/features.csv",header=TRUE,sep=",")
featurename=featuredata[,1][-1]
site=match(featurename,dual.sing)
matrix=SD.matrix[,site]
###figure4B
plasdarna=splsda(matrix,class2,ncomp=3)
performance=perf(plasdarna,validation = "loo")
performance$error.rate
performance$predict
plot(performance)
plotIndiv(plasdarna,group=class2,col=c("blue","grey","red"),lty=2,ind.names=TRUE,ellipse=TRUE,legend=TRUE,X.label='Component 1',
Y.label='Component 2',title='Splicing factors tested',plot.title = element_text(hjust = 0))

head(plasdarna$loadings)##获得每个特征的权重

##tune分析
list.keepX <- c(1:10, seq(20, 150, 10))
tunesplsda <- tune.splsda(matrix, class2, ncomp = 3, validation = 'Mfold', folds = 5, progressBar = TRUE, dist = 'max.dist', measure = "BER",
                          test.keepX = list.keepX, nrepeat = 10, cpus = 2)
error <- tunesplsda$error.rate 
ncomp <- tunesplsda$choice.ncomp$ncomp # optimal number of components based on t-tests
ncomp
select.keepX <- tunesplsda$choice.keepX[1:3]  # optimal number of variables to select
select.keepX
plot(tunesplsda, col = color.jet(3))
splsda.srbct <- splsda(matrix, class2, ncomp = 3, keepX = select.keepX) 
plotIndiv(splsda.srbct, comp = c(1,2),
          group =class2, ind.names = TRUE, 
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA on SRBCT, comp 1 & 2')

head(splsda.srbct$loadings)




