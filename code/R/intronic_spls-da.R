library(ggplot2)       
library(RColorBrewer)
library(MASS)
library(dplyr)
library(stringr)
library(mixOmics)

data=read.table("F:/zhulin/1datasz1/RNA/trainingdataset.csv",header=TRUE,sep=",")
data_1=read.table("F:/zhulin/1datasz1/RNA/Intronic_Splicing_Model_(GUR).csv",header=TRUE,sep=",")
site_no=c(23,30,57,60,68)
name=data_1$Name.New[-site_no]
site=match(name,data$Name)
sequence=data$AA[site]
class=data_1$class[-site_no]

##
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
SD.matrix[i,]=find(sequence[i])
}
rownames(SD.matrix)=name
colnames(SD.matrix)=dual.sing
##


##直接使用文中给的123个特征
featuredata=read.table("F:/zhulin/1datasz1/RNA/features.csv",header=TRUE,sep=",")
featurename=featuredata[,1][-1]
site_feature=match(featurename,dual.sing)
matrix=SD.matrix[,site_feature]
###figure4B
plasdarna=splsda(matrix,class,ncomp=2)
performance=perf(plasdarna,validation = "loo")
performance$error.rate
performance$error.rate.class

plotIndiv(plasdarna,comp=c(1,2),group=class,col=c("blue","grey","red"),lty=2,ind.names=TRUE,ellipse=TRUE,legend=TRUE,X.label='Component 1',
Y.label='Component 2',title='Splicing factors tested',plot.title = element_text(hjust = 0))

head(plasdarna$loadings)##获得每个特征的权重




represention=function(sequence){
single=rep(0,420)
for(i in 1:420){
single[i]=str_count(sequence,pattern=dual.sing[i])
}
single=single/sum(single)
return(single)
}
X=matrix(0,63,420)
for(i in 1:63){
X[i,]=represention(sequence[i])
}
rownames(X)=name
colnames(X)=dual.sing
matrix.scale=X[,site_feature]
plasdarna=plsda(matrix.scale,class,ncomp=3)
plotIndiv(plasdarna,group=class,col=c("blue","grey","red"),lty=2,ind.names=TRUE,ellipse=TRUE,legend=TRUE,X.label='Component 1',
Y.label='Component 2',title='Splicing factors tested',plot.title = element_text(hjust = 0))
head(plasdarna$loadings)
plotVar(plasdarna, comp =c(1,2), title = 'Loadings on comp 1')

