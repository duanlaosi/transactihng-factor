library(ggplot2)       
library(RColorBrewer)
library(MASS)
library(dplyr)
library(stringr)
library(mixOmics)
library(pls)

##63个数据
data=read.table("F:/zhulin/1datasz1/RNA/trainingdataset.csv",header=TRUE,sep=",")
####102个数据
data.valide=read.table("F:/zhulin/1datasz1/RNA/plsvalidedataset.csv",header=TRUE,sep=",")
s=c(5,45,46,57)
data=data.valide[-s,]

ready=function(data)
{
aa=data$AA  
name=data$Entryname                    
class2=log10(data$fold)
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
SD.matrix=matrix(0,length(aa),420)
for(i in 1:length(aa)){
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
plasdarna=plsr(class2~matrix,ncomp=10,scale=TRUE)
plot(plasdarna,line=TRUE)
plot(R2(plasdarna),legendpos = "topright")

plot(RMSEP(plasdarna), legendpos = "topright")
x=RMSEP(plasdarna)
rmse=as.numeric(x$val)
nrmse=as.numeric(x$val)/sd(class2)
plot(nrmse,type='l',main='NRMSE')
pearson=c()
spearman=c()
for(i in 1:10)
{
pearson[i]=cor(plasdarna[["fitted.values"]][,,i],class2)
spearman[i]=cor(plasdarna[["fitted.values"]][,,i],class2,method = "spearman")
}
pearson
spearman

table=data.frame(Fold=class2,predict=plasdarna[["fitted.values"]][,,10])
plot(table)
write.table(table,file="F:/zhulin/1datasz1/RNA/91-123.csv",sep=",",row.names=F,col.names=T,quote=F)

}
ready(data)






