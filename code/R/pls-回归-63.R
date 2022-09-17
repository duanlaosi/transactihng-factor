library(ggplot2)       
library(RColorBrewer)
library(MASS)
library(dplyr)
library(stringr)
library(mixOmics)
library(pls)
##
data1=read.table("F:/zhulin/1datasz1/RNA/trainingdataset.csv",header=TRUE,sep=",")
head(data1)
s=c(5,45,46,57)
data=data1[-s,]
aa=data$AA  
name=data$Name                    
class2=log10(data$Average)
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

##直接使用文中给的123个特征
featuredata=read.table("F:/zhulin/1datasz1/RNA/features.csv",header=TRUE,sep=",")
featurename=featuredata[,1][-1]
site=match(featurename,dual.sing)
matrix=SD.matrix[,site]

pls.63=plsr(class2~matrix,scale = TRUE, ncomp =10)
plot(pls.63,line=TRUE)
plot(R2(pls.63),legendpos = "topright")

plot(RMSEP(pls.63), legendpos = "topright")
x=RMSEP(pls.63)
rmse=as.numeric(x$val)
nrmse=as.numeric(x$val)/sd(class2)
plot(nrmse,type='l',main='NRMSE')
for(i in 1:10)
{
pearson[i]=cor(pls.63[["fitted.values"]][,,i],class2)
spearman[i]=cor(pls.63[["fitted.values"]][,,i],class2,method = "spearman")
}
pearson
spearman