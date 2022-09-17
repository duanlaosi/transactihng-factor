library(ggplot2)       
library(RColorBrewer)
library(MASS)
library(dplyr)
library(stringr)
library(mixOmics)
library(caTools)
library(scorecard)
library(corrplot)
##读取数据
data.valide=read.table("F:/zhulin/1datasz1/RNA/validedataset.csv",header=TRUE,sep=",")
s=c(5,45,46,57)
data=data.valide[-s,]
feature.name182=c(paste(c("x"),c(1:84),sep=""))
#################420个特征######
acid=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
dacid=c("AA","CC","DD","EE","FF","GG","HH","II","KK","LL","MM","NN","PP","QQ","RR","SS","TT","VV","WW","YY")
amino.acid=combn(acid,2)
acid.1=paste(c("^"),acid,c("$"),sep="")
dacid.1=paste(c("^"),dacid,c("$"),sep="")
dual.s=c(paste(c("^"),amino.acid[2,],amino.acid[1,],c("$"),sep=""),paste(c("^"),amino.acid[1,],amino.acid[2,],c("$"),sep=""),dacid.1,acid.1)
dual.sing=c(paste(amino.acid[2,],amino.acid[1,],sep=""),paste(amino.acid[1,],amino.acid[2,],sep=""),dacid,acid)
###############
##########123个特征的位置和名字
featuredata=read.table("F:/zhulin/1datasz1/RNA/features.csv",header=TRUE,sep=",")
featurename=featuredata[,1][-1]
site=match(featurename,dual.sing)

train.aa=data$protein.sequence
train.xname=data$Entry.name                   
train.class2=data$class

################################
################################
#建立数据框
train.matrix=matrix(0,nrow=length(train.aa),ncol=84)
#调用转换代码
for(i in 1:length(train.aa))
{aminoseq1=train.aa[i]
train.matrix[i,c(1:21)]=hydoptransitondata(aminoseq1)
train.matrix[i,c(22:42)]=vandertransition(aminoseq1)
train.matrix[i,c(43:63)]=polaritytransiton(aminoseq1)
train.matrix[i,c(64:84)]=polarizabilitytransiton(aminoseq1)
}

rownames(train.matrix)=train.class2
colnames(train.matrix)=feature.name182

##############sliding windows 计算频率############
############train 矩阵##########
train.420matrix=matrix(0,length(train.aa),420)
for(i in 1:length(train.aa))
{
train.420matrix[i,]=find(train.aa[i])
}
rownames(train.420matrix)=train.class2
colnames(train.420matrix)=dual.sing

##直接使用123个特征
train.matrix123=train.420matrix[,site]
#######合并数据矩阵#########
train.102matrix182=cbind(train.matrix123,train.matrix)

######特征选择##############
feature.selection=founction(train.102matrix){
train.102matrix182.1=train.102matrix182[,-which(colSums(train.102matrix182)<=0)]
train.102matrix182.1=train.102matrix182.1[,-which(colSums(train.102matrix182.1)>=nrow(train.102matrix182))]
corrr=cor(train.102matrix182.1)
corrplot(corrr,tl.cex=0.5)
}
write.table(train.102matrix182.1,file="F:/zhulin/1datasz1/RNA/featurematrix.csv",sep=",",row.names=T,col.names=T,quote=F)

