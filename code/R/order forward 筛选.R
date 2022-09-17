data.mrmr=read.table("F:/zhulin/1datasz1/RNA/mrmrorderfeature.csv",header=TRUE,sep=",")
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

####################################
###############建立420feature的矩阵#####################
matrix420=function(aa)
{
train.420matrix=matrix(0,length(aa),420)
for(i in 1:length(aa))
{
train.420matrix[i,]=find(aa[i])
}
colnames(train.420matrix)=dual.sing
return(train.420matrix)
}
#################################
#######################建立理化性质的矩阵######################
matrix84=function(aa)
{
train.matrix=matrix(0,nrow=length(aa),ncol=84)
for(i in 1:length(aa))
{aminoseq1=train.aa[i]
train.matrix[i,c(1:21)]=hydoptransitondata(aminoseq1)
train.matrix[i,c(22:42)]=vandertransition(aminoseq1)
train.matrix[i,c(43:63)]=polaritytransiton(aminoseq1)
train.matrix[i,c(64:84)]=polarizabilitytransiton(aminoseq1)
} 
colnames(train.matrix)=feature.name182
return(train.matrix)
}
######################剔除mrmr order feature 中的20个compositon feature#########
mrmr.order=function(data.mrmr)
{
mrmr.featurename=data.mrmr$featurename
b=c()
for(i in 1:length(acid.1)){
a=grep(acid.1[i],mrmr.featurename)
b=c(b,a)
}
return(data.mrmr[-b,])
}
########################################################

#############pls###############
pls=function(train.matrix,test.matrix,train.class2,test.class2)
{
spls.train.feature=splsda(train.matrix,train.class2,ncomp=3,scale=TRUE)
test.predict=predict(spls.train.feature,test.matrix,dist="max.dist")
prediction=test.predict$class$max.dist[,2]
confusion.mat=get.confusion_matrix(truth=test.class2,predict=prediction)
accur.ber=1-get.BER(confusion.mat)
accur.all=sum(confusion.mat[row(confusion.mat)==col(confusion.mat)])/sum(confusion.mat)
accuracy=c(accur.ber,accur.all)
return(accuracy)
}

########################################################
########################################################
test.accurarcy=function(n,data)
{
nr=nrow(data)
trainindex=sample(1:nr,n*nr)
train=data[trainindex,]
test=data[-trainindex,]
train.aa=train$protein.sequence
train.xname=train$Entry.name                   
train.class2=train$class
#############################
test.aa=test$protein.sequence  
test.xname=test$Entry.name                     
test.class2=test$class
################################
################################
train.matrix.all=cbind(matrix420(train.aa),matrix84(train.aa))
rownames(train.matrix.all)=train.xname
test.matrix.all=cbind(matrix420(test.aa),matrix84(test.aa))
rownames(test.matrix.all)=test.xname
##############################################################
###############特征###########
mrmr.name=mrmr.order(data.mrmr)$featurename
com=c(401,402)
outdata=data.frame(featurename=0,accuracy.BER=0,accuracy.all=0)
for(i in 403:420)
{
com=c(com,i)
train.matrix=train.matrix.all[,com]
test.matrix=test.matrix.all[,com]
mab=c(colnames(train.matrix.all)[i],pls(train.matrix,test.matrix,train.class2,test.class2))
outdata=rbind(outdata,mab)
}
for(i in 1:length(mrmr.name))
{
a=grep(paste("^",mrmr.name[i],"$",sep=""),colnames(train.matrix.all))
com=c(com,a)
train.matrix=train.matrix.all[,com]
test.matrix=test.matrix.all[,com]
mab=c(mrmr.name[i],pls(train.matrix,test.matrix,train.class2,test.class2))
outdata=rbind(outdata,mab)
}
return(outdata)
}

###########调试重复测试20遍#########
re=20
n=0.8
BER.value=matrix(0,nrow=182,ncol=re)
all.value=matrix(0,nrow=182,ncol=re)
for(i in 1:re)
{outdata=test.accurarcy(n,data)
 BER.value[,i]=as.numeric(outdata$accuracy.BER)
 all.value[,i]=as.numeric(outdata$accuracy.all)
BER=apply(BER.value,1,mean)
all=apply(all.value,1,mean)
}
out.value=data.frame(featurename=outdata$featurename,predict.accuracy.BER=BER,predict.accuracy.all=all)
out.value=rbind(c("A","*","*"),out.value)
out.value[2,]=c("C","*","*")
write.table(out.value,file="F:/zhulin/1datasz1/RNA/predict.accuracy.pls.csv",sep=",",row.names=F,col.names=T,quote=F)


