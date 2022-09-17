##特征筛选######
compositontransition=function(aminoseq1)
{
x=aminoseq1
xamino=strsplit(x,"")
xamino1=xamino[[1]]
aminocompositon=table(xamino1)/length(xamino1)
name=names(aminocompositon)
tdata=as.numeric(aminocompositon)
oname=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
data=rep(0,20)
for(i in 1:length(name))
{j=str_which(oname,name[i])
 data[j]=tdata[i]
}
return(data)
}


hydoptransitondata=function(aminoseq1)
{
x=aminoseq1
#氨基酸序列被字母P,N,H代替
x1=gsub("G","Z",x)
x2=gsub("A","Z",x1)
x3=gsub("S","Z",x2)
x4=gsub("T","Z",x3)
x5=gsub("P","Z",x4)
x6=gsub("H","Z",x5)
x7=gsub("Y","Z",x6)

x8=gsub("R","P",x7)
x9=gsub("K","P",x8)
x10=gsub("E","P",x9)
x11=gsub("Q","P",x10)
x12=gsub("N","P",x11)
x13=gsub("D","P",x12)

x14=gsub("C","H",x13)
x15=gsub("V","H",x14)
x16=gsub("L","H",x15)
x17=gsub("I","H",x16)
x18=gsub("M","H",x17)
x19=gsub("F","H",x18)
x20=gsub("W","H",x19)
x211=gsub("Z","N",x20)

x212=strsplit(x211,"")
x21=x212[[1]]
tlength=length(x21)


#计算分布distribition
p1=grep("P",x21)#提取P的位置
plength=length(p1)#计算P的个数
if(plength!=0)
{
ap=floor(0.25*plength)
if(ap==0)ap=1
bp=floor(0.50*plength)
if(bp==0)bp=1
cp=floor(0.75*plength)
if(cp==0)cp=1
dp1=p1[c(1,ap,bp,cp,plength)]/tlength
}else dp1=c(0,0,0,0,0)

N1=grep("N",x21)
Nlength=length(N1)
if(Nlength!=0)
{
aN=floor(0.25*Nlength)
if(aN==0)aN=1
bN=floor(0.50*Nlength)
if(bN==0)bN=1
cN=floor(0.75*Nlength)
if(cN==0)cN=1
dN1=N1[c(1,aN,bN,cN,Nlength)]/tlength
}else dN1=c(0,0,0,0,0)

H1=grep("H",x21)
Hlength=length(H1)
if(Hlength!=0)
{
aH=floor(0.25*Hlength)
if(aH==0)aH=1
bH=floor(0.50*Hlength)
if(bH==0)bH=1
cH=floor(0.75*Hlength)
if(cH==0)cH=1
dH1=H1[c(1,aH,bH,cH,Hlength)]/tlength
}else dH1=c(0,0,0,0,0)
alld=c(dp1,dN1,dH1)
#计算组成composition
cp1=length(p1)/tlength
cN1=length(N1)/tlength
cH1=length(H1)/tlength
cdata=c(cp1,cN1,cH1)
cname=c("P","N","H")
allc=cdata
#计算转换数据transition
x22=gsub("P","-1",x21)
x23=gsub("N","0",x22)
x24=gsub("H","1",x23)
x25=as.numeric(x24)
x26=x25[-1]
x27=x25[-tlength]
tx=grep("TRUE",x26!=x27)
tx1=x25[tx]+x25[tx+1]
tpn=length(grep("-1",tx1))/tlength
tph=length(grep("0",tx1))/tlength
tnh=length(grep("1",tx1))/tlength
allt=c(tpn,tph,tnh)
data=c(allc,allt,alld)
return(data)
}

vandertransition=function(aminoseq1)
{
x=aminoseq1
#氨基酸序列被字母P,N,H代替
x1=gsub("G","P",x)
x2=gsub("A","P",x1)
x3=gsub("S","P",x2)
x4=gsub("C","P",x3)
x5=gsub("T","P",x4)
x6=gsub("P","P",x5)
x7=gsub("D","P",x6)

x8=gsub("N","N",x7)
x9=gsub("V","N",x8)
x10=gsub("E","N",x9)
x11=gsub("Q","N",x10)
x12=gsub("I","N",x11)
x13=gsub("L","N",x12)

x14=gsub("M","H",x13)
x15=gsub("H","H",x14)
x16=gsub("K","H",x15)
x17=gsub("F","H",x16)
x18=gsub("R","H",x17)
x19=gsub("Y","H",x18)
x20=gsub("W","H",x19)

#计算分布distribition
x212=strsplit(x20,"")
x21=x212[[1]]
tlength=length(x21)
p1=grep("P",x21)
plength=length(p1)
if(plength!=0)
{
ap=floor(0.25*plength)
if(ap==0)ap=1
bp=floor(0.50*plength)
if(bp==0)bp=1
cp=floor(0.75*plength)
if(cp==0)cp=1
dp1=p1[c(1,ap,bp,cp,plength)]/tlength
}else dp1=c(0,0,0,0,0)

N1=grep("N",x21)
Nlength=length(N1)
if(Nlength!=0)
{
aN=floor(0.25*Nlength)
if(aN==0)aN=1
bN=floor(0.50*Nlength)
if(bN==0)bN=1
cN=floor(0.75*Nlength)
if(cN==0)cN=1
dN1=N1[c(1,aN,bN,cN,Nlength)]/tlength
}else dN1=c(0,0,0,0,0)

H1=grep("H",x21)
Hlength=length(H1)
if(Hlength!=0){
aH=floor(0.25*Hlength)
if(aH==0)aH=1
bH=floor(0.50*Hlength)
if(bH==0)bH=1
cH=floor(0.75*Hlength)
if(cH==0)cH=1
dH1=H1[c(1,aH,bH,cH,Hlength)]/tlength
}else dH1=c(0,0,0,0,0)
alld=c(dp1,dN1,dH1)
#计算组成composition
cp1=length(p1)/tlength
cN1=length(N1)/tlength
cH1=length(H1)/tlength
cdata=c(cp1,cN1,cH1)
cname=c("P","N","H")
allc=cdata
#计算转换数据transition
x22=gsub("P","-1",x21)
x23=gsub("N","0",x22)
x24=gsub("H","1",x23)
x25=as.numeric(x24)
x26=x25[-1]
x27=x25[-tlength]
tx=grep("TRUE",x26!=x27)
tx1=x25[tx]+x25[tx+1]
tpn=length(grep("-1",tx1))/tlength
tph=length(grep("0",tx1))/tlength
tnh=length(grep("1",tx1))/tlength
allt=c(tpn,tph,tnh)
data=c(allc,allt,alld)
return(data)
}


polaritytransiton=function(aminoseq1)
{
x=aminoseq1
#氨基酸序列被字母P,N,H代替
x1=gsub("L","P",x)
x2=gsub("I","P",x1)
x3=gsub("F","P",x2)
x4=gsub("W","P",x3)
x5=gsub("C","P",x4)
x6=gsub("M","P",x5)
x7=gsub("Y","P",x6)
x8=gsub("V","P",x7)

x9=gsub("P","N",x8)
x10=gsub("A","N",x9)
x11=gsub("T","N",x10)
x12=gsub("G","N",x11)
x13=gsub("S","N",x12)

x14=gsub("H","H",x13)
x15=gsub("Q","H",x14)
x16=gsub("R","H",x15)
x17=gsub("K","H",x16)
x18=gsub("N","H",x17)
x19=gsub("E","H",x18)
x20=gsub("D","H",x19)

#计算分布distribition
x212=strsplit(x20,"")
x21=x212[[1]]
tlength=length(x21)
p1=grep("P",x21)
plength=length(p1)
if(plength!=0)
{
ap=floor(0.25*plength)
if(ap==0)ap=1
bp=floor(0.50*plength)
cp=floor(0.75*plength)
dp1=p1[c(1,ap,bp,cp,plength)]/tlength
}else dp1=c(0,0,0,0,0)

N1=grep("N",x21)
Nlength=length(N1)
if(Nlength!=0)
{
aN=floor(0.25*Nlength)
if(aN==0)aN=1
bN=floor(0.50*Nlength)
cN=floor(0.75*Nlength)
dN1=N1[c(1,aN,bN,cN,Nlength)]/tlength
}else dN1=c(0,0,0,0,0)

H1=grep("H",x21)
Hlength=length(H1)
if(Hlength!=0)
{
aH=floor(0.25*Hlength)
if(aH==0)aH=1
bH=floor(0.50*Hlength)
cH=floor(0.75*Hlength)
dH1=H1[c(1,aH,bH,cH,Hlength)]/tlength
}else dH1=c(0,0,0,0,0)
alld=c(dp1,dN1,dH1)
#计算组成composition
cp1=length(p1)/tlength
cN1=length(N1)/tlength
cH1=length(H1)/tlength
cdata=c(cp1,cN1,cH1)
cname=c("P","N","H")
allc=cdata
#计算转换数据transition
x22=gsub("P","-1",x21)
x23=gsub("N","0",x22)
x24=gsub("H","1",x23)
x25=as.numeric(x24)
x26=x25[-1]
x27=x25[-tlength]
tx=grep("TRUE",x26!=x27)
tx1=x25[tx]+x25[tx+1]
tpn=length(grep("-1",tx1))/tlength
tph=length(grep("0",tx1))/tlength
tnh=length(grep("1",tx1))/tlength
allt=c(tpn,tph,tnh)
data=c(allc,allt,alld)
return(data)
}


polarizabilitytransiton=function(aminoseq1)
{
x=aminoseq1
#氨基酸序列被字母P,N,H代替
x1=gsub("G","P",x)
x2=gsub("A","P",x1)
x3=gsub("S","P",x2)
x4=gsub("T","P",x3)
x5=gsub("D","P",x4)

x6=gsub("C","N",x5)
x7=gsub("P","N",x6)
x8=gsub("N","N",x7)
x9=gsub("V","N",x8)
x10=gsub("E","N",x9)
x11=gsub("Q","N",x10)
x12=gsub("I","N",x11)
x13=gsub("L","N",x12)

x14=gsub("K","H",x13)
x15=gsub("M","H",x14)
x16=gsub("H","H",x15)
x17=gsub("F","H",x16)
x18=gsub("R","H",x17)
x19=gsub("Y","H",x18)
x20=gsub("W","H",x19)


#计算分布distribition
x212=strsplit(x20,"")
x21=x212[[1]]
tlength=length(x21)
p1=grep("P",x21)
plength=length(p1)
if(plength!=0)
{
ap=floor(0.25*plength)
if(ap==0)ap=1
bp=floor(0.50*plength)
if(bp==0)bp=1
cp=floor(0.75*plength)
if(cp==0)cp=1
dp1=p1[c(1,ap,bp,cp,plength)]/tlength
}else dp1=c(0,0,0,0,0)

N1=grep("N",x21)
Nlength=length(N1)
if(Nlength!=0)
{
aN=floor(0.25*Nlength)
if(aN==0)aN=1
bN=floor(0.50*Nlength)
if(bN==0)bN=1
cN=floor(0.75*Nlength)
if(cN==0)cN=1
dN1=N1[c(1,aN,bN,cN,Nlength)]/tlength
}else dN1=c(0,0,0,0,0)

H1=grep("H",x21)
Hlength=length(H1)
if(Hlength!=0)
{
aH=floor(0.25*Hlength)
if(aH==0)aH=1
bH=floor(0.50*Hlength)
if(bH==0)bH=1
cH=floor(0.75*Hlength)
if(cH==0)cH=1
dH1=H1[c(1,aH,bH,cH,Hlength)]/tlength
}else dH1=c(0,0,0,0,0)
alld=c(dp1,dN1,dH1)
#计算组成composition
cp1=length(p1)/tlength
cN1=length(N1)/tlength
cH1=length(H1)/tlength
cdata=c(cp1,cN1,cH1)
cname=c("P","N","H")
allc=cdata
#计算转换数据transition
x22=gsub("P","-1",x21)
x23=gsub("N","0",x22)
x24=gsub("H","1",x23)
x25=as.numeric(x24)
x26=x25[-1]
x27=x25[-tlength]
tx=grep("TRUE",x26!=x27)
tx1=x25[tx]+x25[tx+1]
tpn=length(grep("-1",tx1))/tlength
tph=length(grep("0",tx1))/tlength
tnh=length(grep("1",tx1))/tlength
allt=c(tpn,tph,tnh)
data=c(allc,allt,alld)
return(data)
}

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
matrix84=function(aa){
train.matrix=matrix(0,nrow=length(aa),ncol=84)
for(i in 1:length(aa))
{aminoseq1=aa[i]
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
spls.train.feature=plsr(train.class2~train.matrix,ncomp=10,scale=TRUE)
test.predict=predict(spls.train.feature,test.matrix,dist="max.dist")
prediction=test.predict$class$max.dist[,2]
confusion.mat=get.confusion_matrix(truth=test.class2,predict=prediction)
accur.ber=1-get.BER(confusion.mat)
three.pre=as.numeric(confusion.mat[row(confusion.mat)==col(confusion.mat)]/apply(confusion.mat,1,sum))
accur.all=sum(confusion.mat[row(confusion.mat)==col(confusion.mat)])/sum(confusion.mat)
accuracy=c(accur.ber,accur.all,three.pre)
return(accuracy)
}

pls.train=function(train.matrix,train.class2,n)
{
pearson=c()
spearman=c()
spls.train.feature=plsr(train.class2~train.matrix,ncomp=n,scale=TRUE)
rmse=as.numeric(RMSEP(spls.train.feature)$val)
nrmse=rmse/sd(train.class2)
r2=as.numeric(R2(spls.train.feature)$val)
for(i in 1:n)
{
pearson[i]=cor(spls.train.feature[["fitted.values"]][,,i],train.class2)
spearman[i]=cor(spls.train.feature[["fitted.values"]][,,i],train.class2,method = "spearman")
}
out=list(rmse=rmse,nrmse=nrmse,r2=r2,pearson=pearson,spearman=spearman)
return(out)
}


library(ggplot2)       
library(RColorBrewer)
library(MASS)
library(dplyr)
library(stringr)
library(mixOmics)
library(caTools)
library(scorecard)
library(corrplot)
library(pls)

data.mrmr=read.table("F:/zhulin/1datasz1/RNA/mrmrorderfeature.csv",header=TRUE,sep=",")
data.valide=read.table("F:/zhulin/1datasz1/RNA/plsvalidedataset.csv",header=TRUE,sep=",")
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



########################################################
#########################  mrmr-order102-trian ###############################
train=data
train.aa=train$AA
train.xname=train$Entryname                   
train.class2=log10(train$fold)
################################
################################
train.matrix.all=cbind(matrix420(train.aa),matrix84(train.aa))
rownames(train.matrix.all)=train.xname
##dim(train.matrix.all)#查看矩阵的维度
##############################################################
###############特征###########
mrmr.name=data.mrmr$featurename
site=match(mrmr.name,colnames(train.matrix.all))
train.matrix182=train.matrix.all[,site]
pls.train(train.matrix182,train.class2,10)
###################mrmr-order添加#####################
mrmr.sequence=c(1:length(mrmr.name))

n=9
com=c(1:n)
outdata.rmse=data.frame(featurename=0,intercept=0,comps1=0,comps2=0,comps3=0,comps4=0,comps5=0,comps6=0,comps7=0,comps8=0,comps9=0,comps10=0)
outdata.nrmse=data.frame(featurename=0,intercept=0,comps1=0,comps2=0,comps3=0,comps4=0,comps5=0,comps6=0,comps7=0,comps8=0,comps9=0,comps10=0)
outdata.r2=data.frame(featurename=0,intercept=0,comps1=0,comps2=0,comps3=0,comps4=0,comps5=0,comps6=0,comps7=0,comps8=0,comps9=0,comps10=0)
outdata.pearson=data.frame(featurename=0,comps1=0,comps2=0,comps3=0,comps4=0,comps5=0,comps6=0,comps7=0,comps8=0,comps9=0,comps10=0)
outdata.spearman=data.frame(featurename=0,comps1=0,comps2=0,comps3=0,comps4=0,comps5=0,comps6=0,comps7=0,comps8=0,comps9=0,comps10=0)

for(i in 10:length(mrmr.sequence))
{
com=c(com,i)
a=match(mrmr.name[com],colnames(train.matrix.all))
train.matrix=train.matrix.all[,a]
one.out=pls.train(train.matrix,train.class2,n)
rmse=one.out$rmse
outdata.rmse=rbind(outdata.rmse,c(mrmr.name[i],rmse))
nrmse=one.out$nrmse
outdata.nrmse=rbind(outdata.nrmse,c(mrmr.name[i],nrmse))
r2=one.out$r2
outdata.r2=rbind(outdata.r2,c(mrmr.name[i],r2))
pearson=one.out$pearson
outdata.pearson=rbind(outdata.pearson,c(mrmr.name[i],pearson))
spearman=one.out$spearman
outdata.spearman=rbind(outdata.spearman,c(mrmr.name[i],spearman))
}

write.table(outdata.rmse,file="F:/zhulin/1datasz1/RNA/pls3-mrmr-order-rmse.csv",sep=",",row.names=F,col.names=T,quote=F)
write.table(outdata.nrmse,file="F:/zhulin/1datasz1/RNA/pls3-mrmr-order-nrmse.csv",sep=",",row.names=F,col.names=T,quote=F)
write.table(outdata.r2,file="F:/zhulin/1datasz1/RNA/pls3-mrmr-order-r2.csv",sep=",",row.names=F,col.names=T,quote=F)
write.table(outdata.pearson,file="F:/zhulin/1datasz1/RNA/pls3-mrmr-order-pearson.csv",sep=",",row.names=F,col.names=T,quote=F)
write.table(outdata.spearman,file="F:/zhulin/1datasz1/RNA/pls3-mrmr-order-spearman.csv",sep=",",row.names=F,col.names=T,quote=F)
