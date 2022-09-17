

#####################
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
################################
################################

library(ggplot2)       
library(RColorBrewer)
library(MASS)
library(dplyr)
library(stringr)
library(mixOmics)
library(caTools)
library(scorecard)
##读取数据
data=read.table("F:/zhulin/1datasz1/RNA/trainingdatasets.csv",header=TRUE,sep=",")
feature.name=c(paste(c("x"),c(1:104),sep=""))
s=c(5,45,46,57)
data=data[-s,]
re=20
n=0.8
pre.value=rep(0,re)
for(i in 1:re)
{pre.value[i]=test78(n,data)
}
pre.value
mean(pre.value)
sd(pre.value)

test78=function(n,data)
{
nr=nrow(data)
trainindex=sample(1:nr,n*nr)
train=data[trainindex,]
test=data[-trainindex,]
train.aa=train$AA  
train.xname=train$Name.1                    
train.class2=train$class2
##############################
test.aa=test$AA  
test.xname=test$Name.1                    
test.class2=test$class2
#建立train数据框
train.matrix=matrix(0,nrow=length(train.aa),ncol=104)
#调用转换代码
for(i in 1:length(train.aa))
{aminoseq1=train.aa[i]
train.matrix[i,c(1:20)]=compositontransition(aminoseq1)
train.matrix[i,c(21:41)]=hydoptransitondata(aminoseq1)
train.matrix[i,c(42:62)]=vandertransition(aminoseq1)
train.matrix[i,c(63:83)]=polaritytransiton(aminoseq1)
train.matrix[i,c(84:104)]=polarizabilitytransiton(aminoseq1)
}
#建立test数据框
test.matrix=matrix(0,nrow=length(test.aa),ncol=104)
#调用转换代码
for(i in 1:length(test.aa))
{aminoseq1=test.aa[i]
test.matrix[i,c(1:20)]=compositontransition(aminoseq1)
test.matrix[i,c(21:41)]=hydoptransitondata(aminoseq1)
test.matrix[i,c(42:62)]=vandertransition(aminoseq1)
test.matrix[i,c(63:83)]=polaritytransiton(aminoseq1)
test.matrix[i,c(84:104)]=polarizabilitytransiton(aminoseq1)
}
rownames(train.matrix)=train.xname
rownames(test.matrix)=test.xname
colnames(test.matrix)=feature.name
colnames(train.matrix)=feature.name
#######去处全部为0的特征变量#########
train.position=nearZeroVar(train.matrix,freqCut=25,uniqueCut=18)$Position
train.matrix78=train.matrix[,-train.position]
test.matrix78=test.matrix[,-train.position]

##PLD-DA分析
spls.train.78feature=splsda(train.matrix78,train.class2,ncomp=3,scale=TRUE)
#plotIndiv(spls.train.78feature,group=train.class2,col=c("blue","grey","red"),lty=2,ind.names=TRUE,ellipse=TRUE,legend=TRUE,X.label='Component 1',
#Y.label='Component 2',title='train-63-new-78features',plot.title = element_text(hjust = 0))
#perf.train.plsda<- perf(spls.train.78feature, validation = "Mfold", folds =3,progressBar = TRUE, auc = TRUE, nrepeat = 10)
#perf.train.plsda$error.rate
#perf.train.plsda$error.rate.class
#plot(perf.train.plsda)
#auc.plsda=auroc(spls.train.78feature, roc.comp=3)
test.predict=predict(spls.train.78feature,test.matrix78,dist="max.dist")
prediction=test.predict$class$max.dist[,2]
confusion.mat=get.confusion_matrix(truth=test.class2,predict=prediction)
return(1-get.BER(confusion.mat))
}














