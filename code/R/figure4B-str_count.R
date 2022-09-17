##figure4B，使用str_count计算频率
data=read.table("F:/zhulin/1datasz1/RNA/trainingdatasets.csv",header=TRUE,sep=",")
aa=data$AA  
xname=data$Name.1                    
s=c(5,45,46,57)
name=xname[-s]
aa=aa[-s]
class2=data$class2
class2=class2[-s]
acid=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
dacid=c("AA","CC","DD","EE","FF","GG","HH","II","KK","LL","MM","NN","PP","QQ","RR","SS","TT","VV","WW","YY")
amino.acid=combn(acid,2)
dual.sing=c(paste(amino.acid[2,],amino.acid[1,],sep=""),paste(amino.acid[1,],amino.acid[2,],sep=""),dacid,acid)

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
X[i,]=represention(aa[i])
}
rownames(X)=name
colnames(X)=dual.sing

featuredata=read.table("F:/zhulin/1datasz1/RNA/feature.csv",header=TRUE,sep=",")
featurename=featuredata[,1][-1]
site=match(featurename,dual.sing)
feature.name=dual.sing[site]
matrix=X[,site]
matrix.scale=scale(matrix,center=TRUE,scale=TRUE)
##主成分分析
pca.srbct = pca(matrix, ncomp = 3, center = TRUE, scale = TRUE)
plot(pca.srbct)
plotIndiv(pca.srbct,group=class2,col=c("blue","grey","red"),lty=2,ind.names=TRUE,ellipse=TRUE,legend=TRUE,X.label='Component 1',
Y.label='Component 2',title='Splicing factors tested',plot.title = element_text(hjust = 0))
##PLD-DA分析
plasdarna=plsda(matrix.scale,class2,ncomp=3)
plotIndiv(plasdarna,group=class2,col=c("blue","grey","red"),lty=2,ind.names=TRUE,ellipse=TRUE,legend=TRUE,X.label='Component 1',
Y.label='Component 2',title='Splicing factors tested',plot.title = element_text(hjust = 0))
head(plasdarna$loadings)
plotVar(plasdarna, comp =c(1,2), title = 'Loadings on comp 1')

background = background.predict(plasdarna, comp.predicted=2, dist = "max.dist") 
plotIndiv(plasdarna, comp = 1:2,
          group = class2, ind.names = TRUE, title = "Maximum distance",
          legend = TRUE,  background = background)

perf.plsda<- perf(plasdarna, validation = "Mfold", folds =10, 
                  progressBar = FALSE, auc = TRUE, nrepeat = 10)
plot(perf.plsda)
perf.plsda$error.rate
