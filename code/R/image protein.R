library(ggplot2)       
library(RColorBrewer)
library(MASS)
library(dplyr)
library(stringr)
library(mixOmics)
library(caTools)
library(scorecard)
library(corrplot)
library(grid)
data.valide=read.table("F:/zhulin/1datasz1/RNA/validedataset.csv",header=TRUE,sep=",")
s=c(5,45,46,57)
protein.name=data$Entry.name
protein.sequence=data$protein.sequence
protein.class=data$class
try.1=protein.sequence[1]

amino.acid=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
image.length=36
try.1.length=36
grid.newpage()
colvec=c('FFCCFF','FFCC99','CCCCFF','CCCCCC','CCCC33','99CCFF','99CC99','990099','FF0000',
         '6666FF','666600','3366FF','00CCFF','0066FF','330000','FFFF33','FF00FF','BBBBBB',
         'FF0033','CC6600')
vp1=viewport(x=0.1,y=0.1,w=0.8,h=0.8,just=c("left","bottom"),name="vp1")
grid.rect(x=20,y=1,height=1,width=1,hjust=0,vjust=0,vp=vp1,gp=gpar(col=1,fill=as.character(colvec)))
grid.show.viewport(vp1)
grid(nx=36,NA,col="lightgray")
find=function(sequence){
k=try.1
str.one=str_split(k,"")[[1]]
kk=c(paste(str.one[-length(str.one)],str.one[-1],sep=""),str.one)
single=rep(0,420)
for(i in 1:420){
single[i]=length(grep(kk,pattern=dual.s[i]))
}
single=single/sum(single)
return(single)
}

plot(1:3)
grid(NA, 5, lwd = 2) # grid only in y-direction



