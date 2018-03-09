#########################################################
#########################################################
#Apply baysian network method from the 'bnlearn' package to the sleep 
#study data
#########################################################
#########################################################


library(mice)
library(bnlearn)
library(igraph)


#data cleaning
#save data in longitudinal form
setwd('/Users/Selene/Desktop/ebdbNet')
data=read.csv("chemodatabase_06252012.csv",header=TRUE)
data=data[data$group==1,]

ind_keep=c(1:3,12,23,231,97,135,154,147,160);names(data[,ind_keep])
data=data[,ind_keep]
dim(data);summary(data)
id=unique(data$sub)
N=length(id)
P=length(ind_keep)-3

data.long=as.data.frame(array(NA,dim=c(N*3,dim(data)[2])))
names(data.long)=names(data)
data.long[,1]=rep(id,each=3)
data.long[,2]=rep(0:2,N)
for(i in 1:N){
	data.long.i=data.long[data.long$sub==id[i],]
	data.i=data[data$sub==id[i],]
	data.long.i[,3]=data.i[1,3]
	for(p in 1:P){
		for(t in 0:2){
			if(length(data.i[data.i$phase==t,p+3]==1)){
				data.long.i[t+1,p+3]=data.i[data.i$phase==t,p+3]
			}	
		}	
	}
	data.long[data.long$sub==id[i],]=data.long.i
}	
for(i in 1:N){
	data.long[data.long$sub==id[i],4]=data.long[data.long$sub==id[i],4][1]
	data.long[data.long$sub==id[i],5]=data.long[data.long$sub==id[i],5][1]
	
}

save(data.long,file='data.long.rdata')

#initial analysis of each variable over different phases
#mle on time
library(nlme)
data.long$sub=as.character(data.long$sub)
data.long$phase=as.factor(data.long$phase)
lis=list()
mod=lme(compscore~phase,random=~1|sub,data.long,method='REML',na.action = na.exclude)
lis[[1]]=summary(mod)$tTable
mod=lme(cesdtotal~phase,random=~1|sub,data.long,method='REML',na.action = na.exclude)
lis[[2]]=summary(mod)$tTable
mod=lme(mfsitotal~phase,random=~1|sub,data.long,method='REML',na.action = na.exclude)
lis[[3]]=summary(mod)$tTable
mod=lme(factbtotal~phase,random=~1|sub,data.long,method='REML',na.action = na.exclude)
lis[[4]]=summary(mod)$tTable
mod=lme(psqitotal~phase,random=~1|sub,data.long,method='REML',na.action = na.exclude)
lis[[5]]=summary(mod)$tTable
mod=lme(fosqtotal~phase,random=~1|sub,data.long,method='REML',na.action = na.exclude)
lis[[6]]=summary(mod)$tTable
names(lis)=c('compscore','cesdtotal','mfsitotal','factbtotal','psqitotal','fosqtotal')
setwd('/Users/Selene/Desktop')
write.csv(lis,file='lis.csv')


#deal with missing data, either impute or delete
imp=FALSE
if(imp==FALSE){
	data.long.keep=rep(1,dim(data.long)[1])
	ind=1:3
	for(i in 1:N){
		if(sum(is.na(data.long[data.long$sub==id[i],(1+3):(P+3)]))>0){
			data.long.keep[ind]=0
		}
		ind=ind+3
	}
	data.long=data.long[data.long.keep==1,]
}
if(imp==TRUE){
	ind_imp=c(3,(1+3):(P+3));names(data.long[,ind_imp])
	data.imp0 <- mice(data.long[data.long$phase==0,ind_imp],m=20)
	data.imp1 <- mice(data.long[data.long$phase==1,ind_imp],m=20)
	data.imp2 <- mice(data.long[data.long$phase==2,ind_imp],m=20)
	data.long[data.long$phase==0,ind_imp]=complete(data.imp0)
	data.long[data.long$phase==1,ind_imp]=complete(data.imp1)
	data.long[data.long$phase==2,ind_imp]=complete(data.imp2)
}

#########################################################
#build a BN model for phase 0 to 1
#########################################################
d1=data.long[!data.long$phase==2,]
Cognition0=d1[d1$phase==0,'compscore']
Cognition1=d1[d1$phase==1,'compscore']
Mood0=d1[d1$phase==0,'cesdtotal']
Mood1=d1[d1$phase==1,'cesdtotal']
Fatigue0=d1[d1$phase==0,'mfsitotal']
Fatigue1=d1[d1$phase==1,'mfsitotal']
QoL0=d1[d1$phase==0,'factbtotal']
QoL1=d1[d1$phase==1,'factbtotal']
Sleep0=d1[d1$phase==0,'psqitotal']
Sleep1=d1[d1$phase==1,'psqitotal']
FOSQ0=d1[d1$phase==0,'fosqtotal']
FOSQ1=d1[d1$phase==1,'fosqtotal']
Age=d1[d1$phase==0,'cb_demage']
Education=d1[d1$phase==0,'cb_demeducation']

d1.s=cbind.data.frame(Cognition0,Cognition1,Mood0,Mood1, Fatigue0,Fatigue1,QoL0,QoL1,Sleep0,Sleep1,FOSQ0,FOSQ1,Age,Education)

#create a blacklist disallowing arrows going certain direction based on
#common sense such as phase constraint
m=(ncol(d1.s)-2)/2
bl1=cbind.data.frame(names(d1.s)[c(2,4,6,8,10,12)],rep("Cognition0",m))
names(bl1)=c('from','to')
bl2=cbind.data.frame(names(d1.s)[c(2,4,6,8,10,12)],rep("Mood0",m))
names(bl2)=c('from','to')
bl3=cbind.data.frame(names(d1.s)[c(2,4,6,8,10,12)],rep("Fatigue0",m))
names(bl3)=c('from','to')
bl4=cbind.data.frame(names(d1.s)[c(2,4,6,8,10,12)],rep("QoL0",m))
names(bl4)=c('from','to')
bl5=cbind.data.frame(names(d1.s)[c(2,4,6,8,10,12)],rep("Sleep0",m))
names(bl5)=c('from','to')
bl6=cbind.data.frame(names(d1.s)[c(2,4,6,8,10,12)],rep("FOSQ0",m))
names(bl6)=c('from','to')
bl7=cbind.data.frame(rep("Cognition0",m*2-2),names(d1.s)[-c(1,2,13,14)])
names(bl7)=c('from','to')
bl8=cbind.data.frame(rep("Cognition1",m-1),names(d1.s)[c(4,6,8,10,12)])
names(bl8)=c('from','to')
bl9=cbind.data.frame(names(d1.s)[-13],rep("Age",ncol(d1.s)-1))
names(bl9)=c('from','to')
bl10=cbind.data.frame(names(d1.s)[-14],rep("Education",ncol(d1.s)-1))
names(bl10)=c('from','to')
bl=rbind(bl1,bl2,bl3,bl4,bl5,bl6,bl7,bl8,bl9,bl10)

for(i in 1:14){
	d1.s[,i]=as.numeric(d1.s[,i])
}
d1.s$Education=as.factor(ifelse(d1.s$Education <=4,0,1))

#run bnlearn
#use boostrapped method to get average network
set.seed(21)
boot=boot.strength(data=d1.s, R=500, algorithm='hc',algorithm.args=list(score='aic-cg',blacklist=bl))
avg.boot=averaged.network(boot)
sig=averaged.network(boot)$learning$args$threshold
avg.boot$arcs
avg.boot1=avg.boot
#equivalent class
cpdag(avg.boot)
cpdag(avg.boot)$arcs

#lm from bn
summary(lm(Cognition0~Age+Education,data=d1.s))
summary(lm(Cognition1~Cognition0+Sleep1,data=d1.s))
summary(lm(Fatigue0~Mood0,data=d1.s))
summary(lm(Fatigue1~Mood1,data=d1.s))
summary(lm(FOSQ0~Fatigue0,data=d1.s))
summary(lm(FOSQ1~Fatigue1+QoL0+FOSQ0+Education,data=d1.s))
summary(lm(Mood1~QoL0,data=d1.s))
summary(lm(QoL0~Fatigue0,data=d1.s))
summary(lm(QoL1~Mood1+Fatigue1+QoL0,data=d1.s))
summary(lm(Sleep0~QoL0,data=d1.s))
summary(lm(Sleep1~Fatigue1+Sleep0,data=d1.s))

#plot network
links=as.data.frame(avg.boot$arcs)
links$strength=boot[(boot$strength>sig)&(boot$direction>0.5),]$strength
links$direction=boot[(boot$strength>sig)&(boot$direction>0.5),]$direction
links$direction=ifelse(links$direction<=0.57,0.6,links$direction)
links=links[order(links$to),]
write.csv(links,file='links.csv')

library(igraph)
library('RColorBrewer')
nodes=data.frame(names(d1.s))
names(nodes)=c('variable')
nodes$type=c(1,2,1,2,1,2,1,2,1,2,1,2,3,3)
links=as.data.frame(avg.boot$arcs)
links$strength=boot[(boot$strength>sig)&(boot$direction>0.5),]$strength
links$direction=boot[(boot$strength>sig)&(boot$direction>0.5),]$direction
links$direction=ifelse(links$direction<=0.57,0.6,links$direction)
links$colind=ceiling((links$direction-0.5)*100/(50/8))
links$sign=c(1,1,1,-1,-1,    -1,-1,1,-1,-1,     1,-1,-1,1,-1,     1,-1,1,-1)

net <- graph.data.frame(links, nodes, directed=T) 
colrs=c("darkolivegreen3","tomato","gray60")
V(net)$color=colrs[V(net)$type]
V(net)$frame.color="white"
E(net)$width = E(net)$strength*5
E(net)$arrow.size=1
l <- layout.fruchterman.reingold(net)
l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
pa1=brewer.pal(8,"Blues")
pa2=brewer.pal(8,"Reds")
plot(net,rescale=F, layout=l*1.25, edge.color=ifelse(links$sign>0,pa1[links$colind],pa2[links$colind]),vertex.label.font=2,vertex.label.color="gray30",vertex.label.cex=0.7) 
library("SDMTools")
pnts = cbind(x =c(1.7,1.9,1.9,1.7)+0.2, y =c(1.4,1.4,1,1))
legend.gradient(pnts,col=pa1,c("Low","High"),title="Positive Direction",cex=0.7,text.col="gray30")
pnts2 = cbind(x =c(1.7,1.9,1.9,1.7)+0.2, y =c(0.9,0.9,0.5,0.5)-0.2)
legend.gradient(pnts2,col=pa2,c("Low","High"),title="Negative Direction",cex=0.7)
legend(-2.75,1.5, legend=c("Pre-chemotherapy","End-of-chemotherapy", "Demographics"),col=c("darkolivegreen3","tomato","gray60"),fill=colrs,ce=0.7,box.col="white")
legend(1.88,-0.8,legend=c("Low",rep("",5),"High"),lwd=1:7,title="Edge Strength",col="grey",box.col="white",cex=0.7)

#plot equivalent class network
nodes=data.frame(names(d1.s))
names(nodes)=c('variable')
nodes$type=c(1,2,1,2,1,2,1,2,1,2,1,2,3,3)
links=as.data.frame(cpdag(avg.boot)$arcs)
links$mode=c(2,0,0,0,0,   2,2,2,2,2,    0,2,0,2,2,   0,2,2,0,2,   2,2,2)
links$col=ifelse(links$mode>0,2,1)

net <- graph.data.frame(links, nodes, directed=T) 
colrs=c("darkolivegreen3","tomato","gray60")
V(net)$color=colrs[V(net)$type]
V(net)$frame.color="white"
l <- layout.fruchterman.reingold(net)
l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
E(net)$arrow.mode=E(net)$mode
colrs.e=c("gray70","black")
plot(net,edge.color=colrs.e[links$col],rescale=F, layout=l*1.25, vertex.label.font=2,vertex.label.color="gray30",vertex.label.cex=0.7)
legend(-2.75,1.5, legend=c("Pre-chemotherapy","End-of-chemotherapy", "Demographics"),col=c("darkolivegreen3","tomato","gray60"),fill=colrs,ce=0.7,box.col="white")
legend(1.8,1.4,legend=c("Directed Edge","Ambiguous Edge"),col=c("black","grey70"),box.col="white",cex=0.7,lty=1)



#########################################################
#build a BN model for phase 1 to 2
#########################################################
d1=data.long[!data.long$phase==0,]
Cognition1=d1[d1$phase==1,'compscore']
Cognition2=d1[d1$phase==2,'compscore']
Mood1=d1[d1$phase==1,'cesdtotal']
Mood2=d1[d1$phase==2,'cesdtotal']
Fatigue1=d1[d1$phase==1,'mfsitotal']
Fatigue2=d1[d1$phase==2,'mfsitotal']
QoL1=d1[d1$phase==1,'factbtotal']
QoL2=d1[d1$phase==2,'factbtotal']
Sleep1=d1[d1$phase==1,'psqitotal']
Sleep2=d1[d1$phase==2,'psqitotal']
FOSQ1=d1[d1$phase==1,'fosqtotal']
FOSQ2=d1[d1$phase==2,'fosqtotal']
Age=d1[d1$phase==1,'cb_demage']
Education=d1[d1$phase==1,'cb_demeducation']

d1.s=cbind.data.frame(Cognition1,Cognition2,Mood1,Mood2, Fatigue1,Fatigue2,QoL1,QoL2,Sleep1,Sleep2,FOSQ1,FOSQ2,Age,Education)

#create a blacklist disallowing arrows going certain direction based on
#common sense such as phase constraint
m=(ncol(d1.s)-2)/2
bl1=cbind.data.frame(names(d1.s)[c(2,4,6,8,10,12)],rep("Cognition1",m))
names(bl1)=c('from','to')
bl2=cbind.data.frame(names(d1.s)[c(2,4,6,8,10,12)],rep("Mood1",m))
names(bl2)=c('from','to')
bl3=cbind.data.frame(names(d1.s)[c(2,4,6,8,10,12)],rep("Fatigue1",m))
names(bl3)=c('from','to')
bl4=cbind.data.frame(names(d1.s)[c(2,4,6,8,10,12)],rep("QoL1",m))
names(bl4)=c('from','to')
bl5=cbind.data.frame(names(d1.s)[c(2,4,6,8,10,12)],rep("Sleep1",m))
names(bl5)=c('from','to')
bl6=cbind.data.frame(names(d1.s)[c(2,4,6,8,10,12)],rep("FOSQ1",m))
names(bl6)=c('from','to')
bl7=cbind.data.frame(rep("Cognition1",m*2-2),names(d1.s)[-c(1,2,13,14)])
names(bl7)=c('from','to')
bl8=cbind.data.frame(rep("Cognition2",m-1),names(d1.s)[c(4,6,8,10,12)])
names(bl8)=c('from','to')
bl9=cbind.data.frame(names(d1.s)[-13],rep("Age",ncol(d1.s)-1))
names(bl9)=c('from','to')
bl10=cbind.data.frame(names(d1.s)[-14],rep("Education",ncol(d1.s)-1))
names(bl10)=c('from','to')
bl=rbind(bl1,bl2,bl3,bl4,bl5,bl6,bl7,bl8,bl9,bl10)

for(i in 1:14){
	d1.s[,i]=as.numeric(d1.s[,i])
}
d1.s$Education=as.factor(ifelse(d1.s$Education <=4,0,1))

#run bnlearn
#use boostrapped method to get average network
set.seed(27)
boot=boot.strength(data=d1.s, R=500, algorithm='hc',algorithm.args=list(score='aic-cg',blacklist=bl))
avg.boot=averaged.network(boot)
sig=averaged.network(boot)$learning$args$threshold
avg.boot$arcs
mb(avg.boot,'compscore0')
mb(avg.boot,'compscore1')
#get equivalent class
cpdag(avg.boot)
cpdag(avg.boot)$arcs

#lm from bn
summary(lm(Cognition1~Age,data=d1.s))
summary(lm(Cognition2~Cognition1+Sleep2+Education,data=d1.s))
summary(lm(Fatigue1~Mood1,data=d1.s))
summary(lm(Fatigue2~Fatigue1+QoL2+Sleep1,data=d1.s))
summary(lm(FOSQ1~Fatigue1,data=d1.s))
summary(lm(FOSQ2~QoL2+FOSQ1,data=d1.s))
summary(lm(Mood2~Mood1+Fatigue1+Fatigue2+Education,data=d1.s))
summary(lm(QoL1~Mood1+Fatigue1,data=d1.s))
summary(lm(QoL2~QoL1,data=d1.s))
summary(lm(Sleep1~Fatigue1,data=d1.s))
summary(lm(Sleep2~Mood2+QoL2,data=d1.s))


#plot network
links=as.data.frame(avg.boot$arcs)
links$strength=boot[(boot$strength>sig)&(boot$direction>0.5),]$strength[-5]
links$direction=boot[(boot$strength>sig)&(boot$direction>0.5),]$direction[-5]
links$direction=ifelse(links$direction<=0.57,0.6,links$direction)
links2=links[order(links$to),]
write.csv(links2,file='links2.csv')

library(igraph)
library('RColorBrewer')
nodes=data.frame(names(d1.s))
names(nodes)=c('variable')
nodes$type=c(1,2,1,2,1,2,1,2,1,2,1,2,3,3)
links=as.data.frame(avg.boot$arcs)
links$strength=boot[(boot$strength>sig)&(boot$direction>0.5),]$strength[-5]
links$direction=boot[(boot$strength>sig)&(boot$direction>0.5),]$direction[-5]
links$direction=ifelse(links$direction<=0.57,0.6,links$direction)
links$colind=ceiling((links$direction-0.5)*100/(50/8))
links$sign=c(1,1,1,-1,1,    -1,1,-1,1,-1,     1,1,-1,-1,1,     1,-1,1,-1,-1,-1)

net <- graph.data.frame(links, nodes, directed=T) 
colrs=c("darkolivegreen3","tomato","gray60")
V(net)$color=colrs[V(net)$type]
V(net)$frame.color="white"
E(net)$width = E(net)$strength*5
E(net)$arrow.size=1
l <- layout.fruchterman.reingold(net)
l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
pa1=brewer.pal(8,"Blues")
pa2=brewer.pal(8,"Reds")
plot(net,rescale=F, layout=l*1.25, edge.color=ifelse(links$sign>0,pa1[links$colind],pa2[links$colind]),vertex.label.font=2,vertex.label.color="gray30",vertex.label.cex=0.7) 
library("SDMTools")
pnts = cbind(x =c(1.7,1.9,1.9,1.7)+0.2, y =c(1.4,1.4,1,1))
legend.gradient(pnts,col=pa1,c("Low","High"),title="Positive Direction",cex=0.7,text.col="gray30")
pnts2 = cbind(x =c(1.7,1.9,1.9,1.7)+0.2, y =c(0.9,0.9,0.5,0.5)-0.2)
legend.gradient(pnts2,col=pa2,c("Low","High"),title="Negative Direction",cex=0.7)
legend(-2.75,1.5, legend=c("Pre-chemotherapy","End-of-chemotherapy", "Demographics"),col=c("darkolivegreen3","tomato","gray60"),fill=colrs,ce=0.7,box.col="white")
legend(1.88,-0.8,legend=c("Low",rep("",5),"High"),lwd=1:7,title="Edge Strength",col="grey",box.col="white",cex=0.7)

#plot equivalent class network
nodes=data.frame(names(d1.s))
names(nodes)=c('variable')
nodes$type=c(1,2,1,2,1,2,1,2,1,2,1,2,3,3)
links=as.data.frame(cpdag(avg.boot)$arcs)
links$mode=c(2,0,2,0,0,   2,0,2,2,0,    0,0,2,0,0,   0,2,0,2,2,   0,2,2,0,2,   0,2,2)
links$col=ifelse(links$mode>0,2,1)

net <- graph.data.frame(links, nodes, directed=T) 
colrs=c("darkolivegreen3","tomato","gray60")
V(net)$color=colrs[V(net)$type]
V(net)$frame.color="white"
l <- layout.fruchterman.reingold(net)
l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
E(net)$arrow.mode=E(net)$mode
colrs.e=c("gray70","black")
plot(net,edge.color=colrs.e[links$col],rescale=F, layout=l*1.25, vertex.label.font=2,vertex.label.color="gray30",vertex.label.cex=0.7)
legend(-2.75,1.5, legend=c("Pre-chemotherapy","End-of-chemotherapy", "Demographics"),col=c("darkolivegreen3","tomato","gray60"),fill=colrs,ce=0.7,box.col="white")
legend(1.8,1.4,legend=c("Directed Edge","Ambiguous Edge"),col=c("black","grey70"),box.col="white",cex=0.7,lty=1)



#########################################################
#analyze performance sensitivity (in terms of aic, bic, loglik) 
#to changes in network
#########################################################

#original model phase0-1
score(avg.boot,data=d1.s,type='aic-cg')
score(avg.boot,data=d1.s,type='bic-cg')
#-1538.393
score(avg.boot,data=d1.s,type='loglik-cg')
#-1469.692

#original model phase1-2
score(avg.boot,data=d1.s,type='aic-cg')
score(avg.boot,data=d1.s,type='bic-cg')
#-1528.772
score(avg.boot,data=d1.s,type='loglik-cg')
#-1452.644


#isolate compscore1 and factbtotal1 in BL-A4
mat=as.data.frame(avg.boot$arcs)
mat1=mat[mat$from!="Cognition1",]
mat1=mat1[mat1$to!="Cognition1",]
dag = empty.graph(names(d1.s))
arcs(dag) = as.matrix(mat1)
score(dag,data=d1.s,type='aic-cg') #-1571.412
score(dag,data=d1.s,type='bic-cg')  #-1571.412 decrease by 33
score(dag,data=d1.s,type='loglik-cg') #-1506.424

mat=as.data.frame(avg.boot$arcs)
mat1=mat[mat$from!="QoL1",]
mat1=mat1[mat1$to!="QoL1",]
dag = empty.graph(names(d1.s))
arcs(dag) = as.matrix(mat1)
score(dag,data=d1.s,type='aic-cg') #-1576.86
score(dag,data=d1.s,type='bic-cg')  #-1576.86 decrease by 38
score(dag,data=d1.s,type='loglik-cg') #-1513.729


#isolate compscore2 and factbtotal2 in A4-Y1
mat=as.data.frame(avg.boot$arcs)
mat1=mat[mat$from!="Cognition2",]
mat1=mat1[mat1$to!="Cognition2",]
dag = empty.graph(names(d1.s))
arcs(dag) = as.matrix(mat1)
score(dag,data=d1.s,type='aic-cg') #-1580.946
score(dag,data=d1.s,type='bic-cg')  #-1580.946 decrease by 52
score(dag,data=d1.s,type='loglik-cg') #-1515.958

mat=as.data.frame(avg.boot$arcs)
mat1=mat[mat$from!="QoL2",]
mat1=mat1[mat1$to!="QoL2",]
dag = empty.graph(names(d1.s))
arcs(dag) = as.matrix(mat1)
score(dag,data=d1.s,type='aic-cg') #-1567.101
score(dag,data=d1.s,type='bic-cg')  #-1567.101 decrease by 38
score(dag,data=d1.s,type='loglik-cg') #-1500.257


#remove sleep psqi to compscore
#0-1
dag=drop.arc(avg.boot,from="Sleep1",to="Cognition1")
score(dag,data=d1.s,type='aic-cg') #-1540.823
score(dag,data=d1.s,type='bic-cg')  #-1540.823 decrease by 2.43
score(dag,data=d1.s,type='loglik-cg') #-1473.979

#1-2
dag=drop.arc(avg.boot,from="Sleep2",to="Cognition2")
score(dag,data=d1.s,type='aic-cg') #-1543.586
score(dag,data=d1.s,type='bic-cg')  #-1543.586 decrease by 15
score(dag,data=d1.s,type='loglik-cg') #-1473.028


#remove mood1-qol1, fatigue1-qol1
dag=drop.arc(avg.boot,from="Mood1",to="QoL1")
score(dag,data=d1.s,type='aic-cg') #-1544.328
score(dag,data=d1.s,type='bic-cg')  #-1544.328 decrease by 6
score(dag,data=d1.s,type='loglik-cg') #-1477.483

dag=drop.arc(avg.boot,from="Fatigue1",to="QoL1")
score(dag,data=d1.s,type='aic-cg') #-1539.517
score(dag,data=d1.s,type='bic-cg')  #-1539.517 decrease by 1.124
score(dag,data=d1.s,type='loglik-cg') #-1472.673


#########################################################
#build a BN model for all phases
#########################################################
d1=data.long
Cognition0=d1[d1$phase==0,'compscore']
Cognition1=d1[d1$phase==1,'compscore']
Cognition2=d1[d1$phase==2,'compscore']
Mood0=d1[d1$phase==0,'cesdtotal']
Mood1=d1[d1$phase==1,'cesdtotal']
Mood2=d1[d1$phase==2,'cesdtotal']
Fatigue0=d1[d1$phase==0,'mfsitotal']
Fatigue1=d1[d1$phase==1,'mfsitotal']
Fatigue2=d1[d1$phase==2,'mfsitotal']
QoL0=d1[d1$phase==0,'factbtotal']
QoL1=d1[d1$phase==1,'factbtotal']
QoL2=d1[d1$phase==2,'factbtotal']
Sleep0=d1[d1$phase==0,'psqitotal']
Sleep1=d1[d1$phase==1,'psqitotal']
Sleep2=d1[d1$phase==2,'psqitotal']
FOSQ0=d1[d1$phase==0,'fosqtotal']
FOSQ1=d1[d1$phase==1,'fosqtotal']
FOSQ2=d1[d1$phase==2,'fosqtotal']
Age=d1[d1$phase==1,'cb_demage']
Education=d1[d1$phase==1,'cb_demeducation']

d1.s=cbind.data.frame(Cognition0,Cognition1,Cognition2,Mood0,Mood1,Mood2,Fatigue0,Fatigue1,Fatigue2,QoL0,QoL1,QoL2,Sleep0,Sleep1,Sleep2,FOSQ0,FOSQ1,FOSQ2,Age,Education)

#create a blacklist disallowing arrows going certain direction based on
#common sense such as phase constraint
m=(ncol(d1.s)-2)*2/3
bl1=cbind.data.frame(names(d1.s)[c(2,3,5,6,8,9,11,12,14,15,17,18)],rep("Cognition0",m))
names(bl1)=c('from','to')
bl2=cbind.data.frame(names(d1.s)[c(2,3,5,6,8,9,11,12,14,15,17,18)],rep("Mood0",m))
names(bl2)=c('from','to')
bl3=cbind.data.frame(names(d1.s)[c(2,3,5,6,8,9,11,12,14,15,17,18)],rep("Fatigue0",m))
names(bl3)=c('from','to')
bl4=cbind.data.frame(names(d1.s)[c(2,3,5,6,8,9,11,12,14,15,17,18)],rep("QoL0",m))
names(bl4)=c('from','to')
bl5=cbind.data.frame(names(d1.s)[c(2,3,5,6,8,9,11,12,14,15,17,18)],rep("Sleep0",m))
names(bl5)=c('from','to')
bl6=cbind.data.frame(names(d1.s)[c(2,3,5,6,8,9,11,12,14,15,17,18)],rep("FOSQ0",m))
names(bl6)=c('from','to')
bl7=cbind.data.frame(names(d1.s)[c(3,6,9,12,15,18)],rep("Cognition1",m/2))
names(bl7)=c('from','to')
bl8=cbind.data.frame(names(d1.s)[c(3,6,9,12,15,18)],rep("Mood1",m/2))
names(bl8)=c('from','to')
bl9=cbind.data.frame(names(d1.s)[c(3,6,9,12,15,18)],rep("Fatigue1",m/2))
names(bl9)=c('from','to')
bl10=cbind.data.frame(names(d1.s)[c(3,6,9,12,15,18)],rep("QoL1",m/2))
names(bl10)=c('from','to')
bl11=cbind.data.frame(names(d1.s)[c(3,6,9,12,15,18)],rep("Sleep1",m/2))
names(bl11)=c('from','to')
bl12=cbind.data.frame(names(d1.s)[c(3,6,9,12,15,18)],rep("FOSQ1",m/2))
names(bl12)=c('from','to')
bl13=cbind.data.frame(rep("Cognition0",ncol(d1.s)-5),names(d1.s)[-c(1,2,3,19,20)])
names(bl13)=c('from','to')
bl14=cbind.data.frame(rep("Cognition1",m-2),names(d1.s)[c(5,6,8,9,11,12,14,15,17,18)])
names(bl14)=c('from','to')
bl15=cbind.data.frame(rep("Cognition2",5),names(d1.s)[c(6,9,12,15,18)])
names(bl15)=c('from','to')
bl16=cbind.data.frame(names(d1.s)[-19],rep("Age",ncol(d1.s)-1))
names(bl16)=c('from','to')
bl17=cbind.data.frame(names(d1.s)[-20],rep("Education",ncol(d1.s)-1))
names(bl17)=c('from','to')
#bl18=cbind.data.frame(rep("Education",ncol(d1.s)-1),names(d1.s)[-20])
#names(bl18)=c('from','to')

bl=rbind(bl1,bl2,bl3,bl4,bl5,bl6,bl7,bl8,bl9,bl10,bl11,bl12,bl13,bl14,bl15,bl16,bl17)#,bl18)

for(i in 1:20){
	d1.s[,i]=as.numeric(d1.s[,i])
}
d1.s$Education=as.factor(ifelse(d1.s$Education <=4,0

#run bnlearn
#use boostrapped method to get average network
set.seed(1)
boot=boot.strength(data=d1.s, R=500, algorithm='hc',algorithm.args=list(score='bic-cg',blacklist=bl))
avg.boot=averaged.network(boot)
sig=averaged.network(boot)$learning$args$threshold
avg.boot$arcs
result=boot[(boot$strength>sig)&(boot$direction>0.5),]
result[order(result$to),]


#check performance
score(avg.boot,data=d1.s,type='bic-cg')
score(avg.boot,data=d1.s,type='aic-cg')
score(avg.boot,data=d1.s,type='loglik-cg')

#apply phase0-1 model to phase1-2 model
mat=as.data.frame(avg.boot$arcs)
mat$from=as.character(mat$from)
mat1=mat[-c(2,8,11,19,20,22,24,26),]
mat1$from=as.character(mat1$from)
mat1$to=as.character(mat1$to)
mat1[21,]=c('Fatigue2','Sleep2')
mat1[22,]=c('QoL1','QoL2')
mat1[23,]=c('QoL1','Mood2')
mat1[24,]=c('QoL1','FOSQ2')
mat1[25,]=c('QoL1','Sleep1')
mat1[26,]=c('Sleep1','Sleep2')
mat1[27,]=c('Sleep2','Cognition2')
dag = empty.graph(names(d1.s))
arcs(dag) = as.matrix(mat1)
score(dag,data=d1.s,type='aic-cg') #-2230.79
score(dag,data=d1.s,type='bic-cg')*(-2)  #4461.581
score(dag,data=d1.s,type='loglik-cg') #-2143.521

#isolate cognition
mat1=mat[mat$from!="Cognition0",]
mat1=mat1[mat1$to!="Cognition0",]
mat1=mat1[mat1$from!="Cognition1",]
mat1=mat1[mat1$to!="Cognition1",]
mat1=mat1[mat1$from!="Cognition2",]
mat1=mat1[mat1$to!="Cognition2",]
dag = empty.graph(names(d1.s))
arcs(dag) = as.matrix(mat1)
score(dag,data=d1.s,type='aic-cg') #-2230.79
score(dag,data=d1.s,type='bic-cg')*(-2)  #4461.581
score(dag,data=d1.s,type='loglik-cg') #-2143.521

#isolate QoL
mat1=mat[mat$from!="QoL0",]
mat1=mat1[mat1$to!="QoL0",]
mat1=mat1[mat1$from!="QoL1",]
mat1=mat1[mat1$to!="QoL1",]
mat1=mat1[mat1$from!="QoL2",]
mat1=mat1[mat1$to!="QoL2",]
dag = empty.graph(names(d1.s))
arcs(dag) = as.matrix(mat1)
score(dag,data=d1.s,type='aic-cg') #-2230.79
score(dag,data=d1.s,type='bic-cg')*(-2)  #4461.581
score(dag,data=d1.s,type='loglik-cg') #-2143.521





















