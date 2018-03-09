rm(list=ls())
list.files()

imp=FALSE
imp=TRUE
##################################
##################################
##  	INSTALL PACKAGES                     ##
##################################
##################################

library(ebdbNet)
library(longitudinal)
library(mice)

###########################
###########################
##  	FUNCTIONS                     ##
###########################
###########################

st=function(x){
	z=(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
	return(z)
}

#######################################
#######################################
##     READ IN & MANIPULATE DATA             ##
#######################################
#######################################

setwd('/Users/Selene/Desktop/ebdbNet')
data=read.csv("chemodatabase_06252012.csv",header=TRUE)
data=data[data$group==1,]
data=data[!data$phase==2,]

names(data)
ind_keep=c(1:3,231,97,135,154,147,160);names(data[,ind_keep])
#ind_keep=c(1:3,97:103);names(data[,ind_keep])
P=length(ind_keep)-3
data=data[,ind_keep]
dim(data);summary(data)
id=unique(data$sub)
N=length(id)

## option: residualize on age, college

data.long=as.data.frame(array(NA,dim=c(N*2,dim(data)[2])))
names(data.long)=names(data)
data.long[,1]=rep(id,each=2)
data.long[,2]=rep(0:1,N)
for(i in 1:N){
	data.long.i=data.long[data.long$sub==id[i],]
	data.i=data[data$sub==id[i],]
	data.long.i[,3]=data.i[1,3]
	for(p in 1:P){
		for(t in 0:1){
			if(length(data.i[data.i$phase==t,p+3]==1)){
				data.long.i[t+1,p+3]=data.i[data.i$phase==t,p+3]
			}	
		}	
	}
	data.long[data.long$sub==id[i],]=data.long.i
}	

#either impute or delete incomplete cases
if(imp==FALSE){
	data.long.keep=rep(1,dim(data.long)[1])
	ind=1:2
	for(i in 1:N){
		if(sum(is.na(data.long[data.long$sub==id[i],(1+3):(P+3)]))>0){
			data.long.keep[ind]=0
		}
		ind=ind+2
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

dim(data.long)
summary(data.long)
id=unique(data.long$sub)
N=length(id)


#get data ready for the ebdbn package
load("/Users/Selene/Desktop/ebdbNet/data.long.rdata")
data.long.1=data.long[!data.long$phase==2,]
data.long.2=data.long[!data.long$phase==0,]

for(p in (1+3):(P+3)){
	data.long.1[,p]=st(data.long.1[,p])
}
for(p in (1+3):(P+3)){
	data.long.2[,p]=st(data.long.2[,p])
}

X.list.1=list()
for(i in 1:N){
	X.i=as.matrix(data.long.1[data.long.1$sub==id[i],])
	X.list.1[[i]]=t(X.i[,(1+3):(P+3)])
}	
X.list.2=list()
for(i in 1:N){
	X.i=as.matrix(data.long.2[data.long.2$sub==id[i],])
	X.list.2[[i]]=t(X.i[,(1+3):(P+3)])
}	


###############################
###############################
##  	RUN ALGORITHM                    ##
###############################
###############################

#estimate number of hidden states
K <- hankel(X.list.1, lag=1, cutoff=0.75)

#run ebdb algorithm
net.1 <- ebdbn(y=X.list.1, K=K$dim)

#visualize results
#pdf("new.edbdnet.results.pdf")
plot(net.1, sig.level=1-.05, clarify=T)#, layout=layout.lgl)


#estimate number of hidden states
K <- hankel(X.list.2, lag=1, cutoff=0.75)

#run ebdb algorithm
net.2 <- ebdbn(y=X.list.2, K=K$dim)

#visualize results
#pdf("new.edbdnet.results.pdf")
plot(net.2, sig.level=1-.05, clarify=T)#, layout=layout.lgl)


dev.off()





