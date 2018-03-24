library(TraMineR)
library(cluster)

seqparents=read.table("F:\\projects\\sequence analysis\\exp_3.csv",sep=",",header=TRUE)
seqkids=read.table(file="F:\\projects\\sequence analysis\\student1.csv",sep=",",fill=TRUE, header= TRUE)
head(seqparents)
head(seqkids)
data2=read.csv("F:\\projects\\sequence analysis\\exp_4.csv",sep=",",header=TRUE)

library(rgl)
library(car)

#scatter3d(pwage~ age+ cohort| pwage ,data=data2, surfacecol = c("blue", "green", "orange", "red", "yellow"), 
 #         ellipsoid=TRUE, parallel=FALSE) 
#qmesh3d(pwage, age)
# to generate the descriptive data (graph)

#sequence defination 
parent.seq<-seqdef(seqparents, 4:14, 
                   labels= c("no work","student","agri","private","public"))
student.seq<-seqdef(seqkids,2:6, 
                    labels=c("live with parents", "live in same commity", "live in a county", "live out"))
names(parent.seq)
# 1st GIMSA step: dissimilarity measures 
parent.cost <- seqsubm(parent.seq, method="CONSTANT", cval=NULL, with.missing = TRUE)
parent.om<-seqdist(head(parent.seq,5000), method="HAM", sm= parent.cost, with.missing = TRUE)
student.cost <- seqsubm(student.seq,method="CONSTANT", cval=NULL, with.missing = TRUE)
student.om <- seqdist(head(student.seq,5000), method = "HAM", 
                      sm=student.cost, with.missing = TRUE)

#2nd GIMSA step: multidimensional scaling
parent.mds <- cmdscale(parent.om,k=5,eig=TRUE)
student.mds <- cmdscale(student.om, k=5, eig=TRUE)
# To choose the number of MDS components to retain 
stress <- function(omres,mdsres){ 
  datadist <- as.dist(omres) 
  res <- numeric(length=ncol(mdsres$points)) 
  for(i in 1:ncol(mdsres$points)) { 
    fitteddist <- dist(mdsres$points[,1:i],diag=TRUE,upper=TRUE) 
    res[i] <- sqrt(sum((datadist-fitteddist)^2)/sum(datadist^2))
  }
  res
}
stress(parent.om,parent.mds) 
parent.mds$eig[1:10]/parent.mds$eig[1] 
seqIplot(parent.seq,sort=parent.mds$points[,1]) 
stress(student.om,student.mds) 
student.mds$eig[1:10]/student.mds$eig[1] 
seqIplot(student.seq,sort=student.mds$points[,1]) 

nbmds.parent <- 5 
nbmds.student <- 4 

## ================================ 
## 3rd GIMSA step:symetrical(ie canonical) PLS 
## ================================ 

a <- parent.mds$points[,1:nbmds.parent] 
b <- student.mds$points[,1:nbmds.student] 

symPLS <- function(a,b){ 
  k <- min(ncol(a),ncol(b),nrow(a),nrow(b)) 
  X <-vector("list", k) 
  
  Y <- vector("list",k) 
  X[[1]] <- scale(a,scale=FALSE) 
  Y[[1]] <- scale(b,scale=FALSE) 
  F <- matrix(nrow=nrow(X[[1]]),ncol=k) 
  G <- matrix(nrow=nrow(X[[1]]),ncol=k) 
  f <- matrix(nrow=nrow(X[[1]]),ncol=k) 
  g <- matrix(nrow=nrow(X[[1]]),ncol=k) 
  vF <- vector(mode="numeric", length=k) 
  vG <- vector(mode="numeric", length=k) 
  corr <- vector(mode="numeric",length=k) 
  for(i in 1:k){ 
    
    u <- eigen(t(X[[i]])%*%Y[[i]]%*%t(Y[[i]])%*%X[[i]])$vectors[,1] 
    F[,i] <- X[[i]]%*%u 
    v <- t(Y[[i]])%*%X[[i]]%*%u 
    v <- v*as.vector(1/((t(v)%*%v)^0.5)) 
    G[,i] <- Y[[i]]%*%v 
    f[,i] <- F[,i]*as.vector(1/((t(F[,i])%*%F[,i])^0.5)) 
    g[,i] <- G[,i]*as.vector(1/((t(G[,i])%*%G[,i])^0.5)) 
    X[[i+1]] <- X[[i]] - f[,i]%*%t(f[,i])%*%X[[i]] 
    Y[[i+1]] <- Y[[i]] - g[,i]%*%t(g[,i])%*%Y[[i]] 
    vF[i] <- var(F[,i]) 
    vG[i] <- var(G[,i]) 
    corr[i] <- cor(x=F[,i],y=G[,i],method="pearson") 
  } 
  
  res <- list(F=F,G=G,vF=vF,vG=vG,corr=corr) 
  rm(k,X,Y,f,g,u,v) 
  return(res) 
} 

pls <- symPLS(a,b) 

## ================================ 
## 4th GIMSA step: distance matrix and clustering 
## ================================ 

# no weighting(w0) 
F <- pls$F 
G <- pls$G 

# weighting by variance of PLS components(w1) 
F <- apply(pls$F,2,scale,center=FALSE) 
G <- apply(pls$G,2,scale,center=FALSE) 

# weighting by number of distinct sequences(w2) 
F <- pls$F/nrow(seqtab(parent.seq,tlim=0)) 
G <- pls$G/nrow(seqtab(student.seq,tlim=0)) 

# weighting by MDS 1st eigenvalue(w3) 
F <- pls$F/parent.mds$eig[1] 
G <- pls$G/parent.mds$eig[1] 

# distance computation 
diff2 <- function(X) return(as.matrix(dist(X,upper=T,diag=T)^2,nrow=nrow(X))) 
D <- (diff2(F)+diff2(G))^0.5 

# clustering 
seq.dist <- as.dist(D) 
seq.agnes <- agnes(seq.dist,method="ward",keep.diss=FALSE) 
seq.part <- cutree(seq.agnes,nbcl <- 10) 


