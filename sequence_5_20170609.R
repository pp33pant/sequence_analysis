
library(TraMineR)
library(cluster)
library(bigmemory.sri)
library(bigmemory)
library(nnet)
library(stargazer)

# put in the tables
# the sequences for the family statuses
all=read.table("H:\\all samples.csv",sep=",",header=TRUE)


fasall.seq <- seqdef(all, 3:32,labels= c("not married, no child","not married, children","married, one child","married, two kids","married, three kids and more"))

dist.fasall <- seqdist(fasall.seq,method="OM",indel = 1, with.missing = TRUE, sm="TRATE")
windows()
seqIplot(fasall.seq, border = NA, group = all$gcse5eq, sortv=dist.fasall, cex.legend=0.45)

# output the transition rates
fasall.trate<-seqtrate(fasall.seq)
round(fasall.trate,2)

#the sequences for the location statuses 
#locall.seq=seqdef(all, 34:63, labels=c("NA", "north", "northeast", "east", "center", "southwest", "northwest","foreign"))
#dist.locall<- seqdist(locall.seq, method="OM", indel = 1, with.missing=TRUE, sm="TRATE")
#seqIplot(locall.seq, border = NA, group = all$gcse5eq, sortv=dist.locall, cex.legend=0.45)

#the sequences for the working status 
wksall.seq<- seqdef(all, 64:93, labels=c("no work", "schooling", "farming","private","state and collective"))
dist.wksall<-seqdist(wksall.seq, method="OM", indel = 1, with.missing=TRUE, sm="TRATE")
windows()
seqIplot(wksall.seq,border = NA, group = all$gcse5eq, sortv=dist.wksall, cex.legend=0.45)

# the general method to scale the multichannels (to sum the total distances)
dist.all<- dist.fasall+dist.wksall#+dist.locall

#all.mds <- cmdscale(dist.all,k=5, eig=TRUE)  #note: command "cmdscale" only matches the OM method with indel=1

#cluster analysis 
all.clusterward<- agnes(dist.all, diss=T, method="ward")
windows()
plot(all.clusterward,ask=F, which.plots=2)
all.cl5<-cutree(all.clusterward, k=5) # the cluster tree
all.pam5<-pam(dist.all, k=5, diss=T)
windows()
plot(all.pam5)
clust.labels<-c("rural+fewer kids","urban+fewer kids","urban+more kids","rural+more kids","missing+fewer kids")
all.cl5.factor<-factor(all.cl5, levels=1:5, labels=clust.labels)

windows()
seqIplot(fasall.seq, group=all.cl5.factor, cex.legend=0.45)
windows()
seqdplot(fasall.seq, group=all.cl5.factor, cex.legend=0.45)
windows()
seqIplot(wksall.seq, group=all.cl5.factor, cex.legend=0.45)
windows()
seqdplot(wksall.seq, group=all.cl5.factor, cex.legend=0.45)

# to define the types using the multi-logistic regression method
status<-all.cl5.factor
dcohort<-class.ind(all$cohort)
dloc<-class.ind(all$lage16)
male<- all$rgender==1
mlogit<-multinom(status~male +cohort1 +cohort2 +cohort3 +cohort4 +red +black +fedu +medu +primary +middle +high +minor +party+ loc1+ loc2+ loc3+ loc4+ loc5+ loc6, data=all)
stargazer(mlogit,type = "text",title = "mlogit")



# scatterplot (MDS) multidimensional scaling 
all2d<- cmdscale(dist.all, k=2)
windows()
plot(all2d, pch=all.cl5, col=all.cl5)
legend("bottomright", col=1:5, pch=1:5, legend=clust.labels)

#first factor of MDS analysis 
windows()
seqIplot(fasall.seq, group=all.cl5.factor, sortv=all2d[,1], cex.legend=0.45)
windows()
seqIplot(wksall.seq, group=all.cl5.factor, sortv=all2d[,1], cex.legend=0.45)

#the index of sequence analysis
dissvar(dist.all)
#ANOVA according to the genders of the interviewers
gender<-dissassoc(dist.all, group = all$rgender, R=1000)
windows()
hist(gender,col="cyan")

fas.diff<-seqdiff(fasall.seq ,group = all.cl5.factor, with.missing=T)
wks.diff<-seqdiff(wksall.seq, group=all.cl5.factor, with.missing = T)
windows()
plot(fas.diff,lwd=3, col="darkred")
windows()
plot(wks.diff, lwd=3, col="yellow")
windows()
plot(fas.diff,lwd=3,stat="Variance",legend.pos = "bottom")
windows()
plot(wks.diff,lwd=3,stat="Variance", legend.pos = "right")
# the discrepancy tree
dt<- disstree(dist.all ~ male +cohort1 +cohort2 +cohort3 +cohort4 +red +black +fedu +medu +primary +middle +high +minor +party+ loc1+ loc2+ loc3+ loc4+ loc5+ loc6, data=all, R=100) #here the dummy variables are forbidden
print(dt)
# to print the discrepancy tree
windows()
seqtree2dot(dt,"fg_allseqtree", seqdata = wksall.seq, type="d", border=NA, withlegend = F, axes = F, ylab="",yaxis=F)
system ("dot -Tsvg-ofg_allseqtree.svg fg_allseqtree.dot")
system ("convert fg_allseqtree.svg fg_allseqtree.jpg")

## GIMSA
wks.mds<-cmdscale(dist.wksall, k=4, eig=T)
fas.mds<-cmdscale(dist.fasall, k=5,eig=T)
stress <- function(omres,mdsres){ 
  datadist <- as.dist(omres) 
  res <- numeric(length=ncol(mdsres$points)) 
  for(i in 1:ncol(mdsres$points)) { 
    fitteddist <- dist(mdsres$points[,1:i],diag=TRUE,upper=TRUE) 
    res[i] <- sqrt(sum((datadist-fitteddist)^2)/sum(datadist^2))
  }
  res
}
stress(dist.fasall,fas.mds)
fas.mds$eig[1:5]/fas.mds$eig[1]
stress(dist.wksall,wks.mds)
wks.mds$eig[1:4]/wks.mds$eig[1]
windows()
seqIplot(fasall.seq, sort=fas.mds$points[,1], cex.legend=0.45)
windows()
seqIplot(wksall.seq, sort=wks.mds$points[,1], cex.legend=0.45)
nbmds.wks<-4
nbmds.fas<-4  ##according to the results of factor analysis above (eig>0.10)

##PLS method 
a<-wks.mds$points[,1:nbmds.fas]
b<-fas.mds$points[,1:nbmds.wks]
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

#matrix clustering 
# no weighting(w0) 
F <- pls$F 
G <- pls$G 

# weighting by variance of PLS components(w1) 
F <- apply(pls$F,2,scale,center=FALSE) 
G <- apply(pls$G,2,scale,center=FALSE) 

# weighting by number of distinct sequences(w2) 
F <- pls$F/nrow(seqtab(wksall.seq,tlim=0)) 
G <- pls$G/nrow(seqtab(fasall.seq,tlim=0)) 

# weighting by MDS 1st eigenvalue(w3) 
F <- pls$F/wks.mds$eig[1] 
G <- pls$G/fas.mds$eig[1] 

# distance computation 
diff2 <- function(X) return(as.matrix(dist(X,upper=T,diag=T)^2,nrow=nrow(X))) 
D <- (diff2(F)+diff2(G))^0.5 

# clustering 
seq.dist <- as.dist(D) 
seq.agnes <- agnes(seq.dist,method="ward",keep.diss=FALSE) 
windows()
plot(seq.agnes, ask=F, which.plots = 2)

seq.part <- cutree(seq.agnes,nbcl <- 10) 
seq.part.factor<- factor(seq.part, levels=1:5, labels=clust.labels)
windows()
seqIplot(fasall.seq, group=seq.part.factor, cex.legend=0.45)
windows()
seqdplot(fasall.seq, group=seq.part.factor, cex.legend=0.45)
windows()
seqIplot(wksall.seq, group=seq.part.factor, cex.legend=0.45)
windows()
seqdplot(wksall.seq, group=seq.part.factor, cex.legend=0.45)

# to define the types using the multi-logistic regression method
status<-seq.part.factor
dcohort<-class.ind(all$cohort)
dloc<-class.ind(all$lage16)
male<- all$rgender==1
mlogit<-multinom(status~male +cohort1 +cohort2 +cohort3 +cohort4 +red +black +fedu +medu +primary +middle +high +minor +party+ loc1+ loc2+ loc3+ loc4+ loc5+ loc6, data=all)
stargazer(mlogit,type = "text",title = "mlogit")


# scatterplot (MDS) multidimensional scaling 
all2d<- cmdscale(dist.all, k=2)
windows()
plot(all2d, pch=seq.part, col=seq.part)
legend("bottomright", col=1:5, pch=1:5, legend=clust.labels)

#first factor of MDS analysis 
windows()
seqIplot(fasall.seq, group=seq.part.factor, sortv=all2d[,1], cex.legend=0.45)
windows()
seqIplot(wksall.seq, group=seq.part.factor, sortv=all2d[,1], cex.legend=0.45)

#the index of sequence analysis
dissvar(dist.all)
#ANOVA according to the genders of the interviewers
gender<-dissassoc(dist.all, group = all$rgender, R=1000)
windows()
hist(gender,col="cyan")

fas.diff<-seqdiff(fasall.seq ,group = seq.part.factor, with.missing=T)
wks.diff<-seqdiff(wksall.seq, group=seq.part.factor, with.missing = T)
windows()
plot(fas.diff,lwd=3, col="darkred")
windows()
plot(wks.diff, lwd=3, col="yellow")
windows()
plot(fas.diff,lwd=3,stat="Variance")
windows()
plot(wks.diff,lwd=3,stat="Variance")

# the discrepancy tree
dt<- disstree(dist.all ~male +cohort1 +cohort2 +cohort3 +cohort4 +red +black +fedu +medu +primary +middle +high +minor +party+ loc1+ loc2+ loc3+ loc4+ loc5+ loc6, data=all, R=100) #here the dummy variables are forbidden
print(dt)
# to print the discrepancy tree
windows()
seqtree2dot(dt,"fg_allseqtree", seqdata = wksall.seq, type="d", border=NA, withlegend = F, axes = F, ylab="",yaxis=F)
system ("dot -Tsvg-ofg_allseqtree.svg fg_allseqtree.dot")
system ("convert fg_allseqtree.svg fg_allseqtree.jpg")