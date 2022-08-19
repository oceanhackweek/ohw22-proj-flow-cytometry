# This document contains functions that we find useful for gene co-expression network analysis
# We recommend to read the tutorial first.
# Steve Horvath, Bin Zhang, Jun Dong, Andy Yip
# To cite this code or the statistical methods please use
# Zhang B, Horvath S (2005) A General Framework for Weighted Gene Co-Expression Network Analysis.
# Statistical Applications in Genetics and Molecular Biology. In Press.
# Technical Report and software code at: www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork.

# Modifications by Peter langfelder: in function ScaleFreePlot1: changed min(kk) and max(kk) to min(kk, rm.na=T) and
# max(kk, na.rm=T).
# in function ModulePrinComps1 made the printed statement a bit more informative and use print.flush.
# Changed titling in ScaleFreePlot1

# Added lots of new functions for "higher level" analysis of several datasets at once

# CONTENTS
# This document contains function for carrying out the following tasks
# A) Assessing scale free topology and choosing the parameters of the adjacency function
#    using the scale free topology criterion (Zhang and Horvath 05)
# B) Computing the topological overlap matrix
# C) Defining gene modules using clustering procedures
# D) Summing up modules by their first principal component (first eigengene)
# E) Relating a measure of gene significance to the modules
# F) Carrying out a within module analysis (computing intramodular connectivity)
#    and relating intramodular connectivity to gene significance.
# G) Miscellaneous other functions, e.g. for computing the cluster coefficient.
# H) Network functions from Bin Zhang for dynamic tree cutting of a hierarchical clustering tree (dendrogram)
# I) General statistical functions, e.g. for scatterplots


#--------------------------------------------------------------------------------------------
# Set global parameters.

# This parameter controls whether my implementation of pmax should use a call to an external function
# (recommended if the requisite library is available). Otherwise an R-only implementation will be used,
# which is significantly slower (but available and stable on all R platforms).

UseCpmax = FALSE;

#--------------------------------------------------------------------------------------------

# Load the requisite libraries

WorkingDirectory = getwd();

#if (!library(MASS, logical.return=TRUE)) { # standard, no need to install
#For some reason, MASS does not seem to be installed on Titan, so we'll try to load it
# from my own library. If that fails as well, stop.
#  if (!library(MASS, logical.return=TRUE, lib.loc="M:/Work/RLibrary")) stop()
#}

library(MASS);
library(class)	# standard, no need to install
library(cluster)
#library(sma)	# install it for the function plot.mat
library(impute)# install it for imputing missing value
library(Hmisc)	# install it for the C-index calculations
library(survival)
library(fields);

#oldwd = getwd();

if (UseCpmax) source("../ComputerDefinition/ComputerDefinition.R");

if (exists("memory.limit"))
{
  # increase the available memory
  memory.limit(size=4000)
}

#setwd(oldwd);


#####################################################################################################
# A) Assessing scale free topology and choosing the parameters of the adjacency function
#    using the scale free topology criterion (Zhang and Horvath 05)
########################################################################################################


# ===================================================
#For hard thresholding, we use the signum (step) function
if(exists("signum1") ) rm(signum1);
signum1=function(corhelp,tau1)  {
adjmat1= as.matrix(abs(corhelp)>=tau1)
dimnames(adjmat1) <- dimnames(corhelp)
diag(adjmat1) <- 0
adjmat1}

# ===================================================
# For soft thresholding, one can use the sigmoid function
# But we have focused on the power adjacency function in the tutorial...
if (exists("sigmoid1") ) rm(sigmoid1); sigmoid1=function(ss,mu1=0.8,alpha1=20) {
1/(1+exp(-alpha1*(ss-mu1)))}





#This function is useful for speeding up the connectivity calculation.
#The idea is to partition the adjacency matrix into consecutive baches of a #given size.
#In principle, the larger the batch size the faster is the calculation. But #smaller batchsizes require #less memory...
# Input: gene expression data set where *rows* correspond to microarray samples #and columns correspond to genes.
# If fewer than MinimumNoSamples contain gene expression information for a given
# gene, then its connectivity is set to missing.
if(exists("SoftConnectivity")) rm(SoftConnectivity);
SoftConnectivity=function(datE, power=6, batchsize=1500, MinimumNoSamples=10) {
no.genes=dim(datE)[[2]]
no.samples=dim(datE)[[1]]
if (no.genes<no.samples | no.genes<10 | no.samples<5 ) {stop("Error: Something seems to be wrong. Make sure that the input data frame has genes as rows and array samples as columns. Alternatively, there could be fewer than 10 genes or fewer than 5 samples. ") } else {
sum1=function(x) sum(x,na.rm=T)
k=rep(NA,no.genes)
no.batches=as.integer(no.genes/ batchsize)
if (no.batches>0) {
for (i in 1:no.batches) {
print(paste("batch number = ", i))
index1=c(1:batchsize)+(i-1)* batchsize
ad1=abs(cor(datE[,index1], datE,use="p"))^power
ad1[is.na(ad1)]=0
k[index1]=apply(ad1,1,sum1)
# If fewer than MinimumNoSamples contain gene expression information for a given
# gene, then we set its connectivity to 0.
NoSamplesAvailable=apply(!is.na(datE[,index1]),2,sum)
k[index1][NoSamplesAvailable< MinimumNoSamples]=NA
} # end of for (i in 1:no.batches
} # end of if (no.batches>0)...
if (no.genes-no.batches*batchsize>0 ) {
restindex=c((no.batches*batchsize+1):no.genes)
ad1=abs(cor(datE[,restindex], datE,use="p"))^power
ad1[is.na(ad1)]=0
k[restindex]=apply(ad1,1,sum1)
NoSamplesAvailable=apply(!is.na(datE[,restindex]),2,sum)
k[restindex][NoSamplesAvailable< MinimumNoSamples]=NA
} # end of if
} # end of else statement
k
} # end of function




# ===================================================
# The function PickHardThreshold can help one to estimate the cut-off value
# when using the signum (step) function.
# The first column lists the threshold ("cut"),
# the second column lists the corresponding p-value based on the Fisher Transform
# of the correlation.
# The third column reports the resulting scale free topology fitting index R^2.
# The fourth column reports the slope of the fitting line, it shoud be negative for
# biologically meaningul networks.
# The fifth column reports the fitting index for the truncated exponential model.
# Usually we ignore it.
# The remaining columns list the mean, median and maximum resulting connectivity.
# To pick a hard threshold (cut) with the scale free topology criterion:
# aim for high scale free R^2 (column 3), high connectivity (col 6) and negative slope
# (around -1, col 4).
# The output is a list with 2 components. The first component lists a sugggested cut-off
# while the second component contains the whole table.
# The removeFirst option removes the first point (k=0, P(k=0)) from the regression fit.
# no.breaks specifies how many intervals used to estimate the frequency p(k) i.e. the no. of points in the
# scale free topology plot.
if (exists("PickHardThreshold")) rm(PickHardThreshold);
PickHardThreshold=function(datExpr1,RsquaredCut=0.85, cutvector=seq(0.1,0.9,by=0.05) ,removeFirst=FALSE,no.breaks=10) {
no.genes   <- dim(datExpr1)[[2]]
no.genes <- dim(datExpr1)[[2]]
no.samples= dim(datExpr1)[[1]]
colname1=c("Cut","p-value", "scale law R^2", "slope="  ,"truncated R^2","mean(k)","median(k)","max(k)")
datout=data.frame(matrix(666,nrow=length(cutvector),ncol=length(colname1) ))
names(datout)=colname1
datout[,1]=cutvector
for (i in c(1:length(cutvector) ) ){
cut1=cutvector[i]
datout[i,2]=2*(1-pt(sqrt(no.samples-1)*cut1/sqrt(1-cut1^2),no.samples-1))}
if(exists("fun1")) rm(fun1)
fun1=function(x) {
corx=abs(cor(x,datExpr1,use="p"))
out1=rep(NA, length(cutvector) )
for (j in c(1:length(cutvector))) {out1[j]=sum(corx>cutvector[j])}
out1
} # end of fun1
datk=t(apply(datExpr1,2,fun1))
for (i in c(1:length(cutvector) ) ){
nolinkshelp <- datk[,i]-1
cut2=cut(nolinkshelp,no.breaks)
binned.k=tapply(nolinkshelp,cut2,mean)
freq1=as.vector(tapply(nolinkshelp,cut2,length)/length(nolinkshelp))
# The following code corrects for missing values etc
breaks1=seq(from=min(nolinkshelp),to=max(nolinkshelp),length=no.breaks+1)
hist1=hist(nolinkshelp,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
binned.k2=hist1$mids
binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
binned.k=ifelse(binned.k==0,binned.k2,binned.k)
freq1=ifelse(is.na(freq1),0,freq1)
xx= as.vector(log10(binned.k))
if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
plot(xx,log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))" )
lm1= lm(as.numeric(log10(freq1+.000000001))~ xx )
lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) )
datout[i,3]=summary(lm1)$adj.r.squared
datout[i,4]=summary(lm1)$coefficients[2,1]
datout[i,5]=summary(lm2)$adj.r.squared
datout[i,6]=mean(nolinkshelp)
datout[i,7]= median(nolinkshelp)
datout[i,8]= max(nolinkshelp)
}
datout=signif(datout,3)
print(data.frame(datout));
# the cut-off is chosen as smallest cut with R^2>RsquaredCut
ind1=datout[,3]>RsquaredCut
indcut=NA
indcut=ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
# this estimates the cut-off value that should be used.
# Don't trust it. You need to consider slope and mean connectivity as well!
cut.estimate=cutvector[indcut][[1]]
list(cut.estimate, data.frame(datout));
} # end of function











# ===========================================================
# The function PickSoftThreshold allows one to estimate the power parameter when using
# a soft thresholding approach with the use of the power function AF(s)=s^Power
# The function PickSoftThreshold allows one to estimate the power parameter when using
# a soft thresholding approach with the use of the power function AF(s)=s^Power
# The removeFirst option removes the first point (k=1, P(k=1)) from the regression fit.
if (exists("PickSoftThreshold")) rm(PickSoftThreshold);
PickSoftThreshold=function(datExpr1,RsquaredCut=0.85, powervector=c(seq(1,10,by=1),seq(12,20,by=2)),
removeFirst=FALSE,no.breaks=10) {
no.genes <- dim(datExpr1)[[2]]
no.samples= dim(datExpr1)[[1]]
colname1=c("Power", "scale law R^2" ,"slope", "truncated R^2","mean(k)","median(k)","max(k)")
datout=data.frame(matrix(666,nrow=length(powervector),ncol=length(colname1) ))
names(datout)=colname1
datout[,1]=powervector
if(exists("fun1")) rm(fun1)
fun1=function(x) {
corx=abs(cor(x,datExpr1,use="p"))
out1=rep(NA, length(powervector) )
for (j in c(1:length(powervector))) {out1[j]=sum(corx^powervector[j])}
out1
} # end of fun1
datk=t(apply(datExpr1,2,fun1))
for (i in c(1:length(powervector) ) ){
nolinkshelp <- datk[,i]-1
cut2=cut(nolinkshelp,no.breaks)
binned.k=tapply(nolinkshelp,cut2,mean)
freq1=as.vector(tapply(nolinkshelp,cut2,length)/length(nolinkshelp))
# The following code corrects for missing values etc
breaks1=seq(from=min(nolinkshelp),to=max(nolinkshelp),length=no.breaks+1)
hist1=hist(nolinkshelp,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
binned.k2=hist1$mids
binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
binned.k=ifelse(binned.k==0,binned.k2,binned.k)
freq1=ifelse(is.na(freq1),0,freq1)

xx= as.vector(log10(binned.k))
if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
plot(xx,log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))" )
lm1= lm(as.numeric(log10(freq1+.000000001))~ xx )
lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) )
datout[i,2]=summary(lm1)$adj.r.squared
datout[i,3]=summary(lm1)$coefficients[2,1]
datout[i,4]=summary(lm2)$adj.r.squared
datout[i,5]=mean(nolinkshelp)
datout[i,6]= median(nolinkshelp)
datout[i,7]= max(nolinkshelp)
}
datout=signif(datout,3)
print(data.frame(datout));
# the cut-off is chosen as smallest cut with R^2>RsquaredCut
ind1=datout[,2]>RsquaredCut
indcut=NA
indcut=ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
# this estimates the power value that should be used.
# Don't trust it. You need to consider slope and mean connectivity as well!
power.estimate=powervector[indcut][[1]]
list(power.estimate, data.frame(datout));
}






# ===================================================
# The function ScaleFreePlot1 creates a plot for checking scale free topology
# when truncated1=T is specificed, it provides the R^2 measures for the following
# degree distributions: a) scale free topology, b) log-log R^2 and c) truncated exponential R^2

# The function ScaleFreePlot1 creates a plot for checking scale free topology
if(exists("ScaleFreePlot1")) rm(ScaleFreePlot1) ;

ScaleFreePlot1=function(kk,no.breaks=10,AF1="" ,truncated1=FALSE, removeFirst=FALSE,cex.lab1=1){

#bin data into no.breaks bins: first create the factor cut1 and code the values
cut1=cut(kk,no.breaks)
#now calculate the mean of each bin
binned.k=tapply(kk,cut1,mean)
freq1=tapply(kk,cut1,length)/length(kk)
# The following code corrects for missing values etc
breaks1=seq(from=min(kk, na.rm=T),to=max(kk, na.rm=T),length=no.breaks+1)
hist1=hist(kk,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
binned.k2=hist1$mids
binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
binned.k=ifelse(binned.k==0,binned.k2,binned.k)
freq1=ifelse(is.na(freq1),0,freq1)
plot(log10(binned.k),log10(freq1+.000000001),xlab=paste(AF1,"log10(k)"),ylab="log10(p(k))",cex.lab=cex.lab1 )
xx= as.vector(log10(binned.k))
if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
lm1=lm(as.numeric(log10(freq1+.000000001))~ xx )
lines(xx,predict(lm1),col=1)
OUTPUT=data.frame(ScaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),Slope=round(lm1$coefficients[[2]],2))
if (truncated1==TRUE) {
lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) );
OUTPUT=data.frame(ScaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),Slope=round(lm1$coefficients[[2]],2),
TruncatedRsquared=round(summary(lm2)$adj.r.squared,2))
print("the red line corresponds to the truncated exponential fit")
lines(xx,predict(lm2),col=2);
title(paste(
"scale free R^2=",as.character(round(summary(lm1)$adj.r.squared,2)),
", slope=", round(lm1$coefficients[[2]],2),
", trunc.R^2=",as.character(round(summary(lm2)$adj.r.squared,2))))} else {
title(paste("R^2=",as.character(round(summary(lm1)$adj.r.squared,2)) ,
" sl=", round(lm1$coefficients[[2]],2)), cex=0.4)
}
OUTPUT
} # end of function









#################################################################################################################
################################################################################################################################
# B) Computing the topological overlap matrix
#################################################################################################################
#################################################################################################################



# ===================================================
#The function TOMdist1 computes a dissimilarity
# based on the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
if(exists("TOMdist1")) rm(TOMdist1);

TOMdist1=function(adjmat1, maxADJ=FALSE) {
diag(adjmat1)=0;
adjmat1[is.na(adjmat1)]=0;
maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) );
if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else {
if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>0   ) {print("ERROR: non-symmetric adjacency matrix!!!") } else {
kk=apply(adjmat1,2,sum)
maxADJconst=1
if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 )))
Dhelp1=matrix(kk,ncol=length(kk),nrow=length(kk))
denomTOM= pmin(as.dist(Dhelp1),as.dist(t(Dhelp1)))   +as.dist(maxADJconst-adjmat1);
gc();gc();
numTOM=as.dist(adjmat1 %*% adjmat1 +adjmat1);
#TOMmatrix=numTOM/denomTOM
# this turns the TOM matrix into a dissimilarity
out1=1-as.matrix(numTOM/denomTOM)
diag(out1)=1
out1
}}
}

#---------------------------------------------------------------------------
# This is a somewhat modified TOMdist1.

SignedTOMdist = function(adjmat1, maxADJ=FALSE)
{
  diag(adjmat1)=0;
  adjmat1[is.na(adjmat1)]=0;
  collect_garbage();
  kk=apply(adjmat1,2,sum)
  collect_garbage();
  maxADJconst=1
  if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 )))
  collect_garbage();
  Dhelp1 = matrix(kk,ncol=length(kk),nrow=length(kk))
  collect_garbage();
  denomTOM = pmin(as.dist(Dhelp1),as.dist(t(Dhelp1))) + as.dist(maxADJconst-adjmat1);
  rm(Dhelp1);
  collect_garbage();
  gc(); gc();
  numTOM=as.dist(adjmat1 %*% adjmat1 +adjmat1);
  collect_garbage();
  #TOMmatrix=numTOM/denomTOM
  # this turns the TOM matrix into a dissimilarity
  out1=1-as.matrix(numTOM/denomTOM)
  rm(numTOM); rm(denomTOM);
  collect_garbage();
  diag(out1)=1
  out1
}



# ===================================================
# This function computes a TOMk dissimilarity
# which generalizes the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
# WARNING:  ONLY FOR UNWEIGHTED NETWORKS, i.e. the adjacency matrix contains binary entries...
# This function is explained in Yip and Horvath (2005)
# http://www.genetics.ucla.edu/labs/horvath/TOM/
if(exists("TOMkdist1")) rm(TOMkdist1);
TOMkdist1 = function(adjmat1,k=1){
    maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) );
    if (k!=round(abs(k))) {
        stop("k must be a positive integer!!!", call.=TRUE);}
    if (maxh1>1 | minh1 < 0 ){
        print(paste("ERROR: entries of the adjacency matrix must be between inclusively 0 and 1!!!, max=",maxh1,", min=",minh1))}
    else {
if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>0   ) {print("ERROR: non-symmetric adjacency matrix!!!") } else {

        B <- adjmat1;
        if (k>=2) {
            for (i in 2:k) {
                diag(B) <- diag(B) + 1;
                B = B %*% adjmat1;}}   # this gives the number of paths with length at most k connecting a pair
        B <- (B>0);   # this gives the k-step reachability from a node to another
        diag(B) <- 0;   # exclude each node being its own neighbor
        B <- B %*% B   # this gives the number of common k-step-neighbor that a pair of nodes share

        Nk <- diag(B);
        B <- B +adjmat1;   # numerator
        diag(B) <- 1;
        denomTOM=outer(Nk,Nk,FUN="pmin")+1-adjmat1;
        diag(denomTOM) <- 1;
        1 - B/denomTOM   # this turns the TOM matrix into a dissimilarity
}}
}


# IGNORE THIS function...
# The function TOMdistROW computes the TOM distance of a gene (node)
# with that of all other genes in the network.
# WhichRow is an integer that specifies which row of the adjacency matrix
# corresponds to the gene of interest.
# Output=vector of TOM distances.
if (exists("TOMdistROW") ) rm(TOMdistROW)
TOMdistROW=function(WhichRow=1, adjmat1, maxADJ=FALSE) {
diag(adjmat1)=0;
maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) );
if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else {
kk=apply(adjmat1,2,sum)
numTOM=adjmat1[WhichRow,] %*% adjmat1 +adjmat1[WhichRow,];
numTOM[WhichRow]=1
maxADJconst=1
if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 )))
denomTOM=pmin(kk[WhichRow],kk)+maxADJconst-adjmat1[WhichRow,]; denomTOM[WhichRow]=1
#TOMmatrix=numTOM/denomTOM
# this turns the TOM matrix into a dissimilarity
1-numTOM/denomTOM
}
}


#####################################################################################################
################################################################################################################################
# C) Defining gene modules using clustering procedures
#####################################################################################################
################################################################################################################################

# ===================================================
#The function modulecolor2 function assigns colors to the observations
# in the branches of a dendrogram
# we use it to define modules....
if (exists("modulecolor2")) rm(modulecolor2);
modulecolor2=function(hier1, h1=0.9,minsize1=50) {
# here we define modules by using a height cut-off for the branches
labelpred= cutree(hier1,h=h1)
sort1=-sort(-table(labelpred))
modulename= as.numeric(names(sort1))
modulebranch= sort1>minsize1
no.modules=sum(modulebranch)
# now we assume that there are fewer than a certain number of colors
#colorcode=c("turquoise","blue","brown","yellow","green","red","black","purple","orange","pink",
#"greenyellow","lightcyan","salmon","midnightblue","lightyellow")
colorcode=c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
                "purple","greenyellow","tan","salmon", "midnightblue", "lightcyan","grey60",
                "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise",
                "darkgrey", "orange", "darkorange", "white" )

# "grey" means not in any module;
colorhelp=rep("grey",length(labelpred))
if ( no.modules==0 | no.modules >length(colorcode)){ print(paste("The number of modules is problematic! \n Number of modules = ", as.character(no.modules)))} else { for (i in c(1:no.modules)) {colorhelp=ifelse(labelpred==modulename[i],colorcode[i],colorhelp)};
colorhelp=factor(colorhelp,levels=c(colorcode[1:no.modules],"grey"))
}
factor(colorhelp, levels=unique(colorhelp[hier1$order] ))
}


#---------------------------------------------------------------------------------
#
# ModuleNumber
#
#---------------------------------------------------------------------------------
# Similar to modulecolor2 above, but returns numbers instead of colors, which is oftentimes more useful.
# 0 means unassigned.
# Return value is a simple vector, not a factor.
# Caution: the module numbers are neither sorted nor sequential, the only guarranteed fact is that grey
# probes are labeled by 0 and all probes belonging to the same module have the same number.
# If size-sorted sequential labels are required, "normalize" the result by calling NormalizeLabels
# (below).

ModuleNumber = function(HierTree, CutHeight = 0.9, MinSize = 50)
{
  Branches = cutree(HierTree, h = CutHeight);
  NOnBranches = table(Branches);
  #NOnBranches[i][[1]] somehow gives the number of elements on branch i.
  TrueBranch = NOnBranches >= MinSize;
  Branches[!TrueBranch[Branches]] = 0;

  #NewLabels = levels(factor(Branches));
  #for (lab in 1:length(NewLabels)) if (NewLabels[lab]!=0)
  #  Branches[Branches==NewLabels[lab]] = lab;

  Branches;

}


# The function hclustplot1 creates a barplot where the colors of the bars are sorted according to
# a hierarchical clustering tree (hclust object)
#if (exists("hclustplot1")) rm(hclustplot1);
#hclustplot1=function(hier1,couleur,title1="Colors sorted by hierarchical clustering")
#{
#if (length(hier1$order) != length(couleur) ) {print("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree")};
#if (length(hier1$order) == length(couleur) ) {
#barplot(height=rep(1, length(couleur)), col= as.character(couleur[hier1$order]),border=F, main=title1,space=0, axes=F)}
#}

if (exists("hclustplot1")) rm(hclustplot1);
hclustplot1=function(hier1,Color1, Color2=NULL,title1="Colors sorted by hierarchical clustering")
{
  options(stringsAsFactors=FALSE);
  if (length(hier1$order) != length(Color1) )
  {
    stop("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree");
  } else {
    if (is.null(Color2))
    {
       barplot(height=rep(1, length(Color1)), col= as.character(Color1[hier1$order]),
               border=F, main=title1,space=0, axes=F)
    } else if (length(Color1)==length(Color2)) {
       # height = matrix(0.5, nrow = 2, ncol = length(Color1));
       C1 = Color1[hier1$order]; C2 = Color2[hier1$order]
       step = 1/length(Color1);
       barplot(height=1, col = "white", border=F, main=title1,space=0, axes=F)
       for (i in 1:(length(Color1)))
       {
         lines(x=rep((i*step), times=2), y=c(0,0.5),  col = as.character(C1[i]));
         lines(x=rep((i*step), times=2), y=c(0.5,1),  col = as.character(C2[i]));
       }
    }
  }
}

if (exists("hclustplotn")) rm(hclustplotn);
hclustplotn=function(hier1, Color, RowLabels=NULL, cex.RowLabels = 0.9, ...)
{
  options(stringsAsFactors=FALSE);
  if (length(hier1$order) != dim(Color)[[1]] )
  {
    stop("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree");
  } else {
       No.Sets = dim(Color)[[2]];
       C = Color[hier1$order, ];
       step = 1/dim(Color)[[1]];
       ystep = 1/No.Sets;
       barplot(height=1, col = "white", border=F,space=0, axes=F, ...)
       for (j in 1:No.Sets)
       {
         ind = (1:(dim(C)[1]));
         xl = (ind-1) * step; xr = ind * step;
         yb = rep(ystep*(j-1), dim(C)[1]); yt = rep(ystep*j, dim(C)[1]);
         rect(xl, yb, xr, yt, col = as.character(C[,j]), border = as.character(C[,j]));
         if (is.null(RowLabels))
         {
             text(as.character(j), pos=2, x=0, y=ystep*(j-0.5), cex=cex.RowLabels, xpd = TRUE);
         } else {
             text(RowLabels[j], pos=2, x=0, y=ystep*(j-0.5), cex=cex.RowLabels, xpd = TRUE);
         }
       }
       for (j in 1:No.Sets) lines(x=c(0,1), y=c(ystep*j,ystep*j));
  }
}


# ===================================================
# The function TOMplotn creates a TOM plot
# Inputs:  distance measure, hierarchical (hclust) object, color label=couleur
if (exists("TOMplotn")) rm(TOMplotn);
TOMplot1=function(disttom,hier1, couleur,terrainColors=FALSE) {
no.nodes=length(couleur)
if (no.nodes != dim(disttom)[[1]] ) {print("ERROR: number of color labels does not equal number of nodes in disttom")} else {
labeltree=as.character(couleur)
labelrow  = labeltree
labelrow[hier1$order[length(labeltree):1]]=labelrow[hier1$order]
options(expressions = 10000)
if (terrainColors) heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none",revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
labRow=F, labCol=F, col = terrain.colors(1000)) else heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none",revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
labRow=F, labCol=F)
}
} #end of function


# ===================================================
# The function TOMplot2 creates a TOM plot where the top and left color bars can be different
# Inputs:  distance measure, hierarchical (hclust) object, color label=couleurTop, couleurLeft
if (exists("TOMplot2")) rm(TOMplot2);
TOMplot2=function(disttom,hier1, couleurTop, couleurLeft) {
no.nodes=length(couleurTop)
if (no.nodes != length(couleurLeft)) {stop("ERROR: number of top color labels does not equal number of left color labels")}
if (no.nodes != dim(disttom)[[1]] ) {stop("ERROR: number of color labels does not equal number of nodes in disttom")} else {
labeltree = as.character(couleurTop)
labelrow  = as.character(couleurLeft)
labelrow[hier1$order[length(labeltree):1]]=labelrow[hier1$order]
options(expressions = 10000)
heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none", revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
labRow=F, labCol=F)
}
} #end of function



# IGNORE THIS FUNCTION...
# The function "BestHeightCut" allows one to find the best height cut-off
# for a hierarchical clustering tree when external gene information is available
# It computes a Kruskal Wallis-test p-value for each height cut-off
# based on determining whether gene significance differs across branch membership.
if(exists("BestHeightCut")) rm(BestHeightCut);
BestHeightCut=function(hier1, GeneSignif, hcut=seq(0.1,.95,by=0.01) ) {
pvalues=rep(NA, length(hcut))
for (i in c(1:length(hcut))) {
colorhelp=modulecolor2(hier1,hcut[i])
if (length(unique(colorhelp))>1 ) {pvalues[i]=kruskal.test(GeneSignif, colorhelp)$p.value}
data.frame(hcut,pvalues)
}}




#####################################################################################################
################################################################################################################################
# D) Summing up modules using their first principal components (first eigengene)
#####################################################################################################
################################################################################################################################

# ===================================================
#The function ModulePrinComps1 finds the first principal component (eigengene) in each
# module defined by the colors of the input vector "couleur" (Pardon my French).
# It also reports the variances explained by the first 5 principal components.
# And it yields a measure of module conformity for each gene,
# which is highly correlated to the within module connectivity.
# The theoretical underpinnings are described in Horvath, Dong, Yip (2005)
# http://www.genetics.ucla.edu/labs/horvath/ModuleConformity/
# This requires the R library impute
# The output is a list with 3 components:
# 1) a data frame of module eigengenes (MEs),
# 2) a data frame that lists the percent variance explained by the first 5 MEs of a module
# 3) a data frame that lists the module conformity for each gene.
# The be used as alternative connectivity measure....
if(exists("ModulePrinComps1")) rm(ModulePrinComps1);
ModulePrinComps1=function(datexpr, couleur, verbose = 1, print.level = 0, Impute = FALSE,
                          GetConformity = TRUE) {
if (is.null(datexpr))
{
  print("ModulePrinComps1: Error: datexpr is NULL. ");
  stop();
}
if (is.null(couleur))
{
  print("ModulePrinComps1: Error: couleur is NULL. ");
  stop()
}
MaxVectors = 5;
#print(paste("datexpr dimensions:", as.character(dim(datexpr))));
spaces = PrintSpaces(print.level);
modlevels=levels(factor(couleur))
PrinComps=data.frame(matrix(666,nrow=dim(datexpr)[[1]],ncol= length(modlevels)))
varexplained= data.frame(matrix(666,nrow= 5,ncol= length(modlevels)))
names(PrinComps)=paste("ME",modlevels,sep="")
for(i in c(1:length(modlevels)) )
{
  if (verbose>0)
    print.flush(paste(spaces, "ModulePrinComps1 : Working on ME for module ", modlevels[i], sep = ""));
  modulename    = modlevels[i]
  restrict1= as.character(couleur)== modulename
  #print(paste("length of couleur:", length(couleur), "; length of restric1:", length(restrict1)));
   # in the following, rows are genes and columns are samples
  datModule=t(datexpr[, restrict1])
  if (Impute)
  {
      saved.seed = .Random.seed;
      datModule=impute.knn(as.matrix(datModule))
      datModule=t(scale(t(datModule)));
      .Random.seed = saved.seed;
  }
  n = dim(datModule)[1]; p = dim(datModule)[2];
  svd1=svd(datModule, nu = min(n, p, MaxVectors), nv = min(n, p, MaxVectors));
  mtitle=paste("MEs of ", modulename," module", sep="");
  varexplained[,i]= (svd1$d[1:5])^2/sum(svd1$d^2)
  # this is the first principal component
      pc1=svd1$v[,1]
  # signh1=sign(sum(cor(pc1,  t(datModule))))
  # if (signh1 != 0)  pc1=signh1* pc1
  PrinComps[,i]= pc1
}
ModuleConformity= rep(666,length=dim(datexpr)[[2]])
if (GetConformity)
{
 for(i in 1:(dim(datexpr)[[2]]))
     ModuleConformity[i] = abs(cor(datexpr[,i], PrinComps[,match(couleur[i], modlevels)],
                                   use="pairwise.complete.obs"))
} else
{
   ModuleConformity = NULL;
}

list(PrinComps=PrinComps, varexplained=varexplained, ModuleConformity=ModuleConformity)
}



#####################################################################################################
################################################################################################################################
# E) Relating a measure of gene significance to the modules
#####################################################################################################
################################################################################################################################

# ===================================================
# The function ModuleEnrichment1 creates a bar plot that shows whether modules are enriched with
# significant genes.
# More specifically, it reports the mean gene significance for each module.
# The gene significance can be a binary variable or a quantitative variable.
# It also plots the 95% confidence interval of the mean (CI=mean +/- 1.96* standard error).
# It also reports a Kruskal Wallis P-value.
if( exists("ModuleEnrichment1") ) rm(ModuleEnrichment1);
ModuleEnrichment1=function(genesignif1,couleur,title1="gene significance across modules",labely="Gene Significance",boxplot=F) {
if (length(genesignif1) != length(couleur) ) print("Error: vectors don\'t have the same lengths") else {
if (boxplot != TRUE) {
mean1=function(x) mean(x,na.rm=T)
means1=as.vector(tapply(genesignif1,couleur,mean1));
se1= as.vector(tapply(genesignif1,couleur,stderr1))
#par(mfrow=c(1,1))
barplot(means1,
names.arg=names(table(couleur) ),col= names(table(couleur) )
,ylab=labely)
err.bp(as.vector(means1), as.vector(1.96*se1), two.side=T)} else {
boxplot(split(genesignif1,couleur),notch=T,varwidth=T, col= names(table(couleur) ),ylab=labely)}

title(paste(title1,", p-value=", signif(kruskal.test(genesignif1,factor(couleur))$p.value,2)))
}
} # end of function


# IGNORE THIS...
# ===================================================
#The function fisherPvector allows one to compute Fisher exact p-values
# Thus it allows one to carry out an EASE analysis
# Output: a table of Fisher��s exact p-values
# Input: annotation1 is a vector of gene annotations
# Input: couleur (French for color) denotes the module color of each gene
# Only those gene functions (levels of annotation1) that occur a certain mininum number of times
# (parameter= minNumberAnnotation) in the data set will be considered.
if (exists("fisherPvector" ) ) rm(fisherPvector);
fisherPvector=function(couleur,annotation1,minNumberAnnotation=50) {
levelsannotation1=levels(annotation1)
levelscouleur=levels(factor(couleur))
no.couleur=length(levelscouleur)
restannotation1=table(annotation1)>minNumberAnnotation
no.annotation=sum( restannotation1)
datoutP=data.frame(matrix(666,nrow=no.annotation,ncol=no.couleur) )
#datoutProp=data.frame(matrix(666,nrow=no.annotation,ncol=2*no.couleur) )
#names(datoutProp)=paste("Prop",paste( rep(levelscouleur ,rep(2, length(levelscouleur))) ) , c("Y","N")  ,sep=".")
datoutProp=data.frame(matrix(666,nrow=no.annotation,ncol=no.couleur) )
names(datoutProp)=paste("Perc",levelscouleur , sep=".")
names(datoutP)=paste("p",levelscouleur,sep=".")
restlevelsannotation1= levelsannotation1[restannotation1]
row.names(datoutP)= restlevelsannotation1
for (i in c(1:no.annotation) ) {
for (j in c(1:no.couleur) ){
tab1=table( annotation1 !=restlevelsannotation1[i], couleur !=levelscouleur[j])
datoutP[i,j]=signif(fisher.test(tab1)$p.value,2)
#datoutProp[i,2*j-1]=signif(tab1[1,1]/sum(tab1[,1] ),2)
#datoutProp[i,2*j]= signif(tab1[1,2]/sum(tab1[,2]) ,2)
}
table2=table(annotation1 !=restlevelsannotation1[i], couleur)
datoutProp[i,]= signif(table2[1,]/apply(table2,2,sum),2)
}
data.frame(datoutP,datoutProp)
} # end of function fisherPvector



#####################################################################################################
################################################################################################################################
# F) Carrying out a within module analysis (computing intramodular connectivity etc)
#####################################################################################################
################################################################################################################################

# ===================================================
#The function DegreeInOut computes for each gene
#a) the total number of connections,
#b) the number of connections with genes within its module,
#c) the number of connections with genes outside its module
# When scaledToOne=TRUE, the within module connectivities are scaled to 1, i.e. the max(K.Within)=1 for each module
if (exists("DegreeInOut")) rm(DegreeInOut); DegreeInOut =function(adj1, couleur,scaledToOne=FALSE) {
no.nodes=length(couleur)
couleurlevels=levels(factor(couleur))
no.levels=length(couleurlevels)
kWithin=rep(-666,no.nodes )
diag(adj1)=0
for (i in c(1:no.levels) ) {
rest1=couleur==couleurlevels[i];
if (sum(rest1) <3 ) { kWithin[rest1]=0 } else {
kWithin[rest1]=apply(adj1[rest1,rest1],2,sum)
if (scaledToOne) kWithin[rest1]=kWithin[rest1]/max(kWithin[rest1])}
}
kTotal= apply(adj1,2,sum)
kOut=kTotal-kWithin
if (scaledToOne) kOut=NA
kDiff=kWithin-kOut
data.frame(kTotal,kWithin,kOut,kDiff)
}


# =======================================================================
# The function WithinModuleCindex1 relates the node measures (e.g. connectivities)
# to "external" node significance information within each  module,
# i.e. it  carries out a by module analysis.
# Output: first column reports the spearman correlation p-value between the network variable and the
# node significance. The next columns contain the Spearman correlations between the variables.
if (exists("WithinModuleAnalysis1")) rm(WithinModuleAnalysis1);
WithinModuleAnalysis1=function(datnetwork,nodesignif, couleur)
{
  cortesthelp=function( x ) {
    len1=dim(x)[[2]]-1
    out1=rep(666, len1);
    for (i in c(1:len1) ) {out1[i]= signif( cor.test(x[,i+1], x[,1], method="s",use="p" )$p.value ,2) }
    data.frame( variable=names(x)[-1] , NS.CorPval=out1, NS.cor=t(signif(cor (x[,1], x[,-1],use="p",method="s"),2)),
                signif(cor(x[,-1],use="p",method="s"),2) )
  } #end of function cortesthelp
  print("IGNORE  the warnings...");
  by( data.frame(nodesignif, datnetwork), couleur, cortesthelp);
} #end of function WithinModuleAnalysis


# =======================================================================
# The function WithinModuleCindex1 relates the node measures (e.g. connectivities)
# to "external" node significance information within each  module,
# i.e. it  carries out a by module analysis.
# BUT it focuses on the C-index also known as area under the ROC curve
# This measure is related to Kendall's Tau statistic and Somer's D,
# see F. Harrel (Regression Modeling Strategies). Springer.
# It requires the following library
library(Hmisc)
# Output: the first column reports the C-index and the second, p-value
if (exists("WithinModuleCindex1")) rm(WithinModuleCindex1);
WithinModuleCindex1=function(datnetwork,nodesignif, couleur) {
CindexFun=function( x ) {
len1=dim(x)[[2]]-1
outC=rep(666, len1);
outP=rep(666, len1);
for (i in c(1:len1) ) {rcor1=rcorr.cens(x[,i+1], x[,1])
outC[i]=rcor1[[1]]
outP[i]=1- pnorm(abs(rcor1[[2]]/rcor1[[3]]))
}
data.frame( variable=names(x)[-1] , C.index=outC, p.value=outP)
} #end of function CindexFun
#print("IGNORE  the warnings...");
by( data.frame(nodesignif, datnetwork),couleur,CindexFun);
} #end of function WithinModuleAnalysis


# The following function allows on to plot a gene (node) significance measure versus
# connectivity.
if(exists("plotConnectivityGeneSignif1") ) rm( plotConnectivityGeneSignif1);
plotConnectivityGeneSignif1=function(degree1,genesignif1,color1="black",
title1="Gene Significance vs Connectivity" , xlab1="Connectivity", ylab1="GeneSignificance") {
lm1=lm(genesignif1~degree1 ,na.action="na.omit")
plot(degree1, genesignif1, col=color1,ylab=ylab1,xlab=xlab1,main=paste(title1, ", cor=",
signif(cor( genesignif1,degree1, method="s",use="p" )   ,2) ))
abline(lm1)
}




#####################################################################################################
################################################################################################################################
# G) Miscellaneous other functions, e.g. for computing the cluster coefficient.
#####################################################################################################
################################################################################################################################



# ===================================================
# The function ClusterCoef.fun computes the cluster coefficients.
# Input is an adjacency matrix
if(exists("ClusterCoef.fun")) rm(ClusterCoef.fun) ; ClusterCoef.fun=function(adjmat1) {
diag(adjmat1)=0
no.nodes=dim(adjmat1)[[1]]
computeLinksInNeighbors <- function(x, imatrix){x %*% imatrix %*% x}
nolinksNeighbors <- c(rep(-666,no.nodes))
total.edge <- c(rep(-666,no.nodes))
maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) );
if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else {
nolinksNeighbors <- apply(adjmat1, 1, computeLinksInNeighbors, imatrix=adjmat1)
plainsum  <- apply(adjmat1, 1, sum)
squaresum <- apply(adjmat1^2, 1, sum)
total.edge = plainsum^2 - squaresum
CChelp=rep(-666, no.nodes)
CChelp=ifelse(total.edge==0,0, nolinksNeighbors/total.edge)
CChelp}
} # end of function



# ===================================================
# The function err.bp  is used to create error bars in a barplot
# usage: err.bp(as.vector(means), as.vector(stderrs), two.side=F)

err.bp<-function(daten,error,two.side=F){
 if(!is.numeric(daten)) {
      stop("All arguments must be numeric")}
 if(is.vector(daten)){
    xval<-(cumsum(c(0.7,rep(1.2,length(daten)-1))))
 }else{
    if (is.matrix(daten)){
      xval<-cumsum(array(c(1,rep(0,dim(daten)[1]-1)),
dim=c(1,length(daten))))+0:(length(daten)-1)+.5
    }else{
      stop("First argument must either be a vector or a matrix") }
 }
 MW<-0.25*(max(xval)/length(xval))
 ERR1<-daten+error
 ERR2<-daten-error
 for(i in 1:length(daten)){
    segments(xval[i],daten[i],xval[i],ERR1[i])
    segments(xval[i]-MW,ERR1[i],xval[i]+MW,ERR1[i])
    if(two.side){
      segments(xval[i],daten[i],xval[i],ERR2[i])
      segments(xval[i]-MW,ERR2[i],xval[i]+MW,ERR2[i])
    }
 }
}

# ===================================================
# this function computes the standard error
if (exists("stderr1")) rm(stderr1)
stderr1 <- function(x){ sqrt( var(x,na.rm=T)/sum(!is.na(x))   ) }



# ===================================================
# The following two functions are for displaying the pair-wise correlation in a panel when using the command "pairs()"
# Typically, we use "pairs(DATA, upper.panel=panel.smooth, lower.panel=panel.cor, diag.panel=panel.hist)" to
# put the correlation coefficients on the lower panel.
panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * r)
}
panel.hist <- function(x, ...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}



# ===================================================
# this function computes the standard error
if (exists("stderr1")) rm(stderr1);
stderr1 <- function(x){ sqrt( var(x,na.rm=T)/sum(!is.na(x))   ) }




# ===================================================
# This function collects garbage
if (exists("collect_garbage")) rm(collect_garbage);
collect_garbage=function(){while (gc()[2,4] != gc()[2,4] | gc()[1,4] != gc()[1,4]){}}
collect_garbage()


# this function is used for computing the Rand index below...
# ===================================================
if (exists("choosenew") ) rm(choosenew)
choosenew <- function(n,k){
  n <- c(n)
  out1 <- rep(0,length(n))
  for (i in c(1:length(n)) ){
    if (n[i]<k) {out1[i] <- 0}
    else {out1[i] <- choose(n[i], k)}}
  out1
}


# ===================================================
# the following function computes the Rand index between 2 clusterings
# assesses how similar two clusterings are
if (exists("Rand1") ) rm(Rand1)
Rand2 <- function(tab,adjust=T) {
  a <- 0; b <- 0; c <- 0; d <- 0; nn <- 0
  m <- nrow(tab);
  n <- ncol(tab);
  for (i in 1:m) {
    c<-0
    for(j in 1:n) {
      a <- a+choosenew(tab[i,j],2)
      nj <- sum(tab[,j])
      c <- c+choosenew(nj,2)
    }
    ni <- sum(tab[i,])
    b <- b+choosenew(ni,2)
    nn <- nn+ni
  }
  if(adjust==T) {
    d <- choosenew(nn,2)
    adrand <- (a-(b*c)/d)/(0.5*(b+c)-(b*c)/d)
    adrand
  } else {
    b <- b-a
    c <- c-a
    d <- choosenew(nn,2)-a-b-c
    rand <- (a+d)/(a+b+c+d)
    rand
  }
}

# ===================================================
# This function is used in "pairs()" function. The problem of the original  panel.cor is that
# when the correlation coefficient is very small, the lower panel will have a large font
# instead of a mini-font in a saved .ps file. This new function uses a format for corr=0.2
# when corr<0.2, but it still reports the original value of corr, with a minimum format.

panel.cor1=function(x, y, digits=2, prefix="", cex.cor){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    txt1=txt
    r1=r
    if (r<0.2) {
        r1=0.2
        txt1 <- format(c(r1, 0.123456789), digits=digits)[1]
        txt1 <- paste(prefix, txt1, sep="")
        }
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt1)
    cex = cex * r1
    r <- round(r, digits)
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    text(0.5, 0.5, txt, cex=cex)
}

merge2Clusters = function(mycolorcode, mainclusterColor, minorclusterColor){
  mycolorcode2 = ifelse(as.character(mycolorcode)==minorclusterColor, mainclusterColor, as.character(mycolorcode) )
  fcolorcode   =factor(mycolorcode2)
  fcolorcode
}



###############################################################################################################################
# I) GENERAL STATISTICAL FUNCTIONS

mean1=function(x) mean(x,na.rm=T)
var1=function(x) var(x,na.rm=T)



if (exists("scatterplot1") ) rm(scatterplot1);
scatterplot1=function(x,y, title = "", ... ){
cor1=signif(cor(x,y,use="p",method="s"),2)
corp=signif(cor.test(x,y,use="p",method="s")$p.value,2)
if (corp<10^(-20) ) corp="<10^{-20}"
plot(x,y, main=paste(title, "cor=", cor1,"p=",corp),...)
}


# Functions for Network Screening



if(exists("HubGeneSignificance") ) rm(HubGeneSignificance)
HubGeneSignificance=function(k,GS,module,NumberHubs=10){
colorlevels=levels(factor(module))
HGS=rep(NA, length(colorlevels) )
for (i in c(1:length(colorlevels) ) ) {
restmodule= module==colorlevels[i]
GSModule= GS[restmodule]
kModule=k[restmodule]
HGS[i]=mean(GSModule[rank(-kModule,  ties.method="first")<=NumberHubs],na.rm=T)
} # end of for loop
barplot(HGS,col=colorlevels, ylab="Mean Hub Gene Significance", main=paste("Hub Gene Significance, top", NumberHubs, "hubs"))
datout=data.frame(matrix(HGS,nrow=1,ncol=length(colorlevels)));
names(datout)=colorlevels;
datout
} # end of function




if(exists("StandardScreening1") ) rm(StandardScreening1);
StandardScreening1=function( GS, LN=c(5,10)  ) {
datout=data.frame(matrix(F, nrow=length(GS),ncol=length(LN) ))
names(datout)=paste("List", LN,sep="")
for ( j in c(1:length(LN))  ) {
datout[,j]=rank(-GS,ties.method="first")<=LN[j]
}
datout
} # end of function StandardScreening1




if(exists("NetworkScreening1") ) rm(NetworkScreening1)
NetworkScreening1=function(datE, GS, MLN=1000, LN=10,beta=6,  minModuleSize=100,powerAllocation=3, NumberHubsHGS=50 , excludegrey=F,excludeturquoise=F, consider.sign.corKGS=F) {
if (length(GS) != dim(datE)[[2]] ) print("Error: length(GS) not compatible with datE. Please check whether length(GS) = dim(datE)[[2]]?")
if ( max(LN)>MLN ) print("Error: requested list number bigger than maximum list number. Please increase MLN or decrease LN")
if ( MLN<minModuleSize) print("Warning: maximum list number smaller than minModuleSize-->every gene is grouped into the turquoise module")
if (length(GS) == dim(datE)[[2]]  &  max(LN) <= MLN  ){
GS=abs(GS)
restMLN=rank(-GS, ties.method="first")<=MLN
GSrest=GS[restMLN]
ADJ= abs(cor(datE[, restMLN], use="p"))^beta
h1=hclust( as.dist(TOMdist1(ADJ)) , method="average")
colorGS=rep("turquoise", dim(ADJ)[[2]] )
if (MLN > minModuleSize ) {colorGS= factor(cutreeDynamic(h1, minModuleSize= minModuleSize)) }
colorGSlevels=levels(factor(colorGS))
if ( length(colorGSlevels)==1   ) { HGS=1; k=apply(ADJ,2,sum) } else {
k= DegreeInOut(ADJ ,colorGS,scaledToOne=FALSE)$kWithin
HGS=as.vector(as.matrix(HubGeneSignificance(k=k, GS=GSrest,module=colorGS,NumberHubs=NumberHubsHGS)[1,]))
# the following definition of HGS uses the correlation between k and GS...
if ( consider.sign.corKGS==T ) {
for (i in c(1:length(colorGSlevels))){
HGS[i]=(cor( k[colorGS==colorGSlevels[i]], GSrest[colorGS==colorGSlevels[i]], use="p" ))}
} # end of for loop
} # end of if (consider.sign  )
if (excludegrey) HGS[colorGSlevels=="grey"]=0
if (excludeturquoise==T) HGS[colorGSlevels=="turquoise"]=0
weightHGS=abs(HGS)^powerAllocation
weightHGS=weightHGS/sum(weightHGS)
datout=data.frame(matrix(F, nrow=length(GS),ncol=length(LN) ))
names(datout)=paste("List", LN,sep="")
for ( j in c(1:length(LN))  ) {
LNcolor=round(LN[j]*weightHGS)
maxNoColorGS=tapply(colorGS,colorGS,length)
# if there are more hubs than the module size, then we pick genes from the next module
# with highest hub gene significance
for (ii in c(1:length(colorGSlevels))) {
ExcessNo=sum(c(LNcolor-maxNoColorGS)[LNcolor>maxNoColorGS])
indexHighestWeightHGS=rank(-weightHGS, ties.method="first")==ii
LNcolor[LNcolor>maxNoColorGS]=maxNoColorGS[LNcolor>maxNoColorGS]
LNcolor[indexHighestWeightHGS  ]= LNcolor[indexHighestWeightHGS  ]+ExcessNo
LNcolor[indexHighestWeightHGS  ]=LN[j]- sum(LNcolor[!indexHighestWeightHGS])
}
LNcolor[LNcolor>maxNoColorGS]=maxNoColorGS[LNcolor>maxNoColorGS]
PickHubs=rep(F,length(GS))
for (i in c(1:length(colorGSlevels))){
# the following takes the sign between K and GS into account
signKGS=1
if ( consider.sign.corKGS==T ) {signKGS=sign(cor( k[colorGS==colorGSlevels[i]], GSrest[colorGS==colorGSlevels[i]], use="p" ))}
PickHubs[restMLN][ colorGS==colorGSlevels[i]] =rank(-signKGS*k[colorGS==colorGSlevels[i]], ties.method="first")<=LNcolor[i]
} # end of for loop
datout[,j]=PickHubs
print(paste("For list", j, "with LN=",LN[j], "the algorithm picks the following number of hubs in each module"))
print( data.frame(colorGSlevels, modulesize=table(colorGS), HubGeneSignif=signif(HGS,2), weightHGS=signif(weightHGS,2),  LNcolor    ))
} # end of for ( j in c(1:length(LN))  )
datout
} # end of if if (length(GS) == dim(datE)[[2]]  &  max(LN) < MLN
} # end of function NetworScreening1


#==========================================================================================================
#
# Peter Langfelder's additions
#
#==========================================================================================================

# PrintFlush.R

if (exists("print.flush")) { remove(print.flush); collect_garbage(); }
print.flush = function(...)
{
  x = print(...)
  if (exists("flush.console")) x=flush.console();
}

if (exists("PrintSpaces")) { remove(PrintSpaces); collect_garbage(); }
PrintSpaces = function(print.level)
{
  if (print.level>0)
    {
      spaces = paste(" ",rep("  ", times=print.level-1), collapse="");
    } else
    {
      spaces = "";
    }
  spaces;
}


#---------------------------------------------------------------------------------------------------------
# HeatmapWithTextLabels.R
#---------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------
#
# ReverseRows = function(Matrix)
#
#--------------------------------------------------------------------------
#


ReverseRows = function(Matrix)
{
  ind = seq(from=dim(Matrix)[1], to=1, by=-1);
  Matrix[ind,];
  #Matrix
}

ReverseVector = function(Vector)
{
  ind = seq(from=length(Vector), to=1, by=-1);
  Vector[ind];
  #Vector
}

#--------------------------------------------------------------------------
#
# HeatmapWithTextLabels = function ( Matrix, xLabels, yLabels, ... ) {
#
#--------------------------------------------------------------------------
# This function plots a heatmap of the specified matrix
# and labels the x and y axes wit the given labels.
# It is assumed that the number of entries in xLabels and yLabels is consistent
# with the dimensions in.
# If ColorLabels==TRUE, the labels are not printed and instead interpreted as colors --
#  -- a simple symbol with the appropriate color is printed instead of the label.
# The x,yLabels are expected to have the form "..color" as in "MEgrey" or "PCturquoise".
# xSymbol, ySymbols are additional markers that can be placed next to color labels

HeatmapWithTextLabels = function ( Matrix, xLabels, yLabels = NULL, xSymbols = NULL, ySymbols = NULL,
                                   InvertColors=FALSE, ColorLabels = NULL, xColorLabels = FALSE,
                                   yColorLabels = FALSE,
                                   SetMargins = TRUE,
                                   colors = NULL, NumMatrix = NULL, cex.Num = NULL, cex.lab = NULL,
                                   plotLegend = TRUE, ... )
{
  if (!is.null(ColorLabels)) {xColorLabels = ColorLabels; yColorLabels = ColorLabels; }

  if (is.null(yLabels) & (!is.null(xLabels)) & (dim(Matrix)[1]==dim(Matrix)[2]))
    yLabels = xLabels;
  if (SetMargins)
  {
    if (ColorLabels)
    {
      par(mar=c(2,2,3,5)+0.2);
    } else {
      par(mar = c(7,7,3,5)+0.2);
    }
  }
  if (is.null(colors))
    #if (IncludeSign)
    #{
      #colors = GreenRedWhite(50);
    #} else {
      colors = heat.colors(30);
    #}
  if (InvertColors) colors = ReverseVector(colors);
  if (plotLegend)
  {
    image.plot(t(ReverseRows(Matrix)), xaxt = "n", xlab="", yaxt="n", ylab="", col=colors, ...);
  } else {
    image(z = t(ReverseRows(Matrix)), xaxt = "n", xlab="", yaxt="n", ylab="", col=colors, ...);
  }
  axis(1, labels = FALSE, tick = FALSE)
  axis(2, labels = FALSE, tick = FALSE)
  axis(3, labels = FALSE, tick = FALSE)
  axis(4, labels = FALSE, tick = FALSE)
  nxlabels = length(xLabels)
#   plot x axis labels using:
#   par("usr")[3] - 0.25 as the vertical placement
#   srt = 45 as text rotation angle
#   adj = 1 to place right end of text at tick mark
#   pd = TRUE to allow for text outside the plot region
  plotbox = par("usr");
  xmin = plotbox[1]; xmax = plotbox[2]; ymin = plotbox[3]; yrange = plotbox[4]-ymin;
  ymax = plotbox[4]; xrange = xmax - xmin;

  # print(paste("plotbox:", plotbox[1], plotbox[2], plotbox[3], plotbox[4]));

  nylabels = length(yLabels)
  axis(2, labels = FALSE, tick = FALSE)
  xspacing = 1/(nxlabels-1); yspacing = 1/(nylabels-1)
  # print(paste("nxlabels:", nxlabels));
  if (!xColorLabels)
  {
    if (is.null(cex.lab)) cex.lab = 1;
    text(((1:nxlabels)-1)*xspacing , ymin - 0.02, srt = 45,
          adj = 1, labels = xLabels, xpd = TRUE, cex = cex.lab)
  } else {
    rect(((1:nxlabels)-1)*xspacing - xspacing/2, ymin-xspacing*1.2,
         ((1:nxlabels)-1)*xspacing + xspacing/2, ymin-xspacing*0.2,
         density = -1,  col = substring(xLabels, 3), border = substring(xLabels, 3), xpd = TRUE)
    if (!is.null(xSymbols))
      text ( ((1:nxlabels)-1)*xspacing, ymin-xspacing*1.3, xSymbols, adj = c(0.5, 0), xpd = TRUE);
  }
  if (!yColorLabels)
  {
    if (is.null(cex.lab)) cex.lab = 1;
    text(xmin - 0.01*xrange, ((1:nylabels)-1)/(nylabels-1), srt = 0,
         adj = 1, labels = ReverseVector(yLabels), xpd = TRUE, cex = cex.lab )
  } else {
    rect(xmin-yspacing*1.2, ((nylabels:1)-1)*yspacing - yspacing/2,
         xmin-yspacing*0.2, ((nylabels:1)-1)*yspacing + yspacing/2,
         density = -1,  col = substring(yLabels, 3), border = substring(yLabels, 3), xpd = TRUE)
    if (!is.null(ySymbols))
      text (xmin-yspacing*1.2, ((nylabels:1)-1)*yspacing, ySymbols, adj = c(1, 0.5), xpd = TRUE);
  }

  if (!is.null(NumMatrix))
  {
    if (is.null(cex.Num)) cex.Num = par("cex");
    #if (dim(NumMatrix)!=dim(Matrix))
    #  stop("HeatmapWithTextLabels: NumMatrix was given, but has dimensions incompatible with Matrix.");
    for (rw in 1:dim(Matrix)[1])
      for (cl in 1:dim(Matrix)[2])
      {
        text((cl-1)*xspacing, (dim(Matrix)[1]-rw)*yspacing,
             as.character(NumMatrix[rw,cl]), xpd = TRUE, cex = cex.Num, adj = c(0.5, 0.5));
      }
  }
  axis(1, labels = FALSE, tick = FALSE)
  axis(2, labels = FALSE, tick = FALSE)
  axis(3, labels = FALSE, tick = FALSE)
  axis(4, labels = FALSE, tick = FALSE)
}

#--------------------------------------------------------------------------
#
# BarplotWithTextLabels = function ( Matrix, Labels, ... ) {
#
#--------------------------------------------------------------------------
#
# Plots a barplot of the Matrix and writes the Labels underneath such that they are readable.

BarplotWithTextLabels = function ( Matrix, Labels, ColorLabels = FALSE, Colored = TRUE,
                                   SetMargins = TRUE, Errors = NULL, ... )
{
  if (SetMargins) par(mar=c(3,3,2,2)+0.2)

  if (Colored)
  {
     colors = substring(Labels, 3);
  } else {
     colors = rep("grey", times = ifelse(length(dim(Matrix))<2, length(Matrix), dim(Matrix)[[2]]));
  }

  mp = barplot(Matrix, col = colors, xaxt = "n", xlab="", yaxt="n", ylab="", ...)

  if (length(dim(Matrix))==2) {
     means = apply(Matrix, 2, sum);
  } else {
     means = Matrix;
  }

  err.bp(means, 1.96*Errors, two.side = T);

  # axis(1, labels = FALSE)
  nlabels = length(Labels)
  plotbox = par("usr");
  xmin = plotbox[1]; xmax = plotbox[2]; ymin = plotbox[3]; yrange = plotbox[4]-ymin;
  ymax = plotbox[4];
  # print(paste("yrange:", yrange));
  if (nlabels>1)
  {
     spacing = (mp[length(mp)] - mp[1])/(nlabels-1);
  } else {
     spacing = (xmax-xmin);
  }
  yoffset = yrange/30
  xshift = spacing/2;
  xrange = spacing * nlabels;
  if (ColorLabels)
  {
    #rect(xshift + ((1:nlabels)-1)*spacing - spacing/2.1, ymin - spacing/2.1 - spacing/8,
    #     xshift + ((1:nlabels)-1)*spacing + spacing/2.1, ymin - spacing/8,
    #     density = -1,  col = substring(Labels, 3), border = substring(Labels, 3), xpd = TRUE)
    rect(mp - spacing/2.1, ymin - 2*spacing/2.1 * yrange/xrange - yoffset,
         mp + spacing/2.1, ymin - yoffset,
         density = -1,  col = substring(Labels, 3), border = substring(Labels, 3), xpd = TRUE)
  } else {
    text(((1:nlabels)-1)*spacing +spacing/2 , ymin - 0.02*yrange, srt = 45,
          adj = 1, labels = Labels, xpd = TRUE)
  }
  axis(2, labels = T)
}

#--------------------------------------------------------------------------
#
# SizeWindow
#
#--------------------------------------------------------------------------
# if the current device isn't of the required dimensions, close it and open a new one.

SizeWindow = function(width, height)
{
  din = par("din");
  if ( (din[1]!=width) | (din[2]!=height) )
  {
    dev.off();
    X11(width = width, height=height);
  }
}

#======================================================================================================
# GreenToRed.R
#======================================================================================================

GreenBlackRed = function(n)
{
  half = as.integer(n/2);
  red = c(rep(0, times = half), 0, seq(from=0, to=1, length.out = half));
  green = c(seq(from=1, to=0, length.out = half), rep(0, times = half+1));
  blue = rep(0, times = 2*half+1);
  col = rgb(red, green, blue, maxColorValue = 1);
}

GreenWhiteRed = function(n)
{
  half = as.integer(n/2);
  red = c(seq(from=0, to=1, length.out = half), rep(1, times = half+1));
  green = c(rep(1, times = half+1), seq(from=1, to=0, length.out = half));
  blue = c(seq(from=0, to=1, length.out = half), 1, seq(from=1, to=0, length.out = half));
  col = rgb(red, green, blue, maxColorValue = 1);
}

RedWhiteGreen = function(n)
{
  half = as.integer(n/2);
  green = c(seq(from=0, to=1, length.out = half), rep(1, times = half+1));
  red = c(rep(1, times = half+1), seq(from=1, to=0, length.out = half));
  blue = c(seq(from=0, to=1, length.out = half), 1, seq(from=1, to=0, length.out = half));
  col = rgb(red, green, blue, maxColorValue = 1);
}

#======================================================================================================
# AverageExpression.R
#======================================================================================================

# Getting a mean of expression data classified by a factor:
# The colors are stored in a vector colorhdataOne.
# Getting a mean of columns of a dataframe juts requires a strightforward application of mean.
# or apply(frame, 2, mean).
# Then can do by(frame, factor, function)
# but it returns an object of class "by" which seems to be a list with elements named by the factors.

AverageExprMatrix = function(NormExprData, colors) {
  no.genes = dim(NormExprData)[2]
  no.samples = dim(NormExprData)[1]
  colorsf = as.factor(colors)
  AverageExpr = matrix(ncol=nlevels(colorsf), nrow = no.samples)
  ExprDataMtrx = as.matrix(NormExprData)
  for (i in (1:no.samples)) AverageExpr[i,] = tapply(ExprDataMtrx[i,], colorsf, mean)
  AverageExpr
}

# A wrapper that turns the average expression matrix into a data frame
# with appropriate column and row names
AverageExprFrame = function(NormExprData, colors)
{
  tmp = as.data.frame(AverageExprMatrix(NormExprData, colors),
                row.names = row.names(NormExprData) )
  names(tmp) = paste("AE", levels(as.factor(colors)), sep="")
  tmp
}

# Calculate correlation of Matrix1[,i] with Matrix2[,i] and return as a vector
ColumnCor = function ( Mtrx1, Mtrx2 )
{
  if ( var(dim(Mtrx1) - dim(Mtrx2))!=0 ) {
    stop("ColumnCorr: Error: Given matrices have different dimensions.")
  }
  corrs = vector(mode="numeric", length=dim(Mtrx1)[2] )
  for (i in (1:dim(Mtrx1)[2]) ) corrs[i] = cor(Mtrx1[,i], Mtrx2[,i], use="p")
  names(corrs) = paste(names(Mtrx1), names(Mtrx2))
  corrs
}


#======================================================================================================
# NetworkAndModules
#======================================================================================================

MaxSets = 20;

#-------------------------------------------------------------------------------------------
#
# CheckSets
#
#-------------------------------------------------------------------------------------------
# Checks sets for consistency and returns some diagnostics.

CheckSets = function(ExprData)
{
  No.Sets = length(ExprData);
  if (No.Sets<=0) stop("No expression data given!");
  if (No.Sets>MaxSets) stop(paste("Number of datasets is excessive. To use more than", MaxSets, "sets,",
                      "change the line MaxSets =", MaxSets,
                      "to your needed value and run the script again."));

  No.Samples = vector(length = No.Sets);

  No.Genes = dim(ExprData[[1]]$data)[2];
  for (set in 1:No.Sets)
  {
    if (No.Genes!=dim(ExprData[[set]]$data)[2])
       stop(paste("Incompatible number of genes in set 1 and", set));
    No.Samples[set] = dim(ExprData[[set]]$data)[1];
  }

  list(ngenes = No.Genes, nsamples = No.Samples);
}

#-------------------------------------------------------------------------------------------
#
# KeepCommonProbes
#
#-------------------------------------------------------------------------------------------
# Filters out probes that are not common to all datasets, and puts probes into the same order in each
# set. Works by creating dataframes of probe names and their indices and merging them all.

KeepCommonProbes = function(ExprData, OrderBy = 1)
{
  No.Sets = length(ExprData);
  if (No.Sets<=0) stop("No expression data given!");
  if (No.Sets>MaxSets) stop(paste("Number of datasets is excessive. To use more than", MaxSets, "sets,",
                      "change the line MaxSets =", MaxSets,
                      "to your needed value and run the script again."));

  Names = data.frame(Names = names(ExprData[[OrderBy]]$data));

  if (No.Sets>1) for (set in (1:No.Sets))
  {
    SetNames = data.frame(Names = names(ExprData[[set]]$data), index = c(1:dim(ExprData[[set]]$data)[2]));
    Names = merge(Names, SetNames, by.x = "Names", by.y = "Names", all = FALSE, sort = FALSE);
  }

  for (set in 1:No.Sets)
    ExprData[[set]]$data = ExprData[[set]]$data[, Names[, set+1]];

  ExprData;
}

#-------------------------------------------------------------------------------------------
#
# GetConnectivity
#
#-------------------------------------------------------------------------------------------
# This function takes expression data (rows=samples, colummns=genes),
# and the SoftPower exponent used in weighting the
# correlations to get the network adjacency matrix, and returns an array of dimensions
# No.Genes * No.Sets containing the connectivities of each gene in each subset.

# Caution: KeepOverlapSign and negative powers are not functional.

GetConnectivity = function(ExprData, SoftPower=6, KeepOverlapSign = FALSE,
                           verbose=1, print.level=0, BatchSize = 1500)
{
  No.Sets = length(ExprData);
  spaces = PrintSpaces(print.level);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;

  if (verbose>1) print.flush(paste(spaces, "GetConnectivity: received", No.Sets, "datasets with No.Genes =",
                           as.character(No.Genes)));
  if (verbose>1) print.flush(paste(spaces, "  Received No.Sets =", as.character(No.Sets)));

  Connectivity = matrix(nrow = No.Genes, ncol = No.Sets);

  for (set in 1:No.Sets)
  {
    if (verbose>1) print.flush(paste(spaces, "  Working on set", set));
    No.Batches = as.integer((No.Genes-1)/BatchSize);
    SetRestrConn = NULL;
    for (batch in 1:(No.Batches+1))
    {
      if (batch<=No.Batches)
      {
         if (verbose>2) print.flush(paste(spaces, "    Working on batch", batch));
         BatchIndex = c(1:BatchSize) + (batch-1)*BatchSize;
      } else {
         BatchIndex = c( (BatchSize*(batch-1)+1):No.Genes)
      }

#      adj_mat = AdjacencyMatrix(ExprData[[set]]$data, ExprData[[set]]$data[, BatchIndex], SoftPower = SoftPower,
#                                KeepSign = KeepOverlapSign,
#                                verbose = verbose-1, print.level = print.level+1);
      adj_mat = AdjacencyMatrixR(ExprData[[set]]$data, ExprData[[set]]$data[, BatchIndex],
                                 SoftPower = SoftPower,
                                 verbose = verbose-1, print.level = print.level+1);
      adj_mat[is.na(adj_mat)] = 0;
      Connectivity[BatchIndex, set] = apply(adj_mat, 2, sum)-1;
    }
  }
  Connectivity;
}

#-------------------------------------------------------------------------------
#
# SelectGenesByConnectivity
#
#-------------------------------------------------------------------------------
#
# Here network genes are selected based on 1. nonzero variance in all subsets, 2.
# ranking in a given subset. The retrun value is a boolean vector of length No.Genes
# signifying whether the corresponding gene is to be included in the network.

SelectGenesByConnectivity = function(ExprData, Connectivity, Subs.Ind=1,
                                     DegreeCut=3600, verbose=1, print.level=0)
{
  spaces = PrintSpaces(print.level);

  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;

  variance = matrix(ncol = No.Sets, nrow = No.Genes)

  for (set in (1:No.Sets))
  {
    variance[ ,set] = as.vector(apply(ExprData[[set]]$data, 2, var, na.rm=T))
  }
  VarProduct = as.vector(apply(variance, 1, prod))
  if ( (DegreeCut>0) & (DegreeCut<No.Genes))
  {
      DegreeRank = rank(-Connectivity[, Subs.Ind], ties.method = "first");
      SelectedGenes = DegreeRank <= DegreeCut & VarProduct > 0
  } else {
      SelectedGenes = VarProduct > 0;
  }
  if (verbose>0)
    print.flush(paste(spaces, "SelectGenesByConnectivity:", as.character(sum(VarProduct==0)),
                           "genes have zero variance in at least one subset, selected",
                           sum(SelectedGenes), "probes."));
  SelectedGenes;
}


#-------------------------------------------------------------------------------
#
# SelectGenesByMinConnectivity
#
#-------------------------------------------------------------------------------
#
# Here network genes are selected based on 1. nonzero variance in all subsets, 2.
# ranking of the minimum (over all subsets) connectivity.

SelectGenesByMinConnectivity = function(ExprData, Connectivity, Subs.Ind=1,
                                     DegreeCut=3600, verbose=1, print.level=0)
{
  spaces = PrintSpaces(print.level);

  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;

  variance = matrix(ncol = No.Sets, nrow = No.Genes)

  for (set in 1:No.Sets)
  {
    variance[ ,set] = as.vector(apply(ExprData[[set]]$data, 2, var, na.rm=T))
  }
  VarProduct = as.vector(apply(variance, 1, prod));
  if ( (DegreeCut>0) & (DegreeCut<No.Genes))
  {
     MinConn = as.vector(apply(Connectivity, 1, min));
     ConnRank = rank(-MinConn, ties.method = "first");
     SelectedGenes = ConnRank <= DegreeCut & VarProduct > 0;
  } else {
     SelectedGenes = VarProduct > 0;
  }

  if (verbose>0)
    print.flush(paste(spaces, "SelectGenesByMinConnectivity:", as.character(sum(VarProduct==0)),
                           "genes have zero variance in at least one subset, selected",
                           sum(SelectedGenes), "probes."));
  SelectedGenes;
}

#-------------------------------------------------------------------------------
#
# SelectGenesByRestrConnectivity
#
#-------------------------------------------------------------------------------
#
# Here network genes are selected based on 1. nonzero variance in all subsets, 2.
# calculating the "N best" gene-gene connectivities, 3. ranking their sum
# Assume NAs have been imputed in the data.

SelectGenesByRestrConnectivity = function(ExprData, NBest = 20, SoftPower = 6,
                                     DegreeCut=3600, BatchSize = 1500, verbose=1, print.level=0,
                                     ReturnDiags = FALSE)
{
  spaces = PrintSpaces(print.level);

  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;

  if (verbose>1)
     print.flush(paste(spaces, "SelectGenesByRestrConnectivity: received a dataset with No.Genes =",
                           as.character(No.Genes), " and No.Samples =", as.character(sum(No.Samples))));

  if (verbose>1) print.flush(paste(spaces, "  Received No.Sets =", as.character(No.Sets)));

  RestrConnectivity = matrix(nrow = No.Genes, ncol = No.Sets);

  for (set in 1:No.Sets)
  {
    if (verbose>1) print.flush(paste(spaces, "  Working on set", set));
    SetExprData = ExprData[[set]]$data;
    No.Batches = as.integer(No.Genes/BatchSize);
    SetRestrConn = NULL;
    for (batch in 1:(No.Batches+1))
    {
      if (verbose>2) print.flush(paste(spaces, "    Working on batch", batch));
      if (batch<=No.Batches)
      {
         BatchIndex = c(1:BatchSize) + (batch-1)*BatchSize;
      } else {
         BatchIndex = c( (BatchSize*(batch-1)+1):No.Genes)
      }
      if (verbose>2) print.flush(paste(spaces, "      Calculating adjacencies..."));
      BatchExprData = SetExprData[, BatchIndex];
      # In the following BatchConn[i,j] = cor(BatchExprData[,i], ExprData[,j])
      BatchConn = AdjacencyMatrixR(SetExprData, BatchExprData, SoftPower,
                                  verbose = verbose - 1, print.level = print.level +1);
      BatchConn[is.na(BatchConn)] = 0;
      if (verbose>2) print.flush(paste(spaces, "      Sorting adjacencies..."));
      SortedBatchConn = apply(BatchConn, 2, sort, decreasing = TRUE);
      if (is.null(SetRestrConn))
      {
         SetRestrConn = SortedBatchConn[c(1:NBest), ];
      } else {
         SetRestrConn = cbind(SetRestrConn, SortedBatchConn[c(1:NBest),]);
      }
    }
    if (verbose>2) print.flush(paste(spaces, "    Summing the highest", NBest, "adjacencies..."));
    RestrConnectivity[, set] = apply(SetRestrConn, 2, sum);
    collect_garbage();
  }


  # The variance checking should strictly spekaing not be necessary, since zero variance would lead to a
  # correlation of NA and hence a zero minimum connectivity... but keep it in to remove such probes in a
  # case when the number of non-zero-variance probes is less than DegreeCut.

  if (verbose>2) print.flush(paste(spaces, " Checking non-zero variance..."));
  variance = matrix(ncol = No.Sets, nrow = No.Genes)

  for (set in 1:No.Sets) {
    variance[ ,set] = as.vector(apply(ExprData[[set]]$data, 2, var, na.rm=T))
  }

  VarProduct = as.vector(apply(variance, 1, prod));

  # This removes zero-variance probes from the ranking (again).

  RestrConnectivity[VarProduct==0, ] = 0;

  # Take the minimum of the restricted sums of connectivities...

  MinConn = as.vector(apply(RestrConnectivity, 1, min));
  ConnRank = rank(-MinConn, ties.method = "first");

  # Select genes based on the rank of the minima.

  if (verbose>2) print.flush(paste(spaces, " Selecting genes..."));
  SelectedGenes = ConnRank <= DegreeCut & VarProduct > 0
  if (verbose>0)
    print.flush(paste(spaces, "  SelectGenesByRestrConnectivity:", as.character(sum(VarProduct==0)),
                           "genes have zero variance in at least one subset, selected",
                           sum(SelectedGenes), "probes."));
  collect_garbage();
  if (ReturnDiags)
  {
     RetVal = list(SelectedGenes = SelectedGenes, RestrConnectivity = RestrConnectivity);
  } else {
     RetVal = SelectedGenes;
  }

  RetVal;
}

#---------------------------------------------------------------------
#
# TrafoCorMatrix
#
#---------------------------------------------------------------------
# Transforms correlation matrix to filter noise. Caution, this function discards sign of the correlations
# in the matrix.
#

TrafoCorMatrix = function(Matrix, method, Power, No.Samples)
{
  if (method=="power")
  {
     return (abs(Matrix)^Power);
  } else if (method=="probability")
  {
     abs_M = abs(Matrix);
     norm = pnorm(1, sd = 1/sqrt(No.Samples));
     raw_weight = pnorm(abs_M, sd = 1/sqrt(No.Samples));
     weight = (raw_weight - 0.5)/(norm - 0.5);
     return (abs_M * weight^Power);
  } else stop("Unrecognized \'method\' given.");
}


#---------------------------------------------------------------------
#
# AdjacencyMatrix
#
#---------------------------------------------------------------------
# Computes the adjacency from the expression data: takes cor, transforms it as appropriate and possibly
# adds a sign if requested. No subselection on ExprData is performed.

AdjacencyMatrix = function(ExprData, ExprData2=NULL, SoftPower, KeepSign = FALSE,
                           verbose=1, print.level = 0)
{
  spaces = PrintSpaces(print.level);
  No.Samples = dim(ExprData)[1];

  if (SoftPower<=0)
  {
    method = "probability"; Power = -SoftPower;
  } else {
    method = "power"; Power = SoftPower;
  }
  if (verbose>2) print.flush(paste(spaces, "Transforming the correlation matrix using", method,
                             "method with power", Power));
  cor_mat = cor(ExprData, ExprData2, use="p");
  # trafoed_cor = TrafoCorMatrix(cor_mat, method, SoftPower, No.Samples); collect_garbage();
  if (method=="power")
  {
    trafoed_cor = abs(cor_mat)^Power;
  } else
  {
     abs_M = abs(cor_mat);
     norm = pnorm(1, sd = 1/sqrt(No.Samples));
     raw_weight = pnorm(abs_M, sd = 1/sqrt(No.Samples));
     weight = (raw_weight - 0.5)/(norm - 0.5);
     trafoed_cor = abs_M * weight^Power;
  }
  if (KeepSign)
  {
     sign_cm = sign(cor_mat);
     adj_mat = sign_cm * trafoed_cor;
     #print("Keeping overlap sign");
     #print(paste("No. of negative entries: ", sum(gtom0<0)));
  } else {
     adj_mat = trafoed_cor;
  }
  rm(trafoed_cor, cor_mat);
  collect_garbage();
  adj_mat;
}

# A faster and less memory-intensive version, but it only works for method==power and it won't keep
# overlap sign

AdjacencyMatrixR = function(ExprData, ExprData2=NULL, SoftPower,
                           verbose=1, print.level = 0)
{
  abs(cor(ExprData, ExprData2, use="p"))^SoftPower;
}

#---------------------------------------------------------------------
#
# GetModules
#
#---------------------------------------------------------------------

# Get the vector of module colors for each gene in the given dataset.
# The following code computes the topological overlap matrices for the selected genes
# in each subset:

# CAUTION:  KeepOverlapSign and negative soft powers are not functional in this version.

GetModules = function(ExprData, SelectedGenes, SoftPower = 6,
                      BranchHeightCutoff = 0.94, ModuleMinSize = 125,
                      GetDissTOM = NULL, DissimilarityLevel = 1, KeepOverlapSign = FALSE,
		      ClusterOnlyOne = NULL,
                      verbose = 1, print.level = 0 )
{
  spaces = PrintSpaces(print.level);

  if (verbose>0) print.flush(paste(spaces, "Entering GetModules:"));

  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;

  No.NetworkGenes = sum(SelectedGenes==TRUE)

  dissTOM = vector(mode="list", length = No.Sets)
  hierTOM = vector(mode="list", length = No.Sets);
  colordata = array(dim = c(No.NetworkGenes, No.Sets));

  if (!is.null(GetDissTOM))
  {
    if (!is.null(ClusterOnlyOne))
    {
      if (GetDissTOM!=ClusterOnlyOne)
      {
        stop(paste("GetModules: Error: Incompatible parameters given:",
                   "if both given, GetDissTOM and ClusterOnlyOne must equal. Their values are",
                   GetDissTOM, ClusterOnlyOne));
      }
    } else
    {
      ClusterOnlyOne = GetDissTOM;
    }
  }

  if (is.null(ClusterOnlyOne))
  {
    ClusterOnlyOne_n = 0;
  } else
  {
     ClusterOnlyOne_n = ClusterOnlyOne;
  }

  if (is.null(GetDissTOM))
  {
    GetDissTOM_n = 0;
  } else
  {
     GetDissTOM_n = ClusterOnlyOne;
  }

  if ((!is.null(GetDissTOM)) && (!is.element(GetDissTOM_n, c(1:No.Sets))))
    {
      if (verbose > 0) print.flush(paste(spaces,
                  "  Calculation of DissTOM is disabled by GetDissTOM out of valid range."));
    } else {
        if (verbose>0) print.flush(paste(spaces,
                      "  Calculating TOM and clustering trees..."))
    }
  for (set in 1:No.Sets)
  {
    if (is.null(GetDissTOM) | (GetDissTOM_n==set))
    {
      if (verbose>1) print.flush(paste(spaces,
                    "  Calculating TOM in subset", as.character(set)))
#      gtom0 = AdjacencyMatrix(ExprData[[set]]$data[, SelectedGenes], SoftPower = SoftPower,
#                              KeepSign = KeepOverlapSign,
#                              verbose = verbose, print.level = print.level+1);
      gtom0 = AdjacencyMatrixR(ExprData[[set]]$data[, SelectedGenes], SoftPower = SoftPower,
                              verbose = verbose, print.level = print.level+1);
      if (DissimilarityLevel==0)
      {
         dissTOM[[set]] = list(data = 1-gtom0)
      } else
      {
         dissTOM[[set]] = list(data = SignedTOMdist(gtom0))
      }
      rm(gtom0); collect_garbage();
    }
    collect_garbage()
  }
  if ((!is.null(ClusterOnlyOne)) && (!is.element(ClusterOnlyOne_n, c(1:No.Sets))))
    {
      if (verbose > 0) print.flush(paste(spaces,
                  "  Module detection is disabled by ClusterOnlyOne out of valid range."));
    }
  for (set in 1:No.Sets)
  {
    if (is.null(ClusterOnlyOne) | (ClusterOnlyOne_n==set))
    {
      if (verbose>1) print.flush(paste(spaces,
                    "  Calculating clustering tree in set", set))
      hierTOM[[set]] =
         list(data = hclust(as.dist(dissTOM[[set]]$data),method="average"));
      collect_garbage()
      colordata[,set] = as.character(modulecolor2(hierTOM[[set]]$data,
                      h1=BranchHeightCutoff, minsize1=ModuleMinSize))
      collect_garbage();
    }
  }
  if (!is.null(ClusterOnlyOne) & is.element(ClusterOnlyOne_n, c(1:No.Sets)))
  {
    for (set in 1:No.Sets)
    {
      hierTOM[[set]] = hierTOM[[ClusterOnlyOne]];
      dissTOM[[set]] = dissTOM[[ClusterOnlyOne]];
      colordata[,set] = colordata[,ClusterOnlyOne];
    }
  }

  Modules = list( SelectedGenes = SelectedGenes, Dissimilarity = dissTOM,
                  DissimilarityLevel = DissimilarityLevel,
                  ClusterTree = hierTOM, Colors = colordata,
                  BranchHeightCutoff = BranchHeightCutoff, ModuleMinSize = ModuleMinSize,
                  ClusterOnlyOne = ClusterOnlyOne);
  Modules;
}


#---------------------------------------------------------------------
#
#    GetNetwork
#
#---------------------------------------------------------------------

# ExprData should be a matrix or a data.frame; columns containing valid
# expression data can be specified in ExpressionColumns. By default all columns
# are taken. DegreeCut is the number of genes to be taken into the network
# construction (dissimilarity and clustering). If set to NULL or 0, all genes
# will be taken.
# The return value is a list with the following components:

#network = list(SoftPower = SoftPower, Connectivity = Connectivity, IsInNetwork = SelectedGenes,
#               Dissimilarity = dissTOM, ClusterTree = hierTOM, Colors = colordata);

# This function needs the functions in
# ../CommonFunctions/NetworkFunctions-PL.txt and ../CommonFunctions/PrintFlush.R to be loaded.

GetNetwork= function(ExprData, ProbeSelection = "set connectivity rank",
                      DegreeCut = 3600, NBest = 20, SoftPower = 6,
                      BranchHeightCutoff = 0.94, ModuleMinSize = 125,
                      GetDissTOM = NULL, DissimilarityLevel = 1,
                      KeepOverlapSign = FALSE, ClusterOnlyOne = NULL,
                      verbose = 1, print.level = 0 )
{
  spaces = PrintSpaces(print.level);
  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;

  if (verbose>1) print.flush(paste(spaces, "GetNetwork: received a dataset with No.Genes =",
                           as.character(No.Genes), " and No.Samples =", as.character(sum(No.Samples))));
  if (verbose>1) print.flush(paste(spaces, "  Received No.Sets =", as.character(No.Sets)));

  if (ProbeSelection=="set connectivity rank")
  {
     if ((DegreeCut>0) & (DegreeCut < No.Genes))
     {
         Connectivity = GetConnectivity(ExprData, SoftPower=SoftPower,
                                verbose=verbose-1, print.level = print.level+1);
     } else {
         Connectivity = NULL;
     }
     SelectedGenes = SelectGenesByConnectivity(ExprData, Connectivity, Subs.Ind=1,
                   DegreeCut, verbose, print.level+1);
  } else if (ProbeSelection=="min connectivity rank")
  {
     if ((DegreeCut>0) & (DegreeCut < No.Genes))
     {
         Connectivity = GetConnectivity(ExprData, SoftPower=SoftPower,
                                verbose=verbose-1, print.level = print.level+1);
     } else {
         Connectivity = NULL;
     }
     SelectedGenes = SelectGenesByMinConnectivity(ExprData, Connectivity, Subs.Ind=1,
                   DegreeCut, verbose, print.level+1);
  } else if (ProbeSelection=="min restricted connectivity rank")
  {
     SelectedGenes =  SelectGenesByRestrConnectivity(ExprData, NBest = NBest,
                            SoftPower, DegreeCut, verbose = verbose-1, print.level = print.level+1);
     Connectivity = NULL;
  } else if (ProbeSelection=="all")
  {
     SelectedGenes = rep(TRUE, times = No.Genes);
     Connectivity = NULL;
  } else
      stop(paste("GetNetwork: unrecognized ProbeSelection:", ProbeSelection));

  collect_garbage();collect_garbage();collect_garbage();

  No.NetworkGenes = sum(SelectedGenes==TRUE)

# The following code computes the topological overlap matrices for the selected genes
# in each subset:

  Modules = GetModules(ExprData, SelectedGenes = SelectedGenes, SoftPower = SoftPower,
                       BranchHeightCutoff, ModuleMinSize, GetDissTOM, DissimilarityLevel, KeepOverlapSign,
                       ClusterOnlyOne, verbose, print.level);

  network = list(SoftPower = SoftPower, Connectivity = Connectivity,
                 SelectedGenes = SelectedGenes, Dissimilarity = Modules$Dissimilarity,
                 DissimilarityLevel = Modules$DissimilarityLevel,
                 ClusterTree = Modules$ClusterTree, Colors = Modules$Colors,
                 BranchHeightCutoff = Modules$BranchHeightCutoff,
                 ModuleMinSize = Modules$ModuleMinSize,
                 SignedDiss = KeepOverlapSign);
  network
}

kWithinModule = function(TOM, Colors, gene)
{
  if ((dim(TOM)[[1]]!=dim(TOM)[[2]]) || (length(Colors)!=dim(TOM)[[2]]))
  {
    stop("kWithinModule: Error: TOM is not a square matrix or the dimensions of Colors does not match.");
  }
  color = Colors[gene];
  GeneModuleTOM = TOM[Colors==color, gene];
  sum(GeneModuleTOM);
}

#---------------------------------------------------------------------
#
# PowerConnectivities
#
#---------------------------------------------------------------------
# Computes the connectivities for a given vector of SoftPowers. It is assumed that all entries
# of SoftPowers have the same sign for simplicity.
# No subselection is performed on ExprData

PowerConnectivities = function(SetExprData, SoftPowers, KeepSign = FALSE,
                           verbose=1, print.level = 0, BatchSize = 1536)
{
  spaces = PrintSpaces(print.level);

  No.Genes = dim(SetExprData)[2];
  No.Samples = dim(SetExprData)[1];

  No.Powers = length(SoftPowers);

  Connectivities = matrix(0, nrow = No.Genes, ncol = No.Powers);

  if (SoftPowers[1]<=0)
  {
    method = "probability"; Powers = -SoftPowers;
  } else {
    method = "power"; Powers = SoftPowers;
  }
  if (verbose>2) print.flush(paste(spaces, "Transforming the correlation matrix using", method,
                             "method with powers", paste(Powers, collapse = ", ")));
  No.Batches = as.integer(No.Genes/BatchSize);
  for (batch in 1:(No.Batches+1))
  {
     if (verbose>2) print.flush(paste(spaces, "    Working on batch", batch));
     if (batch<=No.Batches)
     {
        BatchIndex = c(1:BatchSize) + (batch-1)*BatchSize;
     } else {
        BatchIndex = c( (BatchSize*(batch-1)+1):No.Genes)
     }
     if (verbose>2) print.flush(paste(spaces, "      Calculating adjacencies..."));
     BatchExprData = SetExprData[, BatchIndex];
     # In the following BatchConn[i,j] = cor(BatchExprData[,i], ExprData[,j])
     BatchCorr = cor(SetExprData, BatchExprData, use = "p")
     BatchCorr[is.na(BatchCorr)] = 0;

     for (power in 1:No.Powers)
     {
       trafoed_cor = TrafoCorMatrix(BatchCorr, method, Powers[power], No.Samples); collect_garbage();
       if (KeepSign)
       {
          sign_cm = sign(cor_mat);
          adj_mat = sign_cm * trafoed_cor;
       } else {
          adj_mat = trafoed_cor;
       }
       # adj_mat[is.na(adj_mat)] = 0;
       Connectivities[BatchIndex, power] = apply(adj_mat, 2, sum) - 1;
     }
  }
  Connectivities;
}

#---------------------------------------------------------------------
#
# ScaleFreeAnalysis
#
#---------------------------------------------------------------------
# Analyzes the scale-free criterion for a given SetExprData.
# Here SetExprData is assumed to come from only one set and be in a plain matrix format.

ScaleFreeAnalysis = function(SetExprData, Powers, KeepSign = FALSE,
                             No.HistBins = 10,
                             verbose=1, print.level = 0, BatchSize = 1536)
{
  spaces = PrintSpaces(print.level);
  if (verbose>0) print.flush(paste(spaces, "Calculating scale-free diagnostics"));

  # Get connectivities

  if (length(Powers[Powers>0]))
  {
    PositiveConnects = PowerConnectivities(SetExprData, Powers[Powers>0],  KeepSign,
                            verbose = verbose - 1, print.level = print.level + 1);
  } else
     PositiveConnects = NULL;

  if (length(Powers[Powers<=0]))
  {
    NegativeConnects = PowerConnectivities(SetExprData, Powers[Powers<=0],  KeepSign,
                            verbose = verbose - 1, print.level = print.level + 1);
  } else
     NegativeConnects = NULL;


  # Reorganize powers if necessary

  Powers = c(Powers[Powers>0], Powers[Powers<=0]);
  No.Powers = length(Powers);

  Connectivities = cbind(PositiveConnects, NegativeConnects);

  DiagNames = c("Power", "scale law R^2" ,"slope", "truncated R^2","mean(k)","median(k)","max(k)",
                "intercept", "a", "b", "c" );
  Diags = matrix(0, nrow = No.Powers, ncol = length(DiagNames));

  Diags[, 1] = Powers;

  Histos = vector(mode = "list", length = No.Powers);

  # Calculate the fits and diagnostics for each power separately

  for (power in 1:No.Powers)
  {
    minC = min(Connectivities[, power]);
    maxC = max(Connectivities[, power]);
    #breaks = seq(from = 0, to = maxC * (1+1/(2*No.HistBins)) , length.out = No.HistBins +1);
    breaks = seq(from = log10(minC), to = log10(maxC), length.out = No.HistBins +1);
    h = hist(log10(Connectivities[, power]), breaks = breaks, plot = FALSE,
             include.lowest = TRUE, right = TRUE);
    x = h$mids;
    y = log10(h$counts+0.1);
    fit = lm(y ~ x);
    xx = 10^x;
    fit2 = lm(y ~ x + xx);
    Diags[power, 2] = summary(fit)$adj.r.squared
    Diags[power, 3] = summary(fit)$coefficients[2,1]
    Diags[power, 4] = summary(fit2)$adj.r.squared
    Diags[power, 5] = mean(Connectivities[, power]);
    Diags[power, 6] = median(Connectivities[, power]);
    Diags[power, 7] = max(Connectivities[, power]);
    Diags[power, 8] = summary(fit)$coefficients[1,1];
    Diags[power, 9:11] = summary(fit2)$coefficients[1:3,1];
    Histos[[power]] = list(h = h);
  }

  # adjust the format...

  Diags = data.frame(Diags);
  names(Diags) = DiagNames;

  list(Diags = Diags, Histos = Histos);
}


#======================================================================================================
# ModulePrincipalComponents-03.R
#======================================================================================================

#-------------------------------------------------------------------------------------
#
#  ModulePrincipalComponents
#
#-------------------------------------------------------------------------------------

# Calculates the principal components of modules of a given network.
#   - - should multiply them by -1 wherever appropriate.
# Input: Data: expression data, module colors. AlignPCs can take the values "", "along average".
# output : a dataframe of principal components.

ModulePrincipalComponents = function(Data, ModuleColors, AlignPCs = "along average", Impute = FALSE,
   verbose = 1, print.level=0)
{
  spaces = PrintSpaces(print.level);

  AlignPCsRecognizedValues =  c("", "along average");
  if (!is.element(AlignPCs, AlignPCsRecognizedValues)) {
    print.flush(paste("ModulePrincipalComponents: Error:",
                "parameter AlignPCs has an unrecognised value:",
                AlignPCs, "; Recognized values are ", AlignPCsRecognizedValues));
    stop()
  }

  if (verbose>0) print.flush(paste(spaces, "ModulePrincipalComponents: Calculating PCs"));
  FullPCs = ModulePrinComps1(Data, ModuleColors, verbose = verbose-1, print.level = print.level+1,
                             GetConformity = FALSE, Impute = Impute);
  PCs = FullPCs$PrinComps;

  if (AlignPCs == "") AlignedPCs = PCs;
  if (AlignPCs == "along average")
  {
    if (verbose>0) print.flush(paste(spaces,"ModulePrincipalComponents:",
                     "Aligning PCs with average expression for each module."))
    if (verbose>1) print.flush(paste(spaces,"  ++ Calculating averages..."));
    NormData = scale(Data);
    AverageModuleExpr = data.frame(AverageExprMatrix(NormData, ModuleColors));
    if (verbose>1) print.flush(paste(spaces,"  ++ Aligning principal components..."));
    AverageAndPCCor = diag(cor(PCs, AverageModuleExpr, use = "pairwise.complete.obs"));
    sign.matrix = matrix(0, nrow = dim(PCs)[2], ncol = dim(PCs)[2]);
    diag(sign.matrix) = sign(AverageAndPCCor);
    AlignedPCs = as.data.frame(as.matrix(PCs) %*% sign.matrix);
    names(AlignedPCs) = names(PCs);
    rownames(AlignedPCs) = rownames(PCs);
    names(AverageModuleExpr) = names(PCs);
    rownames(AverageModuleExpr) = rownames(PCs);
    if (verbose>1) print.flush(paste(spaces,"  ++ done."));
  }
  RetPCs = list(data = AlignedPCs, VarExplained = FullPCs$varexplained,
                ModuleConformity = FullPCs$ModuleConformity, AverageExpr = AverageModuleExpr);
  RetPCs;
}

#--------------------------------------------------------------------------------------
#
# ModulePCs
# NetworkModulePCs
#
#--------------------------------------------------------------------------------------

ModulePCs = function(ExprData, SelectedGenes, ModuleColors, UniversalModuleColors = NULL,
                     OnlySet = NULL,
                     AlignPCs="along average", Impute = FALSE,  verbose=1, print.level=0)
{
  spaces = PrintSpaces(print.level)
  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;

  if (verbose>0) print.flush(paste(spaces,"ModulePCs: Looking for module PCs."));
  PCs = vector(mode="list", length=No.Sets);
  if (is.null(OnlySet))
  {
      CalculatedSubsets = c(1:No.Sets);
  } else {
      CalculatedSubsets = c(OnlySet);
  }
  for (set in CalculatedSubsets) {
    if (verbose>0) print.flush(paste(spaces,"  Working on subset", as.character(set), "..."));
    if (is.null(UniversalModuleColors))
      {
        SubsetColors = ModuleColors[,set];
      } else {
      SubsetColors = UniversalModuleColors; }
    SubsetPCs = ModulePrincipalComponents(Data = ExprData[[set]]$data[,SelectedGenes==1],
                          ModuleColors = SubsetColors, AlignPCs = AlignPCs, Impute = Impute,
                          verbose = verbose-1, print.level = print.level+1);
    PCs[[set]] = list(data = SubsetPCs$data, AverageExpr = SubsetPCs$AverageExpr,
                           ModuleConformity = SubsetPCs$ModuleConformity,
                           VarExplained = SubsetPCs$VarExplained);
    rm(SubsetColors); rm(SubsetPCs); collect_garbage();
  }
  PCs;
}

# The Network is the same list that is returned by GetModules. It does not contain expression
# data, so those must be given separately.

NetworkModulePCs = function(ExprData, Network, UniversalModuleColors = NULL,
                            OnlySet = NULL,
                            AlignPCs="along average", Impute = FALSE, verbose=1, print.level=0)
{
  PCs = ModulePCs(ExprData, SelectedGenes = Network$SelectedGenes,
            ModuleColors = Network$Colors, UniversalModuleColors, OnlySet = OnlySet,
            AlignPCs, Impute = Impute, verbose, print.level);

  PCs;
}

#--------------------------------------------------------------------------------------
#
# AddTraitToPCs
#
#--------------------------------------------------------------------------------------

# Adds a trait vector to a set of eigenvectors.
# Caution: Traits is assumed to be a vector of lists with each list having an entry data which is
# a No.Samples x No.Traits data frame with an appropriate column name, not a vector.

AddTraitToPCs = function(PCs, Traits, verbose=0, print.level=0)
{
  spaces = PrintSpaces(print.level);

  No.Sets = length(Traits);
  setsize = CheckSets(Traits);
  No.Traits = setsize$ngenes;
  No.Samples = setsize$nsamples;

  if (length(PCs)!=No.Sets)
    stop("Numbers of sets encoded in the length of PCs and Traits parameters differ - must be the same.");

  if (verbose>0) print.flush(paste(spaces, "AddTraitToPCs: Adding traits to principal components."));

  PCTs = vector(mode="list", length=No.Sets);
  for (set in 1:No.Sets)
  {
    trait.subs = Traits[[set]]$data;
    PCT = cbind(PCs[[set]]$data, trait.subs);
    AET = cbind(PCs[[set]]$AverageExpr, trait.subs);
    names(PCT) = c(names(PCs[[set]]$data), names(trait.subs));
    names(AET) = c(names(PCs[[set]]$AverageExpr), names(trait.subs));
    PCTs[[set]] = list(data=PCT, AverageExpr = AET,
                           ModuleConformity = PCs[[set]]$ModuleConformity,
                           VarExplained = PCs[[set]]$VarExplained);
  }
  PCTs;
}


#--------------------------------------------------------------------------------------
#
# OrderPCs
#
#--------------------------------------------------------------------------------------
#
# performs hierarchical clustering on PCs and returns the order suitable for plotting.

OrderNames = function(Names, Order)
{
  if (length(Names)!=length(Order))
  {
    cat("OrderNames: Error: Length of names is different from the length of the order vector.");
    stop();
  }
  OrderedNames = Names;
  for (i in 1:length(Order))
  {
    OrderedNames[i] = Names[Order[i]];
  }
  OrderedNames;
}

OrderPCs = function(PCs, GreyLast = TRUE, GreyName = "PCgrey", OrderBy = 1, Order = NULL, OnlySet = NULL)
{
  if (!is.null(OnlySet)) OrderBy = OnlySet;

  if (is.null(Order))
  {
    print.flush(paste("OrderPCs: order not given, clustering given data in subset", OrderBy));
    corPC = cor(PCs[[OrderBy]]$data, use="p")
    disPC = as.dist(1-corPC);
    clust = hclust(disPC, method = "average");
    Order = clust$order;
  }

  if (length(Order)!=dim(PCs[[OrderBy]]$data)[2])
    stop("OrderPCs: given PCs and Order have incompatible dimensions.");

  if (GreyLast)
  {
    print.flush("OrderPCs:: Putting grey module last");
    ind.grey = 0;
    PCNames = names(PCs[[OrderBy]]$data);
    for (i in 1:length(PCNames))
    {
      if (PCNames[Order[i]]==GreyName) {order.grey = i; ind.grey = Order[i]; }
    }
    if (ind.grey==0)
    {
      print(paste("OrderPCs:: Error: The grey ME name", GreyName,
                  "was not found among the names of the given MEs:", PCNames));
      stop();
    }

    if (order.grey<length(Order))
    {
      for (i in order.grey:(length(Order)-1))
      {
        Order[i] = Order[i+1];
      }
    }
    Order[length(Order)] = ind.grey;
  }
  No.Sets = length(PCs);
  OrderedPCs = vector(mode="list", length = No.Sets);
  if (is.null(OnlySet))
  {
      CalculatedSubsets = c(1:No.Sets);
  } else {
      CalculatedSubsets = c(OnlySet);
  }
  for (subset in CalculatedSubsets)
  {
    OrderedData = PCs[[subset]]$data;
    OrderedAE = PCs[[subset]]$AverageExpr;
    for (col in (1:dim(OrderedData)[2]))
    {
      OrderedData[col] = PCs[[subset]]$data[Order[col]];
      OrderedAE[col] = PCs[[subset]]$AverageExpr[Order[col]];
    }
    names(OrderedData) = OrderNames(names(PCs[[subset]]$data), Order);
    names(OrderedAE) = names(OrderedData);
    OrderedPCs[[subset]] = list(data = OrderedData, VarExplained = PCs[[subset]]$VarExplained,
                                ModuleConformity = PCs[[subset]]$ModuleConformity, AverageExpr = OrderedAE,
                                order = Order);
  }
  OrderedPCs;
}

#---------------------------------------------------------------------------------------------
#
# ConsensusOrderPCs
#
#---------------------------------------------------------------------------------------------
# Orders PCs by the dendrogram of their consensus dissimilarity.

ConsensusOrderPCs = function(PCs, UseAbs = FALSE, OnlySet = NULL, GreyLast = TRUE, GreyName = "MEgrey")
{
  Diss = ConsensusPCDissimilarity(PCs, UseAbs = UseAbs, OnlySet = OnlySet);
  h = hclust(as.dist(Diss));

  OrderPCs(PCs, GreyLast = GreyLast, GreyName = GreyName, Order = h$order);
}

#--------------------------------------------------------------------------------------
#
# CorrelationPreservation
#
#--------------------------------------------------------------------------------------
#
# Given a set of PCs (or OrderedPCs), calculate the preservation values for each module in each pair
# of datasets and return them as a matrix

CorrelationPreservation = function(PCs, Set.Labels, LeaveGreyOut = TRUE)
{
  No.Sets = length(PCs);
  if (No.Sets!=length(Set.Labels)) stop("The lengths of PCs and Set.Labels must equal.");
  if (No.Sets<=1) stop("Something is wrong with argument PCs: its length is 0 or 1");
  Names = names(PCs[[1]]$data);
  if (LeaveGreyOut)
  {
      Use = substring(Names, 3)!="grey";
  } else {
      Use = rep(TRUE, times = length(Names));
  }
  No.Mods = ncol(PCs[[1]]$data[, Use]);
  CP = matrix(0, nrow = No.Mods, ncol = No.Sets*(No.Sets-1)/2);
  diag(CP) = 1;
  CPInd = 1;
  CPNames = NULL;
  for (i in 1:(No.Sets-1))
    for (j in (i+1):No.Sets)
    {
      corPC1 = cor(PCs[[i]]$data[, Use], use="p");
      corPC2 = cor(PCs[[j]]$data[, Use], use="p");
      d = 1-abs(tanh((corPC1 - corPC2) / (abs(corPC1) + abs(corPC2))^2));
      CP[ ,CPInd] = apply(d, 1, sum)-1;
      CPNames = c(CPNames, paste(Set.Labels[i], "::", Set.Labels[j], collapse = ""));
      CPInd = CPInd + 1;
    }
  CPx = as.data.frame(CP);
  names(CPx) = CPNames;
  rownames(CPx) = Names[Use];
  CPx;
}


#--------------------------------------------------------------------------------------
#
# CorrelationPreservation
#
#--------------------------------------------------------------------------------------
#
# Given a set of PCs (or OrderedPCs), calculate the preservation values for each each pair
# of datasets and return them as a matrix.

SetCorrelationPreservation = function(PCs, Set.Labels, LeaveGreyOut = TRUE, method = "absolute")
{
  m = charmatch(method, c("absolute", "hyperbolic"));
  if (is.na(m))
  {
    stop("Unrecognized method given. Recognized methods are absolute, hyperbolic. ");
  }
  No.Sets = length(PCs);
  if (No.Sets!=length(Set.Labels)) stop("The lengths of PCs and Set.Labels must equal.");
  if (No.Sets<=1) stop("Something is wrong with argument PCs: its length is 0 or 1");
  Names = names(PCs[[1]]$data);
  if (LeaveGreyOut)
  {
      Use = substring(Names, 3)!="grey";
  } else {
      Use = rep(TRUE, times = length(Names));
  }
  No.Mods = ncol(PCs[[1]]$data[, Use]);
  SCP = matrix(0, nrow = No.Sets, ncol = No.Sets);
  diag(SCP) = 0;
  for (i in 1:(No.Sets-1))
    for (j in (i+1):No.Sets)
    {
      corPC1 = cor(PCs[[i]]$data[, Use], use="p");
      corPC2 = cor(PCs[[j]]$data[, Use], use="p");
      if (m==1) {
        d = 1 - abs(corPC1 - corPC2)/2;
      } else {
        d = 1-abs(tanh((corPC1 - corPC2) / (abs(corPC1) + abs(corPC2))^2));
      }
      SCP[i,j] = sum(d[upper.tri(d)])/sum(upper.tri(d));
      SCP[j,i] = SCP[i,j];
    }
  SCPx = as.data.frame(SCP);
  names(SCPx) = Set.Labels;
  rownames(SCPx) = Set.Labels;
  SCPx;
}

#--------------------------------------------------------------------------------------
#
# PlotCorPCs
#
#--------------------------------------------------------------------------------------
# Plots a matrix plot of the PC(T)s. On the diagonal the heatmaps show correlation of PCs in the
# particular subset; off-diagonal are differences in the correlation matrix.
# Titles is a vector of titles for the diagonal diagrams; the off-diagonal will have no title
# for now.

# Now using the d = tanh( (C1 - C2)/(|C1| + |C2|)^2 ) measure.

#PlotCorPCs = function(PCs, Titles, Powers = c(1,2))
PlotCorPCs = function(PCs, Titles, ColorLabels = FALSE, colors = NULL, IncludeSign = FALSE,
                      ColoredBarPlot = TRUE, LetterSubPlots = TRUE, Letters = NULL, IncludeGrey = TRUE,
                      setMargins = TRUE, plotCPMeasure = TRUE, plotMeans = FALSE, CPzlim = c(0,1),
                      printCPVals = FALSE, CPcex = 0.9, PlotDiagAdj = FALSE,
                      ...)
{
  #letters = "abcdefghijklmnopqrstuvwxyz";
  if (is.null(Letters)) Letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

  print(PlotDiagAdj);

  if (is.null(colors))
    if (IncludeSign)
    {
      colors = RedWhiteGreen(50);
    } else {
      colors = heat.colors(30);
    }
  No.Sets = length(PCs);
  cex = par("cex");
  mar = par("mar");
  par(mfrow = c(No.Sets, No.Sets));
  par(cex = cex);
  if (!IncludeGrey)
  {
    for (set in 1:No.Sets)
      PCs[[set]]$data = PCs[[set]]$data[ , substring(names(PCs[[set]]$data),3)!="grey"]
  }
  for (i.row in (1:No.Sets))
  {
    for (i.col in (1:No.Sets))
    {
      letter.ind = (i.row-1) * No.Sets + i.col;
      if (LetterSubPlots)
      {
         #letter = paste("(", substring(letters, first = letter.ind, last = letter.ind), ")", sep = "");
         letter = paste( substring(Letters, first = letter.ind, last = letter.ind), ".  ", sep = "");
      } else {
         letter = NULL;
      }
      par(cex = cex);
      if (setMargins) {
        if (ColorLabels) {
          par(mar = c(1,2,3,4)+0.2);
        } else {
          par(mar = c(6,7,3,5)+0.2);
        }
      } else {
        par(mar = mar);
      }
      No.Modules = dim(PCs[[i.col]]$data)[2]
      if (i.row==i.col)
      {
        corPC = cor(PCs[[i.col]]$data, use="p")
        if (IncludeSign)
        {
           if (PlotDiagAdj) {
              HeatmapWithTextLabels((1+corPC)/2, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                 main=paste(letter, Titles[[i.col]]), InvertColors=TRUE, zlim=c(0,1.0),
                                 ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
           } else {
              HeatmapWithTextLabels(corPC, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                    main=paste(letter, Titles[[i.col]]), InvertColors=TRUE, zlim=c(-1,1.0),
                                    ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
           }
        } else {
           HeatmapWithTextLabels(abs(corPC), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                 main=paste(letter, Titles[[i.col]]), InvertColors=TRUE, zlim=c(0,1.0),
                                 ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
        }
      } else
      {
        corPC1 = cor(PCs[[i.col]]$data, use="p");
        corPC2 = cor(PCs[[i.row]]$data, use="p");
        cor.dif = (corPC1 - corPC2)/2;
        d = tanh((corPC1 - corPC2) / (abs(corPC1) + abs(corPC2))^2);
        # d = abs(corPC1 - corPC2) / (abs(corPC1) + abs(corPC2));
        dispd = cor.dif;
        if (plotCPMeasure) dispd[upper.tri(d)] = d[upper.tri(d)];
        if (i.row>i.col)
        {
          if (IncludeSign)
          {
            half = as.integer(length(colors)/2);
            halfColors = colors[1:half+1];
          } else {
            halfColors = colors;
          }
          if (printCPVals) {
            printMtx = matrix(paste(".", as.integer((1-abs(dispd))*100), sep = ""),
                               nrow = nrow(dispd), ncol = ncol(dispd));
            printMtx[printMtx==".100"] = "1";
          } else {
            printMtx = NULL;
          }
          if (sum( (1-abs(dispd)<CPzlim[1]) | (1-abs(dispd)>CPzlim[2]) )>0)
            warning("PlotCorPCs: Correlation preservation data out of zlim range!");
          if (plotCPMeasure) {
             HeatmapWithTextLabels(1-abs(dispd), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                   main=paste(letter, "UT: Cor.Pres\nLT: 1-Cor.Diff"), InvertColors=TRUE,
                                   ColorLabels = ColorLabels, zlim = CPzlim, colors = halfColors,
                                   SetMargins = FALSE,
                                   NumMatrix = printMtx, cex.Num = CPcex, ...);
          } else {
             HeatmapWithTextLabels(1-abs(dispd), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                   main=paste(letter, "Preservation"), InvertColors=TRUE,
                                   ColorLabels = ColorLabels, zlim = CPzlim, colors = halfColors,
                                   SetMargins = FALSE,  NumMatrix= printMtx, cex.Num = CPcex, ...);
          }
        } else {
          if (plotCPMeasure) {
             dp = 1-abs(d);
             method = "";
          } else {
             dp = 1-abs(cor.dif);
             method = "Preservation:";
          }
          diag(dp) = 0;
          if (plotMeans) {
            sum_dp = mean(dp[upper.tri(dp)]);
            means = apply(dp, 2, sum)/(ncol(dp)-1);
            BarplotWithTextLabels(means, names(PCs[[i.col]]$data),
                                 main=paste(letter, "D=", signif(sum_dp,2)),
                                 ylim=c(0,1),
                                 ColorLabels = ColorLabels, Colored = ColoredBarPlot,
                                 SetMargins = FALSE, ... )
          } else {
            sum_dp = sum(dp[upper.tri(dp)]);
            BarplotWithTextLabels(dp, names(PCs[[i.col]]$data),
                                 main=paste(letter, method, "sum = ", signif(sum_dp,3)),
                                 ylim=c(0,dim(dp)[[1]]),
                                 ColorLabels = ColorLabels, Colored = ColoredBarPlot,
                                 SetMargins = FALSE, ... )
          }
        }
      }
    }
  }
}


#--------------------------------------------------------------------------------------
#
# PlotCorPCsAndDendros
#
#--------------------------------------------------------------------------------------
# Plots a matrix plot of the PC(T)s. On the diagonal the heatmaps show correlation of PCs in the
# particular subset; off-diagonal are differences in the correlation matrix.
# Titles is a vector of titles for the diagonal diagrams; the off-diagonal will have no title
# for now.

# Now using the d = tanh( (C1 - C2)/(|C1| + |C2|)^2 ) measure.

#PlotCorPCs = function(PCs, Titles, Powers = c(1,2))
PlotCorPCsAndDendros = function(PCs, Titles, ColorLabels = FALSE, colors = NULL, IncludeSign = FALSE,
                      ColoredBarPlot = TRUE, LetterSubPlots = TRUE, Letters = NULL, IncludeGrey = TRUE,
                      setMargins = TRUE, plotCPMeasure = TRUE, plotMeans = FALSE, CPzlim = c(0,1),
                      printCPVals = FALSE, CPcex = 0.9, plotErrors = FALSE, marDendro = NULL,
                      marHeatmap = NULL, PlotDiagAdj = FALSE,
                      ...)
{
  #Letters = "abcdefghijklmnopqrstuvwxyz";
  if (is.null(Letters)) Letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

  if (is.null(colors))
    if (IncludeSign)
    {
      colors = RedWhiteGreen(50);
    } else {
      colors = heat.colors(30);
    }
  No.Sets = length(PCs);
  cex = par("cex");
  mar = par("mar");
  par(mfrow = c(No.Sets+1, No.Sets));
  par(cex = cex);
  if (!IncludeGrey)
  {
    for (set in 1:No.Sets)
      PCs[[set]]$data = PCs[[set]]$data[ , substring(names(PCs[[set]]$data),3)!="grey"]
  }
  letter.ind = 1;
  for (set in 1:No.Sets)
  {
    #par(cex = StandardCex/1.4);
    par(mar = marDendro);
    labels = names(PCs[[set]]$data);
    uselabels = labels[substring(labels,3)!="grey"];
    corPC = cor(PCs[[set]]$data[substring(labels,3)!="grey",
                                 substring(labels,3)!="grey"], use="p");
    disPC = as.dist(1-corPC);
    clust = hclust(disPC, method = "average");
    if (LetterSubPlots) {
      main = paste(substring(Letters, letter.ind, letter.ind), ". ", Set.Labels[set], sep="");
    } else {
      main = Set.Labels[set];
    }
    plot(clust, main = main, sub="", xlab="",
         labels = substring(uselabels, 3), ylab="", ylim=c(0,1));
    letter.ind = letter.ind + 1;
  }

  for (i.row in (1:No.Sets))
  {
    for (i.col in (1:No.Sets))
    {
      letter.ind = i.row * No.Sets + i.col;
      if (LetterSubPlots)
      {
         #letter = paste("(", substring(Letters, first = letter.ind, last = letter.ind), ")", sep = "");
         letter = paste( substring(Letters, first = letter.ind, last = letter.ind), ".  ", sep = "");
      } else {
         letter = NULL;
      }
      par(cex = cex);
      if (setMargins) {
        if (ColorLabels) {
          par(mar = c(1,2,3,4)+0.2);
        } else {
          par(mar = c(6,7,3,5)+0.2);
        }
      } else {
        par(mar = marHeatmap);
      }
      No.Modules = dim(PCs[[i.col]]$data)[2]
      if (i.row==i.col)
      {
        corPC = cor(PCs[[i.col]]$data, use="p")
        if (IncludeSign)
        {
          if (PlotDiagAdj) {
           HeatmapWithTextLabels((1+corPC)/2, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                 main=paste(letter, Titles[[i.col]]), InvertColors=TRUE, zlim=c(0,1.0),
                                 ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
          } else {
           HeatmapWithTextLabels(corPC, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                 main=paste(letter, Titles[[i.col]]), InvertColors=TRUE, zlim=c(-1,1.0),
                                 ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
          }
        } else {
           HeatmapWithTextLabels(abs(corPC), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                 main=paste(letter, Titles[[i.col]]), InvertColors=TRUE, zlim=c(0,1.0),
                                 ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
        }
      } else
      {
        corPC1 = cor(PCs[[i.col]]$data, use="p");
        corPC2 = cor(PCs[[i.row]]$data, use="p");
        cor.dif = (corPC1 - corPC2)/2;
        d = tanh((corPC1 - corPC2) / (abs(corPC1) + abs(corPC2))^2);
        # d = abs(corPC1 - corPC2) / (abs(corPC1) + abs(corPC2));
        dispd = cor.dif;
        if (plotCPMeasure) dispd[upper.tri(d)] = d[upper.tri(d)];
        if (i.row>i.col)
        {
          if (IncludeSign)
          {
            half = as.integer(length(colors)/2);
            halfColors = colors[1:half+1];
          } else {
            halfColors = colors;
          }
          if (printCPVals) {
            printMtx = matrix(paste(".", as.integer((1-abs(dispd))*100), sep = ""),
                               nrow = nrow(dispd), ncol = ncol(dispd));
            printMtx[printMtx==".100"] = "1";
          } else {
            printMtx = NULL;
          }
          if (sum( (1-abs(dispd)<CPzlim[1]) | (1-abs(dispd)>CPzlim[2]) )>0)
            warning("PlotCorPCs: Correlation preservation data out of zlim range!");
          if (plotCPMeasure) {
             HeatmapWithTextLabels(1-abs(dispd), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                   main=paste(letter, "UT: Cor.Pres\nLT: 1-Cor.Diff"), InvertColors=TRUE,
                                   ColorLabels = ColorLabels, zlim = CPzlim, colors = halfColors,
                                   SetMargins = FALSE,
                                   NumMatrix = printMtx, cex.Num = CPcex, ...);
          } else {
             HeatmapWithTextLabels(1-abs(dispd), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                   main=paste(letter, "Preservation"), InvertColors=TRUE,
                                   ColorLabels = ColorLabels, zlim = CPzlim, colors = halfColors,
                                   SetMargins = FALSE,  NumMatrix= printMtx, cex.Num = CPcex, ...);
          }
        } else {
          if (plotCPMeasure) {
             dp = 1-abs(d);
             method = "Cor.Pres.:";
          } else {
             dp = 1-abs(cor.dif);
             method = "Preservation:";
          }
          diag(dp) = 0;
          if (plotMeans) {
            sum_dp = mean(dp[upper.tri(dp)]);
            means = apply(dp, 2, sum)/(ncol(dp)-1);
            if (plotErrors) {
               Errors = sqrt( (apply(dp^2, 2, sum)/(ncol(dp)-1) - means^2)/(ncol(dp)-2));
            } else {
               Errors = NULL;
            }
            BarplotWithTextLabels(means, names(PCs[[i.col]]$data),
                                 main=paste(letter, "D=", signif(sum_dp,2)),
                                 ylim=c(0,1),
                                 ColorLabels = ColorLabels, Colored = ColoredBarPlot,
                                 SetMargins = FALSE, Errors = Errors, ... )
          } else {
            sum_dp = sum(dp[upper.tri(dp)]);
            BarplotWithTextLabels(dp, names(PCs[[i.col]]$data),
                                 main=paste(letter, method, "sum = ", signif(sum_dp,3)),
                                 ylim=c(0,dim(dp)[[1]]),
                                 ColorLabels = ColorLabels, Colored = ColoredBarPlot,
                                 SetMargins = FALSE, ... )
          }
        }
      }
    }
  }
}

#---------------------------------------------------------------------------------------------
#
# PlotCorPCs_1
#
#---------------------------------------------------------------------------------------------
# This function plots the correlations and differences/preservation only, no barplot

PlotCorPCs_1 = function(PCs, Titles, ColorLabels = FALSE, colors = NULL, IncludeSign = FALSE,
                      ColoredBarPlot = TRUE, ...)
{
  if (is.null(colors))
    if (IncludeSign)
    {
      colors = GreenWhiteRed(50);
    } else {
      colors = heat.colors(30);
    }
  No.Sets = length(PCs);
  par(mfrow = c(No.Sets, No.Sets));
  for (i.row in (1:No.Sets))
  {
    for (i.col in (1:No.Sets))
    {
      No.Modules = dim(PCs[[i.col]]$data)[2]
      if (i.row==i.col)
      {
        corPC = cor(PCs[[i.col]]$data, use="p")
        if (IncludeSign)
        {
           HeatmapWithTextLabels(corPC, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                 main=Titles[[i.col]], InvertColors=TRUE, zlim=c(-1,1.0),
                                 ColorLabels = ColorLabels, colors = colors, ...);
        } else {
           HeatmapWithTextLabels(abs(corPC), names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                 main=Titles[[i.col]], InvertColors=TRUE, zlim=c(0,1.0),
                                 ColorLabels = ColorLabels, colors = colors, ...);
        }
      } else
      {
        corPC1 = cor(PCs[[i.col]]$data, use="p");
        corPC2 = cor(PCs[[i.row]]$data, use="p");
        cor.dif = corPC1 - corPC2;
        d = tanh((corPC1 - corPC2) / (abs(corPC1) + abs(corPC2))^2);
        # d = abs(corPC1 - corPC2) / (abs(corPC1) + abs(corPC2));
        if (IncludeSign)
        {
          half = as.integer(length(colors)/2);
          halfColors = colors[1:half+1];
        } else {
          halfColors = colors;
        }
        if (i.row>i.col)
        {
             HeatmapWithTextLabels(d, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                   main="Correlation NONPreservation", InvertColors=TRUE,
                                   ColorLabels = ColorLabels, zlim=c(-1,1.0), colors = colors, ...);
        } else {
             HeatmapWithTextLabels(cor.dif, names(PCs[[i.col]]$data), names(PCs[[i.col]]$data),
                                   main="Correlation Difference", InvertColors=TRUE,
                                   ColorLabels = ColorLabels, zlim=c(-1,1.0), colors = colors, ...);
        }
      }
    }
  }
}

#---------------------------------------------------------------------------------------------
#
# kME
#
#---------------------------------------------------------------------------------------------
# This function calculates the ME-based connectivity.

kME = function(SetExprData, PCs, ColorF, gene)
{
  color.ind = as.integer(ColorF[gene]);
  #print(paste(gene, color.ind));
  #if (is.na(color.ind)) stop(paste("color.ind is NA for gene", gene));
  #if (is.null(PCs) | is.null(dim(PCs))) stop(paste("PCs is NULL for gene", gene));
  #
  #if (color.ind>ncol(PCs)) stop(paste("For gene", gene, "have color factor", ColorF[gene],
                                      #"whose index", color.ind, "exceeds the dimension", ncol(PCs),
                                      #"of PCs."));
  kme = cor(SetExprData[, gene], PCs[, color.ind], use="p");
  kme;
}


#---------------------------------------------------------------------------------------------
#
# ConsensusPCDissimilarity
#
#---------------------------------------------------------------------------------------------
# This function calcualtes a consensus dissimilarity (i.e., correlation) among sets of PCs (more generally,
# any sets of vectors).
# CAUTION: when not using absolute value, the minimum similarity will favor the large negative values!

ConsensusPCDissimilarity = function(PCs, UseAbs = FALSE, OnlySet = NULL)
{
  MEDiss = vector(mode="list", length = No.Sets);
  if (is.null(OnlySet))
  {
      CalculatedSubsets = c(1:No.Sets);
  } else {
      CalculatedSubsets = c(OnlySet);
  }
  for (set in CalculatedSubsets)
  {
    if (UseAbs)
    {
        diss = 1-abs(cor(PCs[[set]]$data, use="p"));
    } else
    {
        diss = 1-cor(PCs[[set]]$data, use="p");
    }
    MEDiss[[set]] = list(Diss = diss);
  }

  if (is.null(OnlySet))
  {
    ConsDiss = (MEDiss[[1]]$Diss)
    if (No.Sets>1) for (set in 2:No.Sets)
      ConsDiss = pmax(ConsDiss, MEDiss[[set]]$Diss);
  } else {
    ConsDiss = MEDiss[[OnlySet]]$Diss;
  }

  ConsDiss = as.data.frame(ConsDiss);
  names(ConsDiss) = names(PCs[[1]]$data);
  rownames(ConsDiss) = names(PCs[[1]]$data);

  ConsDiss;
}


#---------------------------------------------------------------------------------------------
#
# MergeCloseModules
#
#---------------------------------------------------------------------------------------------
# This function merges modules whose PCs fall on one branch of a hierarchical clustering tree

MergeCloseModules = function(ExprData, Network, SmallModuleColors, CutHeight, OnlySet = NULL,
                             StandardColors = NULL,
                             OrderedPCs = NULL, UseAbs = FALSE, IncludeGrey = FALSE,
                             Relabel = TRUE, Impute = FALSE,
                             verbose = 1, print.level=0)
{
  spaces = PrintSpaces(print.level);
  if (verbose>0) print.flush(paste(spaces,
            "MergeCloseModules: Merging modules whose distance is less than", CutHeight));

  # If ordered PCs were not given, calculate them

  if (is.null(OrderedPCs))
  {
    PCs = NetworkModulePCs(ExprData, Network, UniversalModuleColors = SmallModuleColors,
                           OnlySet = OnlySet, Impute = Impute,
                           verbose = verbose-1, print.level = print.level+1);
    OrderedPCs = OrderPCs(PCs, GreyLast=TRUE, GreyName = "MEgrey", OnlySet = OnlySet );
  } else if (nlevels(as.factor(SmallModuleColors))!=dim(OrderedPCs[[1]]$data)[2])
  {
    if (verbose>0) print.flush(paste(spaces, "MergeCloseModules: Number of goven module colors",
              "does not match number of given MEs => recalculating the MEs."))
    PCs = NetworkModulePCs(ExprData, Network, UniversalModuleColors = SmallModuleColors,
                           OnlySet = OnlySet, Impute = Impute,
                           verbose = verbose-1, print.level = print.level+1);
    OrderedPCs = OrderPCs(PCs, GreyLast=TRUE, GreyName = "MEgrey", OnlySet = OnlySet);
  }

  # Cluster the found module eigengenes and merge ones that are too close to one another _in both sets_.

  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = setsize$ngenes;
  No.Samples = setsize$nsamples;

  MEDiss = vector(mode="list", length = No.Sets);
  if (is.null(OnlySet))
  {
      CalculatedSubsets = c(1:No.Sets);
  } else {
      CalculatedSubsets = c(OnlySet);
  }
  for (set in CalculatedSubsets)
  {
    if (IncludeGrey)
    {
        IndexRange = c(1:(nlevels(as.factor(SmallModuleColors))));
    } else {
        IndexRange = c(1:(nlevels(as.factor(SmallModuleColors))-1));
    }
    if (UseAbs)
    {
        diss = 1-abs(cor(OrderedPCs[[set]]$data[, IndexRange], use="p"));
    } else
    {
        diss = 1-cor(OrderedPCs[[set]]$data[, IndexRange], use="p");
    }
    MEDiss[[set]] = list(Diss = diss);
  }

  if (is.null(OnlySet))
  {
    ConsDiss = (MEDiss[[1]]$Diss)
    if (No.Sets>1) for (set in 2:No.Sets)
      ConsDiss = pmax(ConsDiss, MEDiss[[set]]$Diss);
  } else {
    ConsDiss = MEDiss[[OnlySet]]$Diss;
  }

  METree = hclust(as.dist(ConsDiss), method = "average");
  METreeBranches = as.factor(ModuleNumber(HierTree = METree, CutHeight = CutHeight, MinSize = 1));

  # Analyze the branches: look for the ones that contain more than one original module

  MEUniqueBranches = levels(METreeBranches);
  MENo.Branches = nlevels(METreeBranches)
  MENumberOnBranch = rep(0, times = MENo.Branches);
  for (branch in 1:MENo.Branches)
  {
    MENumberOnBranch[branch] = sum(METreeBranches==MEUniqueBranches[branch]);
  }

  MergedColors = SmallModuleColors;

  # Merge modules on the same branch

  for (branch in 1:MENo.Branches) if (MENumberOnBranch[branch]>1)
  {
    if (verbose>3) print.flush(paste(spaces, "   Working on branch", branch, "having",
                               MENumberOnBranch[branch], "original modules"));
    ModulesOnThisBranch = names(METreeBranches)[METreeBranches==MEUniqueBranches[branch]];
    ColorsOnThisBranch = substring(ModulesOnThisBranch, 3);
    if (verbose>3) print.flush(paste("    Original colors on this branch:", paste(ColorsOnThisBranch,
                                         collapse=", ")));
    for (color in 2:length(ColorsOnThisBranch))
      MergedColors[MergedColors==ColorsOnThisBranch[color]] = ColorsOnThisBranch[1];
  }

  No.Mods = nlevels(as.factor(MergedColors));
  RawModuleColors = levels(as.factor(MergedColors));

  # print(paste("No. of new modules: ", No.Mods));
  # print(paste("Merged module colors:"));
  # print(table(as.factor(MergedColors)));

  MergedNewColors = MergedColors;
  if (Relabel)
  {
    # Relabel the merged colors to the usual order based on the number of genes in each module
    if (is.null(StandardColors))
    {
      StandardColors = c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
                  "purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60",
                  "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise",
                  "darkgrey", "orange", "darkorange", "white" );

    }
    No.GenesInModule = rep(0, No.Mods);
    for (mod in 1:No.Mods) No.GenesInModule[mod] = sum(MergedColors==RawModuleColors[mod]);

    SortedRawModuleColors = RawModuleColors[order(-No.GenesInModule)]

    # Change the color names to the standard sequence, but leave grey grey (that's why rank in general does
    # not equal color)

    if (verbose>3) print(paste(spaces, "   Changing original colors:"));
    rank = 0;
    for (color in 1:length(SortedRawModuleColors)) if (SortedRawModuleColors[color]!="grey")
    {
      rank = rank + 1;
      if (verbose>3) print(paste(spaces, "      ", SortedRawModuleColors[color],
                                 "to ", StandardColors[rank]));
      MergedNewColors[MergedColors==SortedRawModuleColors[color]] = StandardColors[rank];
    }
  }

  list(Colors = MergedNewColors, ClustTree = METree, CutHeight = CutHeight);
}


#======================================================================================================
# ConsensusModules.R
#======================================================================================================

#

#-------------------------------------------------------------------------------------
#
# ConsensusModules
#
#-------------------------------------------------------------------------------------

# This function uses a call to a precompiled function written in C.

SetConsensusModules2 = function(Type = "consensus",
                            ExprData, Network, PCs = NULL, ConsBranchHeightCut = 0.1, ConsModMinSize = 20,
                            Impute = FALSE,
                            verbose = 2, print.level = 0)
{
  if (Type != "consensus" & Type != "majority")
    stop("The Type parameter must be either \'consensus\' or \'majority\'.");

  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Samples = setsize$nsamples;

  spaces = PrintSpaces(print.level);
  if (is.null(PCs))
  {
    print.flush(paste(spaces, "ConsensusModules: PCs not given, calculating them."));
    PCs = NetworkModulePCs(ExprData, Network, UniversalModuleColors = NULL, Impute = Impute,
          verbose = verbose-1, print.level = print.level+1);
  }

  ModuleCor = vector(mode="list", length = No.Sets);
  No.Modules = vector(length = No.Sets);
  GreyIndex = vector(length = No.Sets);

  ModuleIndex = vector(mode="list", length = No.Sets);

  for (i in 1:No.Sets)
  {
    No.Modules[i] = dim(PCs[[i]]$data)[[2]];
    ModuleCor[[i]] = list(cor = data.frame(cor(PCs[[i]]$data, use="p")), GreyIndex = NA,
                           No.Modules = No.Modules[i]);
    names(ModuleCor[[i]]$cor) = names(PCs[[i]]$data);
    module_names = names(PCs[[i]]$data);

    # Find and penalize grey module
    for (j in (1:No.Modules[i])) if (substring(module_names[j], 3)=="grey")
    {
      GreyIndex[i] = j;
      ModuleCor[[i]]$GreyIndex = j;
      ModuleCor[[i]]$cor[j,] = -1;
      ModuleCor[[i]]$cor[,j] = -1;
    }

    ColorF = factor(Network$Colors[,i])
    ModuleIndex[[i]] = list(Ind = as.integer(ColorF));
    rm(ColorF);
  }

  No.Genes = sum(Network$SelectedGenes);

  GeneModuleDiss = array(0, dim = c(No.Genes, No.Genes));
  GeneModuleDiss1 = array(0, dim = c(No.Genes, No.Genes));

#void GeneModuleDiss(double * GeneDiss, int * NGenes, double * ModuleCor, int * NModules,
#                    int * ModuleMembership)

  wd = getwd();
  setwd("../ConsensusModules");
  if (Computer=="genetics-Windows")
  {
    gmd.lib = "GeneModuleDiss.dll";
  } else {
    gmd.lib = "GeneModuleDiss.so";
  }

  if (is.loaded(gmd.lib)) dyn.unload(gmd.lib);

  dyn.load(gmd.lib);
  setwd(wd);

  if (verbose)
  {
    print.flush(paste(spaces, "Calculating module-based gene dissimilarity"));
  }

  for (set in 1:No.Sets)
  {
    if (verbose>1) print.flush(paste(spaces, "++ Working on set ", set));
    ModuleMembership = as.integer(ModuleIndex[[set]]$Ind);
    ModuleCors = as.matrix(ModuleCor[[set]]$cor);
    GeneModuleDiss1 = .C("GeneModuleDiss", GeneModuleDiss1, as.integer(No.Genes), ModuleCors,
                         as.integer(No.Modules[set]), ModuleMembership)[[1]];
    dim(GeneModuleDiss1) = c(No.Genes, No.Genes);
    if (Type=="majority")
    {
       GeneModuleDiss = GeneModuleDiss + GeneModuleDiss1;
    } else {
       GeneModuleDiss = My_pmax(GeneModuleDiss, GeneModuleDiss1);
    }
    rm(ModuleMembership); rm(ModuleCors);
    collect_garbage();
  }

  rm(GeneModuleDiss1);
  rm(ModuleCor); rm(ModuleIndex); rm(GreyIndex);
  collect_garbage();

  GeneModuleDiss = GeneModuleDiss + t(GeneModuleDiss);

  # Cluster genes based on the GeneModuleDiss dissimilarity measure

  if (verbose>0)
  {
    print.flush(paste(spaces, "Clustering the gene module-based dissimilarity..."));
  }

  Cluster = hclust(as.dist(GeneModuleDiss), method = "average");

  if (!is.null(ConsBranchHeightCut))
  {
    ConsensusColor = as.character(modulecolor2(Cluster, h1 = ConsBranchHeightCut, minsize1 = ConsModMinSize));
  } else {
    ConsensusColor = NULL;
  }

  list(ClustTree = Cluster, Colors = ConsensusColor, ConsBranchHeightCut = ConsBranchHeightCut,
       ConsModMinSize = ConsModMinSize, Dissimilarity = GeneModuleDiss);
}

# --------------------------------------------------------------------------------------------------------
#
# Same thing one more time, without the C call of GeneModuleDiss.
#
# --------------------------------------------------------------------------------------------------------

SetConsensusModules = function(Type = "consensus",
                            ExprData, Network, PCs = NULL, ConsBranchHeightCut = 0.1, ConsModMinSize = 20,
                            Impute = FALSE, verbose = 2, print.level = 0)
{
  if (Type != "consensus" & Type != "majority")
    stop("The Type parameter must be either \'consensus\' or \'majority\'.");

  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Samples = setsize$nsamples;

  spaces = PrintSpaces(print.level);
  if (is.null(PCs))
  {
    print.flush(paste(spaces, "ConsensusModules: PCs not given, calculating them."));
    PCs = NetworkModulePCs(ExprData, Network, UniversalModuleColors = NULL, Impute = Impute,
          verbose = verbose-1, print.level = print.level+1);
  }

  ModuleCor = vector(mode="list", length = No.Sets);
  No.Modules = vector(length = No.Sets);
  GreyIndex = vector(length = No.Sets);

  ModuleIndex = vector(mode="list", length = No.Sets);

  for (i in 1:No.Sets)
  {
    No.Modules[i] = dim(PCs[[i]]$data)[[2]];
    ModuleCor[[i]] = list(cor = data.frame(cor(PCs[[i]]$data, use="p")), GreyIndex = NA,
                           No.Modules = No.Modules[i]);
    names(ModuleCor[[i]]$cor) = names(PCs[[i]]$data);
    module_names = names(PCs[[i]]$data);

    # Find and penalize grey module
    for (j in (1:No.Modules[i])) if (substring(module_names[j], 3)=="grey")
    {
      GreyIndex[i] = j;
      ModuleCor[[i]]$GreyIndex = j;
      ModuleCor[[i]]$cor[j,] = -1;
      ModuleCor[[i]]$cor[,j] = -1;
    }

    ColorF = factor(Network$Colors[,i])
    ModuleIndex[[i]] = list(Ind = as.integer(ColorF));
    rm(ColorF);
  }

  No.Genes = sum(Network$SelectedGenes);

  GeneModuleDiss = array(0, dim = c(No.Genes, No.Genes));
  GeneModuleDiss1 = array(0, dim = c(No.Genes, No.Genes));

  if (verbose)
  {
    print.flush(paste(spaces, "Calculating module-based gene dissimilarity"));
  }

  Index = c(1:No.Genes)
  for (set in 1:No.Sets)
  {
    if (verbose>1) print.flush(paste(spaces, "++ Working on set ", set));
    ModuleMembership = as.integer(ModuleIndex[[set]]$Ind);
    ModuleCors = as.matrix(ModuleCor[[set]]$cor);
    GeneModuleDiss1 = ModuleCors[ Modulemebership[Index], ModuleMembership[Index] ];
    dim(GeneModuleDiss1) = c(No.Genes, No.Genes);
    if (Type=="majority")
    {
       GeneModuleDiss = GeneModuleDiss + GeneModuleDiss1/No.Sets;
    } else {
       GeneModuleDiss = My_pmax(GeneModuleDiss, GeneModuleDiss1);
    }
    rm(ModuleMembership); rm(ModuleCors);
    collect_garbage();
  }

  rm(GeneModuleDiss1);
  rm(ModuleCor); rm(ModuleIndex); rm(GreyIndex);
  collect_garbage();

  # GeneModuleDiss = GeneModuleDiss + t(GeneModuleDiss);

  # Cluster genes based on the GeneModuleDiss dissimilarity measure

  if (verbose>0)
  {
    print.flush(paste(spaces, "Clustering the gene module-based dissimilarity..."));
  }

  Cluster = hclust(as.dist(GeneModuleDiss), method = "average");

  if (!is.null(ConsBranchHeightCut))
  {
    ConsensusColor = as.character(modulecolor2(Cluster, h1 = ConsBranchHeightCut, minsize1 = ConsModMinSize));
  } else {
    ConsensusColor = NULL;
  }

  list(ClustTree = Cluster, Colors = ConsensusColor, ConsBranchHeightCut = ConsBranchHeightCut,
       ConsModMinSize = ConsModMinSize, Dissimilarity = GeneModuleDiss);
}

ConsensusModules = function(ExprData, Network, PCs = NULL, ConsBranchHeightCut = 0.1, ConsModMinSize = 20,
                            verbose = 2, print.level = 0)
{
  SetConsensusModules("majority", ExprData, Network, PCs, ConsBranchHeightCut, ConsModMinSize,
                            verbose, print.level);
}
#-------------------------------------------------------------------------------------
#
# AverageModules
#
#-------------------------------------------------------------------------------------

AverageModules = function(Network, ConsBranchHeightCut = NULL, ConsModMinSize = 40,
                            verbose = 2, print.level = 0)
{
  No.Sets = length(Network$Dissimilarity);
  No.Genes = ncol(Dissimilarity[[1]]$data);

  spaces = PrintSpaces(print.level);

  AverageDissTOM = matrix(0, nrow = No.Genes, ncol = No.Genes);

  for (set in (1:No.Sets))
    AverageDissTOM = AverageDissTOM + Network$Dissimilarity[[set]]$data/No.Sets;


  # Cluster genes based on the AverageDissTOM dissimilarity measure

  if (verbose>0)
  {
    print.flush(paste(spaces, "Clustering the average gene dissimilarity..."));
  }

  Cluster = hclust(as.dist(AverageDissTOM), method = "average");

  if (!is.null(ConsBranchHeightCut))
  {
    ConsensusColor = as.character(modulecolor2(Cluster, h1 = ConsBranchHeightCut, minsize1 = ConsModMinSize));
  } else {
    ConsensusColor = NULL;
  }

  list(ClustTree = Cluster, Colors = ConsensusColor, ConsBranchHeightCut = ConsBranchHeightCut,
       ConsModMinSize = ConsModMinSize, Dissimilarity = AverageDissTOM);
}


#-------------------------------------------------------------------------------------
#
# My_pmax
#
#-------------------------------------------------------------------------------------
# Re-implementation of R's pmax, because the latter eats up way too much memory.
# Two versions, one that uses an external function written in C; the other one uses only internal R calls,
# but is slower. Which one is chosen depends on the option UseCpmax ___at the time the function is
# loaded___.

if (UseCpmax)
{
   My_pmax = function(a,b)
   {
     if (Computer=="genetics-Windows")
     {
       pmax.lib = "pmax.dll";
     } else {
       pmax.lib = "pmax.so";
     }

     if (!is.loaded(pmax.lib))
     {
       path = getwd();
       setwd("../CommonFunctions");
       dyn.load(pmax.lib)
       setwd(path);
     }
     if (length(a)==length(b))
     {
       result = .C("pmax", as.double(a), as.double(b), as.integer(length(a)), DUP=FALSE);
       dim(result[[1]]) = dim(a);
     } else {
       stop(paste("My_pmax: the length of parameters a and b differ:", length(a), length(b)));
     }
     result[[1]];
   }
} else
{
   My_pmax = function(a,b, batch=200000, verbose = 0)
   {
     dim = dim(a);
     a = as.vector(a); b = as.vector(b);
     if (length(a)==length(b))
     {
       start = 1; stop = min(start+batch-1, length(a));
       while (start <= length(a))
       {
         if (verbose>1) print.flush(paste("  My_pmax:", start, "through", stop, "of", length(a)));
         c = cbind(a[start:stop], b[start:stop])
         a[start:stop] = apply(c, 1, max);
         start = stop+1; stop = min(start+batch-1, length(a));
       }
       dim(a) = dim;
     } else {
       stop(paste("My_pmax: the length of parameters a and b differ:", length(a), length(b)));
     }
     a;
   }
}

#-------------------------------------------------------------------------------------
#
# IntersectModules
#
#-------------------------------------------------------------------------------------
IntersectModules = function(Network, ConsBranchHeightCut = NULL, ConsModMinSize = 40,
                            verbose = 2, print.level = 0)
{
  No.Sets = length(Network$Dissimilarity);
  No.Genes = sum(Network$SelectedGenes);

  spaces = PrintSpaces(print.level);
  if (verbose>0) print.flush(paste(spaces, "IntersectModules: Calculating minimum",
                                           "gene TOM overlap of given sets"));

  spaces = PrintSpaces(print.level);

  IntersectDissTOM = matrix(0, nrow = No.Genes, ncol = No.Genes);

  #print(paste(memory.size(), memory.size(TRUE)));

  #for (set in 1:No.Sets)
  #{
    #print.flush(paste("Checking on set", set,
                  #paste(dim(Network$Dissimilarity[[set]]$data), collapse=",")));
  #}

  #print(paste(memory.size(), memory.size(TRUE)));
  for (set in 1:No.Sets)
  {
    #print.flush(paste("Working on set", set, paste(dim(IntersectDissTOM), collapse=",") ,
                  #paste(dim(Network$Dissimilarity[[set]]$data), collapse="")));
    IntersectDissTOM = My_pmax(IntersectDissTOM, Network$Dissimilarity[[set]]$data);
    collect_garbage();
  }

  #print(paste(memory.size(), memory.size(TRUE)));
  collect_garbage()
  #print(paste(memory.size(), memory.size(TRUE)));
  # Cluster genes based on the IntersectDissTOM dissimilarity measure

  if (verbose>1)
  {
    print.flush(paste(spaces, "IntersectModules: Clustering the maximum gene dissimilarity..."));
  }

  Cluster = hclust(as.dist(IntersectDissTOM), method = "average");

  if (!is.null(ConsBranchHeightCut))
  {
    ConsensusColor = as.character(modulecolor2(Cluster, h1 = ConsBranchHeightCut, minsize1 = ConsModMinSize));
  } else {
    ConsensusColor = NULL;
  }
  collect_garbage()
  #print(paste(memory.size(), memory.size(TRUE)));

  list(ClustTree = Cluster, Colors = ConsensusColor, ConsBranchHeightCut = ConsBranchHeightCut,
       ConsModMinSize = ConsModMinSize, Dissimilarity = IntersectDissTOM);
}

#-----------------------------------------------------------------------------------------
#
# ConsensusDiagnostics
#
#------------------------------------------------------------------------------------------

# a big WARNING: the order of PCs must be the same as the order of the colors when they are converted to
# a factor. This means PCs must NOT be the OrderedPCs!

ConsensusDiagnostics = function(ExprData, Network, ConsColors, PCs = NULL, Impute = FALSE,
                                verbose=2, print.level=0)
{
  No.Sets = length(ExprData);
  setsize = CheckSets(ExprData);
  No.Genes = sum(Network$SelectedGenes)
  No.Samples = setsize$nsamples;


  spaces = PrintSpaces(print.level);
  if (verbose) print.flush(paste(spaces, "Calculating Consensus Module Diagnostics."));
  if (is.null(PCs))
  {
    PCs = NetworkModulePCs(ExprData, Network, UniversalModuleColors = ConsColors, Impute = Impute,
                           verbose = verbose-1, print.level = print.level+1);
  }
  Modules = factor(ConsColors);
  No.Mods = nlevels(Modules);
  No.AssdGenes = sum(ConsColors!="grey");
  if (verbose>1) print.flush(paste(spaces, "  Calculating IM and ME-based connectivities"));
  # k_me contains ME-based connectivities for every gene in each set
  k_me = matrix(0, ncol = No.Sets, nrow = No.Genes);
  # k_im contains intra-modular connectivities for every gene in each set
  k_im = matrix(0, ncol = No.Sets, nrow = No.Genes);
  for (set in 1:No.Sets)
  {
    TOM = 1-Network$Dissimilarity[[set]]$data;
    for (gene in 1:No.Genes)
    {
      k_me[gene, set] = kME(ExprData[[set]]$data[, Network$SelectedGenes], PCs[[set]]$data,
                            Modules, gene);
      k_im[gene, set] = kWithinModule(TOM, ConsColors, gene);
    }
  }
  if (verbose>1) print.flush(paste(spaces, "  Calculating correlations of IM and ME-based connectivities"));
  k_me_cor = vector(mode="list", length = No.Mods);
  k_im_cor = vector(mode="list", length = No.Mods);
  for (mod in 1:No.Mods)
  {
    color = levels(Modules)[mod];
    kme_c = matrix(0, nrow = No.Sets, ncol = No.Sets);
    kim_c = matrix(0, nrow = No.Sets, ncol = No.Sets);
    kme_p = matrix(0, nrow = No.Sets, ncol = No.Sets);
    kim_p = matrix(0, nrow = No.Sets, ncol = No.Sets);
    for (s1 in 1:No.Sets)
      for (s2 in 1:No.Sets)
      {
        kme1 = k_me[ConsColors==color, s1];
        kme2 = k_me[ConsColors==color, s2];
        kme_c[s1,s2] = cor(kme1, kme2, use="pairwise.complete.obs");
        kme_p[s1,s2] = cor.test(kme1, kme2, use="p", method="s")$p.value;
        kim1 = k_im[ConsColors==color, s1];
        kim2 = k_im[ConsColors==color, s2];
        kim_c[s1,s2] = cor(kim1, kim2, use="p");
        kim_p[s1,s2] = cor.test(kim1, kim2, use="p", method="s")$p.value;
      }
    k_me_cor[[mod]] = list(cor = kme_c, p.value = kme_p);
    k_im_cor[[mod]] = list(cor = kim_c, p.value = kim_p);
  }

  list(No.AssdGenes = No.AssdGenes, k_im = k_im, k_me = k_me, k_im_cor = k_im_cor,
       k_me_cor = k_me_cor, No.Mods = No.Mods);

}

#-----------------------------------------------------------------------------------------
#
# ProbePermutationTest
#
#------------------------------------------------------------------------------------------

ProbePermutationTest = function(#ExprData,
                           Network, No.Perms, method, ConsTreeCut, ConsModMinSize = 40,
                           No.Plots = 0,
                           verbose = 2, print.level = 0)
{
  spaces = PrintSpaces(print.level);

  # No.Sets = length(ExprData);
  # setsize = CheckSets(ExprData);
  # No.Genes = setsize$ngenes;
  # No.Samples = setsize$nsamples;

  No.Sets = length(Network$Dissimilarity);

  # Initialize some useful numbers
  if (verbose)
  {
     print.flush(paste(spaces, "PermutationTest: Will generate", No.Perms,
                           "permutations and consensus modules via", method, "method."));
     print.flush(paste(spaces, "  Using cuts:", paste(ConsTreeCut, collapse = ",")));
  }

  No.Probes = sum(Network$SelectedGenes);
  No.Cuts = length(ConsTreeCut);

  plot = 1;

  # Prepare some statistics

  No.Modules = matrix(NA, nrow = No.Cuts, ncol = No.Perms+1);
  ModuleSizes = vector(mode = "list", length = No.Perms+1);

  # Here we go...
  for (iperm in (1:(No.Perms+1)))
  {
    if ((verbose==2 && iperm/100 == as.integer(iperm/100)) || (verbose>2))
      print.flush(paste(spaces, "Working on permutation", iperm));
    # PermExprData = vector(mode="list", length = No.Sets);
    # PermExprData[[1]] = ExprData[[1]];

    PermNetwork = Network;
    # Permute expression data and the corresponding dissimilarity in Network

    for (set in 1:No.Sets)
    {
      if (iperm!=1)
      {
          perm = sample(No.Probes);
      } else {
          perm = c(1:No.Probes);
      }
      # copy original data
      # PermExprData[[set]] = ExprData[[set]];
      # permute selected samples
      # PermExprData[[set]]$data[, Network$SelectedGenes]
         # = PermExprData[[set]]$data[,Network$SelectedGenes][, perm];
      # Permute dissimilarities in the Network
      PermNetwork$Dissimilarity[[set]]$data = PermNetwork$Dissimilarity[[set]]$data[perm, perm];
    }
    collect_garbage();

    # Get the consensus modules

    if (method=="Average")
    {
      Consensus = AverageModules(Network = PermNetwork,
                            verbose = verbose-1, print.level = print.level+1);
    } else if (method=="Intersection")
    {
      Consensus = IntersectModules(Network = PermNetwork,
                            verbose = verbose-1, print.level = print.level+1);
    } else {
      stop(paste("PermutationTest: The method parameter has an unrecognized value:", method));
    }
    if (plot<=No.Plots)
    {
       if (plot==1) plot(Consensus$ClustTree, labels=FALSE, main = "Original consensus", sub = "",
                         xlab = "", ylab = "");
       if (plot>1) plot(Consensus$ClustTree, labels=FALSE, main = paste("Permutated data", plot-1), sub = "",
                         xlab = "", ylab = "");
       plot = plot + 1;
    }

    collect_garbage();

    # Count modules and their sizes

    PermModuleSizes = vector(mode = "list", length = No.Cuts);
    for (cut in 1:No.Cuts)
    {
      ConsensusColor = ModuleNumber(Consensus$ClustTree, CutHeight = ConsTreeCut[cut],
                                    MinSize = ConsModMinSize);
      PermModuleSizes[[cut]] = list(Sizes = table(ConsensusColor));
      No.Modules[cut, iperm] = length(PermModuleSizes[[cut]]$Sizes);
      #print.flush(paste("Table of ConsensusColor for cut", ConsTreeCut[cut]));
      #print.flush(table(ConsensusColor));
    }
    rm(PermNetwork);
    collect_garbage();
    collect_garbage();
    ModuleSizes[[iperm]] = list(ForCut = PermModuleSizes);
    if (verbose>2) print.flush(paste(spaces, "Permutation", iperm, ": No.Modules for all cuts =",
                               paste(No.Modules[, iperm], collapse = ",")));
  }

  list(No.Modules = No.Modules, ModuleSizes = ModuleSizes);

}


#-----------------------------------------------------------------------------------------------
# Trait selection based on independence and significance
# Assumes that, just like with the PCs, the Traits have the same columns in each dataset (though the
# sample sets need not be the same).

SelectTraits = function(Traits, BranchCut = 0.25, SelectOnSignificance = FALSE, PCs = NULL,
                        SignifThres = 0.03, Impute = FALSE, verbose = 1, print.level = 0)
{
  spaces = PrintSpaces(print.level);
  if (verbose>0) print.flush(paste(spaces, "SelectTraits: Selecting from ", dim(Traits[[1]]$data)[2],
                                   "traits."));
  No.Sets = length(Traits);
  TDiss = 1-cor(Traits[[1]]$data, use = "pairwise.complete.obs");
  if (No.Sets>1) for (set in 2:No.Sets)
  {
     TDiss = pmax(TDiss, 1-cor(Traits[[set]]$data, use = "pairwise.complete.obs"));
  }
  h = hclust(as.dist(TDiss), method = "average");
  TMods = ModuleNumber(h, CutHeight = BranchCut, MinSize = 1);
  No.TMods = nlevels(as.factor(TMods));
  SelTraits = vector(mode="list", length = No.Sets);
  for (set in 1:No.Sets)
  {
    TData = Traits[[set]]$data;
    TData[is.na(TData)] = 0;
    TPCs = ModulePrinComps1(TData, as.factor(TMods), Impute = Impute, verbose = 0,
                            GetConformity = TRUE);
    SelTraits[[set]] = list(data = TPCs$PrinComps);
    for (tmod in 1:No.TMods)
    {
      if (sum(TMods==tmod)>1)
      {
        rnk = order(-TPCs$ModuleConformity[TMods==tmod]);
        SelTraits[[set]]$data[, tmod] = (TData[, TMods==tmod])[, rnk[1]];
        names(SelTraits[[set]]$data)[tmod] = (names(TData)[TMods==tmod])[rnk[1]];
      } else {
        SelTraits[[set]]$data[, tmod] = TData[, TMods==tmod];
        names(SelTraits[[set]]$data)[tmod] = names(TData)[TMods==tmod];
      }
    }
  }

  if (verbose>0) print.flush(paste(spaces, "SelectTraits: Clustering led to ", dim(SelTraits[[1]]$data)[2],
                                   "traits."));

  if (SelectOnSignificance)
  {
    # Reduce further: calculate cor.tests for each ME with each trait in each set;
    # keep only traits that have at least one cor.test$p.value below a threshold

    if (is.null(PCs)) stop("PCs must be given when SelectOnSignificance is requested.");

    No.Mods = dim(PCs[[1]]$data)[2];
    if (is.null(No.Mods))
       stop("Given PCs do not appear to have the correct structure (vector of list",
            "with \'data\' component being a matrix whose columns are PC vectors");

    No.Traits = dim(SelTraits[[1]]$data)[2];

    SelectTrait = rep(FALSE, times = No.Traits);

    for (trait in 1:No.Traits)
      for (mod in 1:No.Mods)
      {
        Significant = TRUE;
        for (set in (1:No.Sets))
        {
          ct = cor.test(PCs[[set]]$data[, mod], SelTraits[[set]]$data[, trait]);
          if (ct$p.value>SignifThres) Significant = FALSE;
        }
        if (Significant) SelectTrait[trait] = TRUE;
      }

    for (set in 1:No.Sets)
    {
      SelTraits[[set]]$data = SelTraits[[set]]$data[, SelectTrait];
    }

    # print(paste("No. of selected traits expected by chance:", No.Mods * No.Traits * TraitThres));
  }

  # Re-cluster the significant traits for diagnostic purposes
  if (sum(SelectTrait)>1)
  {
    TDiss = 1-cor(SelTraits[[1]]$data, use = "pairwise.complete.obs");
    if (No.Sets>1) for (set in 2:No.Sets)
    {
       TDiss = pmax(TDiss, 1-cor(SelTraits[[set]]$data, use = "pairwise.complete.obs"));
    }
    newh = hclust(as.dist(TDiss), method = "average");
  } else {
    newh = NULL;
  }

  if (verbose>0)
  {
    print.flush(paste(spaces, "SelectTraits: Selected", sum(SelectTrait), "traits: "));
    print.flush(paste(spaces, paste(names(SelTraits[[1]]$data), collapse = ", ")));
  }


  list(No.SelectedTraits = sum(SelectTrait), Traits = SelTraits, ClusterTree = h, NewClusterTree = newh);
}

#======================================================================================================
# ColorHandler.R
#======================================================================================================

# A set of global variables and functions that should help handling color names for some 400+ modules.
# A vector called GlobalStandardColors is defined that holds color names with first few entries
# being the well-known and -loved colors. The rest is randomly chosen from the color names of R,
# excluding grey colors.

#---------------------------------------------------------------------------------------------------------
#
# GlobalStandardColors
#
#---------------------------------------------------------------------------------------------------------
# This code forms a vector of color names in which the first entries are given by BaseColors and the rest
# is "randomly" chosen from the rest of R color names that do not contain "grey" nor "gray".

BaseColors = c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
                "purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan",
                "grey60", "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen",
                "darkturquoise", "darkgrey",
                "orange", "darkorange", "white", "skyblue", "saddlebrown", "steelblue",
                "paleturquoise", "violet", "darkolivegreen", "darkmagenta", "white" );

RColors = colors()[-grep("grey", colors())];
RColors = RColors[-grep("gray", RColors)];
InBase = match(BaseColors, RColors);
ExtraColors = RColors[-c(InBase[!is.na(InBase)])];
No.Extras = length(ExtraColors);

# Here is the vector of colors that should be used by all functions:

GlobalStandardColors = c(BaseColors, ExtraColors[rank(sin(13*c(1:No.Extras) +sin(13*c(1:No.Extras))) )] );

rm(BaseColors, RColors, ExtraColors, No.Extras);

#---------------------------------------------------------------------------------------------------------
#
# NormalizeLabels
#
#---------------------------------------------------------------------------------------------------------
# "Normalizes" numerical labels such that the largest group is labeled 1, the next largest 2 etc.
# If KeepZero == TRUE, label zero is preserved.

NormalizeLabels = function(Labels, KeepZero = TRUE)
{
  if (KeepZero)
  {
    NonZero = (Labels!=0);
  }
  else
  {
    NonZero = rep(TRUE, length(Labels));
  }
  f = as.numeric(factor(Labels[NonZero]));
  t = table(Labels[NonZero]);
  # print(t)
  r = rank(-as.vector(t));
  norm_labs = rep(0, times = length(Labels));
  norm_labs[NonZero] = r[f];
  norm_labs;
}

#---------------------------------------------------------------------------------------------------------
#
# ColorsFromLabels
#
#---------------------------------------------------------------------------------------------------------
# This function converts integer numerical labels Labels into color names in the order either given by
# StandardColors,
# or (if StandardColors==NULL) by GlobalStandardColors. If GreyIsZero == TRUE, labels 0 will be assigned
# the color grey; otherwise presence of labels below 1 will trigger an error.

ColorsFromLabels = function(Labels, ZeroIsGrey = TRUE, StandardColors = NULL)
{
  if (is.null(StandardColors)) StandardColors = GlobalStandardColors;
  if (ZeroIsGrey) MinLabel = 0 else MinLabel = 1
  if (sum( (Labels>=MinLabel) & (Labels <= length(StandardColors)) )!= length(Labels))
    stop(paste("Input error: something's wrong with Labels. Either they are not a numeric vector,",
               "or some values are below", MinLabel,
               "or some values are above the maximum number of colors", length(StandardColors)));
  Colors = rep("grey", length(Labels));
  Colors[Labels!=0] = StandardColors[Labels[Labels!=0]];
  Colors;
}

#---------------------------------------------------------------------------------------------------------
#
# DisplayColors
#
#---------------------------------------------------------------------------------------------------------
# This function plots a barplot with all colors given. If Colors are not given, GlobalStandardColors are
# used, i.e. if you want to see the GlobalStandardColors, just call this function without parameters.

DisplayColors = function(Colors = NULL)
{
  if (is.null(Colors)) Colors = GlobalStandardColors;
  barplot(rep(1, length(Colors)), col = Colors, border = Colors);
}



#---------------------------------------------------------------------------------------------------------
#
# ImputeExprData
#
#---------------------------------------------------------------------------------------------------------

# Impute values for NAs. This handles all the t() and as.matrix() operations.

ImputeExprData = function(ExprData)
{
  imputed = data.frame(scale(t(impute.knn(t(as.matrix(ExprData))))))
  names(imputed) = names(ExprData);
  rownames(imputed) = rownames(ExprData);
  imputed;
}

#----------------------------------------------------------------------------------------------
#
# EvaluateClusters
#
#----------------------------------------------------------------------------------------------
# Input: distance matrix and cluster labels. 0 in Labels means the particular point is a
# singleton (not assigned to a cluster).
# Output: a list of cluster quality indicators.

EvaluateClusters = function(DistM, Labels, RefLabels = NULL, Sample = FALSE, SampleProp = 1,
                            verbose = 2, print.level = 0)
{

  spaces = PrintSpaces(print.level);

  d = dim(DistM);
  No.Points = length(Labels);

  if (d[1]!=d[2]) stop("Distance matrix DistM must be square.");

  if (d[1]!=No.Points) stop("Dimension of DistM incompatible with number of cluster labels");

  if (verbose>0) print(paste(spaces, "Evaluating clusters..."));

  # Calculate indices relating to the internal and external links

  if (Sample)
  {
    if ( (SampleProp<=0) | (SampleProp>=1) )
      stop(paste("Incorrect SampleProp parameter given:", SampleProp));
    No.Sampled = as.integer(No.Points * SampleProp);
    Sampled = sample(x = No.Points, size = No.Sampled);
    DistM = DistM[Sampled, Sampled];
    Labels = Labels[Sampled];
    RefLabels = RefLabels[Sampled];
    No.Points = No.Sampled;
  }

  InternalLinks = matrix(0, nrow = No.Points, ncol = No.Points);
  for (point in 1:No.Points) if (Labels[point]!=0)
    InternalLinks[, point] = ifelse(Labels==Labels[point], Labels, 0);

  IntLinks = as.vector(InternalLinks[upper.tri(InternalLinks)]);
  Dist = as.vector(DistM[upper.tri(DistM)]);
  No.Internal = sum(IntLinks!=0);
  No.External = length(IntLinks) - No.Internal;

  if ( (No.Internal > 0) & (No.External>0) )
  {
     ord = order(Dist, c(1:length(Dist)));

     # ...number of "misplaced" internal and external links

     s1 = Dist[ord[No.Internal]];
     s2 = Dist[ord[No.Internal+1]];
     if (s1==s2)
     {
       MaxInt = s1 + 1e-10;
       MinExt = s1 - 1e-10;
     } else {
       MaxInt = s1;
       MinExt = s2;
     }
     No.Misplaced = sum( ( (IntLinks>0) & (Dist>MaxInt) ) | ( (IntLinks==0) & (Dist<MinExt) ) );

     # ...weight index of the internal and external links

     DistInt = sum(Dist[IntLinks>0]);
     DistExt = sum(Dist[IntLinks==0]);
     DistSm = sum(Dist[ord <= No.Internal]);
     DistLg = sum(Dist[ord > No.Internal]);
     DistIndex = DistInt/DistSm * DistLg/DistExt;
  } else {
     No.Misplaced = 0;
     DistIndex = 0;
  }

  if (!is.null(RefLabels))
  {
    if (length(RefLabels)!=No.Points)
       stop("Length of given RefLabels incompatible with number of points.");

    RefInternalLinks = matrix(0, nrow = No.Points, ncol = No.Points);
    for (point in 1:No.Points) if (RefLabels[point]!=0)
      RefInternalLinks[, point] = ifelse(RefLabels==RefLabels[point], RefLabels, 0);

    RefIntLinks = as.vector(RefInternalLinks[upper.tri(RefInternalLinks)]);

    No.Agreed = sum( xor(IntLinks>0, RefIntLinks==0));
    No.Disagreed = sum( xor(IntLinks==0, RefIntLinks==0));

    RandIndex = No.Agreed/(No.Agreed+No.Disagreed);

    # Make each unassigned label unique
    XLabels = Labels;
    NUnassigned = sum(XLabels==0);
    start = max(XLabels);
    XLabels[XLabels==0] = start + c(1:NUnassigned);

    # Same for the reference labels
    XRefLabels = RefLabels;
    NRefUnassigned = sum(XRefLabels==0);
    start = max(XRefLabels);
    XRefLabels[XRefLabels==0] = start + c(1:NRefUnassigned);

# To use the adjustedRandIndex, must load the mclust library that requires a license...
#    AdjRandIndex = adjustedRandIndex(XLabels, XRefLabels);  # This requires package mclust
    AdjRandIndex = NA;
  } else {
    RandIndex = 0;
    AdjRandIndex = 0;
  }

  list(No.Misplaced = No.Misplaced, DistIndex = DistIndex, RandIndex = RandIndex,
       AdjRandIndex = AdjRandIndex);
}



setwd(WorkingDirectory);

#----------------------------------------------------------------------------
#
# CausalChildren
#
#----------------------------------------------------------------------------
# Note: The returned vector may contain multiple occurences of the same child.

CausalChildren = function(Parents, Cause)
{
  No.Nodes = dim(Cause)[[1]];

  # print(paste("Length of Parents: ",length(Parents)));
  if (length(Parents)==0) return(NULL);

  Child_ind = apply(as.matrix(abs(Cause[, Parents])), 1, sum)>0;
  if (sum(Child_ind)>0)
  {
     children = c(1:No.Nodes)[Child_ind]
  } else {
     children = NULL;
  }
  children;
}

# Old version:
#CausalChildren = function(Parents, Cause)
#{
  #children = NULL;
#
  #No.Nodes = dim(Cause)[[1]];
#
  #for (i in 1:length(Parents))
  #{
    #for (j in 1:No.Nodes)
    #{
      #if (Cause[j, Parents[i]]!=0)
      #{
        #children = c(children, j);
      #}
    #}
  #}
  #children;
#}


#----------------------------------------------------------------------------
#
# CreateSeedVectors
#
#----------------------------------------------------------------------------
#
# Given a set of causal anchors, this function creates a network of vectors that should satisfy the
# causal relations encoded in the causal matrix Cause, i.e. Cause[j,i] is the causal effect of vector i on
# vector j. (This is Jason's convention.)

# The function starts by initializing all vectors to noise given in the noise specification. (The noise
# can be specified for each vector separately.) Then it runs the standard causal network signal
# propagation and returns the resulting vectors.

CreateSeedVectors = function(Cause, AnchorInd, AnchorVecs, Noise, verbose = 2, print.level = 0)
{
  spaces = PrintSpaces(print.level);

  if (verbose>0) print.flush(paste(spaces, "Creating seed vectors..."));
  No.Nodes = dim(Cause)[[1]];
  No.Samples = dim(AnchorVecs)[[1]];

  if (length(AnchorInd)!=dim(AnchorVecs)[[2]])
  {
    print(paste("CreateSeedVectors: Error: Length of AnchorInd must equal",
            "the number of vectors in AnchorVecs."));
    stop();
  }
  if (length(Noise)!=No.Nodes)
  {
    print(paste("CreateSeedVectors: Error: Length of Noise must equal",
            "the number of nodes as given by the dimension of the Cause matrix."));
    stop();
  }

  # Initialize all node vectors to noise with given standard deviation

  NodeVectors = matrix(0, nrow = No.Samples, ncol = No.Nodes);
  for (i in 1:No.Nodes)
  {
    NodeVectors[,i] = rnorm(n=No.Samples, mean=0, sd=Noise[i]);
  }

  Levels = rep(0, times = No.Nodes);

  # Calculate levels for all nodes: start from anchors and go through each successive level of children

  level = 0;
  Parents = AnchorInd;
  Children = CausalChildren(Parents = Parents, Cause = Cause);
  if (verbose>1) print.flush(paste(spaces, "..Determining level structure..."));
  while (!is.null(Children))
  {
    # print(paste("level:", level));
    # print(paste("   Parents:", Parents));
    # print(paste("   Children:", Children));
    level = level + 1;
    if ((verbose>1) & (level/10 == as.integer(level/10)))
          print.flush(paste(spaces, "  ..Detected level", level));
    #print.flush(paste("Detected level", level));
    Levels[Children] = level;
    Parents = Children;
    Children = CausalChildren(Parents = Parents, Cause = Cause);
  }

  HighestLevel = level;

  # Generate the whole network

  if (verbose>1) print.flush(paste(spaces, "..Calculating network..."));
  NodeVectors[,AnchorInd] = NodeVectors[,AnchorInd] + AnchorVecs;
  for (level in (1:HighestLevel))
  {
    if ( (verbose>1) & (level/10 == as.integer(level/10)) )
      print.flush(paste(spaces, " .Working on level", level));
    #print.flush(paste("Working on level", level));
    LevelChildren = c(1:No.Nodes)[Levels==level]
    for (child in LevelChildren)
    {
      LevelParents = c(1:No.Nodes)[Cause[child, ]!=0]
      for (parent in LevelParents)
        NodeVectors[, child] = scale(NodeVectors[, child] + Cause[child, parent]*NodeVectors[,parent]);
    }
  }

  Nodes = list(Vectors = NodeVectors, Cause = Cause, Levels = Levels, AnchorInd = AnchorInd);
  Nodes;
}

# -04:
#    Adding variable maximum correlation determined as runif between MIN and MAX, with MAX a
#    constant and MIN depending on number of genes N as MIN = (AbsMax-AbsMin) * exp(N/N0) + AbsMin,
#    that is MIN = AbsMin for N=0 and MIN = AbsMax for N->infinity, with characteristic crossover at
#    N0. Note that minimum correlation in a module is scaled down proportionally to max correlation
#    so the correlation bounds make sense.
#
#    SmallLayer modules have similarly varied maximum correlation, but there N0 equals the
#    AveragModuleSize.


#--------------------------------------------------------------------------------------------
#
# RandomMaxCorr
#
#--------------------------------------------------------------------------------------------

RandomMaxCorr = function(NGenes, N0 = 100, AbsMin = 0.5, AbsMax = 1)
{
  runif(n=1, min = (AbsMax - AbsMin) * (1-exp(-NGenes/N0)) + AbsMin, max = AbsMax);
}

#--------------------------------------------------------------------------------------------
#
# SimulateModule
#
#--------------------------------------------------------------------------------------------
# The resulting data is normalized.
# CorPower controls how fast the correlation drops with index i in the module; the curve is roughly
# x^{1/CorPower} with x<1 and x~0 near the "center", so the higher the power, the faster the curve rises.

SimulateModule = function(ME, NModGenes, NNearGenes, MinCorr, MaxCorr = 1, CorPower = 1,
                          verbose = 1, print.level = 0)
{
    #spaces = PrintSpaces(print.level);
    #print.flush(paste(spaces, "SimulateModule: MinCorr = ", MinCorr, "MaxCorr = ", MaxCorr));
    NSamples = length(ME);

    ExprData = matrix(rnorm((NModGenes+NNearGenes)*NSamples), nrow = NSamples,
                            ncol = NModGenes+NNearGenes)

    VarME = var(ME)

    # generate the in-module genes
    CorME = MaxCorr - (c(1:NModGenes)/NModGenes)^(1/CorPower) * (MaxCorr-MinCorr);
    noise = sqrt(VarME * (1-CorME^2)/CorME^2);
    for (gene in 1:NModGenes)
    {
      ExprData[, gene] = ME + rnorm(NSamples, sd = noise[gene]);
    }
    # generate the near-module genes

    CorME = c(1:NNearGenes)/NNearGenes * MinCorr;
    noise = sqrt(VarME * (1-CorME^2)/CorME^2);
    for (gene in 1:NNearGenes)
    {
      ExprData[, gene] = ME + rnorm(NSamples, sd = noise[gene]);
    }

    ExprData = scale(ExprData);

    ExprData;

}

#--------------------------------------------------------------------------------------------
#
# SimulateSmallLayer
#
#--------------------------------------------------------------------------------------------
# Simulates a bunch of small and weakly expressed modules.

SimulateSmallLayer = function(Order, NSamples,
                              MinCorr = 0.5, MinMaxCorr = 0.5, MaxMaxCorr = 1,
                              AverageModuleSize, AverageExpr, InvDensity,
                              verbose = 4, print.level = 0)
{
  spaces = PrintSpaces(print.level);
  NGenes = length(Order)
  ExprData = matrix(0, nrow = NSamples, ncol = NGenes);

  MaxCorrN0 = AverageModuleSize;

  if (verbose>0) print.flush(paste(spaces, "SimulateSmallLayer: simulating modules with min corr",
          MinCorr, ", average expression", AverageExpr, ", average module size", AverageModuleSize,
          ", inverse density", InvDensity));

  index = 0;
  while (index < NGenes)
  {
     ModSize = as.integer(rexp(1, 1/AverageModuleSize));
     if (ModSize<3) ModSize = 3;
     if (index + ModSize>NGenes) ModSize = NGenes - index;
     if (ModSize>2)   # Otherwise don't bother :)
     {
       ModuleExpr = rexp(1, 1/AverageExpr);
       if (verbose>4) print.flush(paste(spaces, "  Module of size", ModSize, ", expression", ModuleExpr,
                                  ", min corr", MinCorr,
                                  "inserted at index", index+1));
       ME = rnorm(NSamples, sd = ModuleExpr);
       NInModule = as.integer(ModSize*2/3);
       NNearModule = ModSize - NInModule;
       MaxCorr = RandomMaxCorr(NInModule, N0 = AverageModuleSize, AbsMin = MinMaxCorr,
                              AbsMax = MaxMaxCorr);
       EffMinCorr = MinCorr * MaxCorr;
       ExprData[, (index+1):(index + ModSize)] =
           ModuleExpr * SimulateModule(ME, NInModule, NNearModule, EffMinCorr, MaxCorr, CorPower);
     }
     index = index + ModSize * InvDensity;
  }
  ExprData;
}


#--------------------------------------------------------------------------------------------
#
# SimulateExprData
#
#--------------------------------------------------------------------------------------------
#
# Caution: the last Mod.Props entry gives the number of "true grey" genes;
# the corresponding MinCorr entry must be absent (i.e. length(MinCor) = length(ModProps)-1

# OrderedSmallLayers: layers of small modules with weaker correlation, ordered in the same order as the
# genes in the big modules. Needs average number of genes in a module (exponential distribution),
# average expression strength (exponential density) and inverse density.

# UnorderedSmallLayers: Layers of small modules whose order is random.

if(exists("SimulateExprData") ) rm(SimulateExprData);

SimulateExprData=function(SeedVectors, NGenes, ModProps,
                          MinCorr = 0.5, MaxCorr = 1, MinMaxCorr = 0.5, MaxMaxCorr = 1, MaxCorrN0 = 100,
                          CorPower = 1, BackgrNoise = 0.1, LeaveOut = NULL,
                          ModuleColors = NULL,
			  NOrderedSmallLayers = 0, NUnorderedSmallLayers = 0,
                          AverageNGenesInSmallModule = 10, AverageExprInSmallModule = 0.2,
                          InvDensityOfSmallModules = 2,
                          verbose = 1, print.level = 0)
{
    spaces = PrintSpaces(print.level);

    NMods=length(ModProps)-1;

    NSamples = dim(SeedVectors)[[1]];

    if (length(MinCorr)==1) MinCorr = rep(MinCorr, NMods);
    if (length(MaxCorr)==1) MaxCorr = rep(MaxCorr, NMods);

    if (length(MinCorr)!=NMods)
      stop(paste("SimulateExprData: Input error: MinCorr is an array of different lentgh than",
                  "the length-1 of ModProps array."));

    if (length(MaxCorr)!=NMods)
      stop(paste("SimulateExprData: Input error: MaxCorr is an array of different lentgh than",
                  "the length-1 of ModProps array."));

    if (dim(SeedVectors)[[2]]!=NMods)
       stop(paste("SimulateExprData: Input error: Number of seed vectors must equal the",
                  "length of ModProps."));

    if(is.null(ModuleColors))
    {
        ModuleColors=c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
                "purple","greenyellow","tan","salmon", "midnightblue", "lightcyan","grey60",
                "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise",
                "darkgrey", "orange", "darkorange", "white" )

    }

    grey = "grey";

    if(sum(ModProps)>=1) stop("SimulateExprData: Input error: the sum of Mod.Props must be less than 1");
    #if(sum(ModProps[c(1:(length(ModProps)-1))])>=0.5)
           #print(paste("SimulateExprData: Input warning: the sum of ModProps for proper modules",
                                         #"should ideally be less than 0.5."));

    no.in.modules = as.integer(NGenes*ModProps);
    no.in.proper.modules = no.in.modules[c(1:(length(ModProps)-1))];
    no.near.modules = as.integer((NGenes - sum(no.in.modules)) *
                           no.in.proper.modules/sum(no.in.proper.modules));

    simulate.module = rep(TRUE, times = NMods);
    if (!is.null(LeaveOut)) simulate.module[LeaveOut] = FALSE;

    no.in.modules[NMods+1] = NGenes - sum(no.in.proper.modules[simulate.module]) -
                                            sum(no.near.modules[simulate.module]);

    SortedModColors = ModuleColors[rank(-no.in.proper.modules, ties.method = "first")];
    SortedModColors = c(SortedModColors, grey);

    if (verbose>0) print.flush(paste(spaces, "SimulateExprData: simulating", NGenes, "genes in",
                          NMods, "modules."));

    if (verbose>1)
    {
    #  print.flush(paste(spaces, "    Minimum correlation in a module is", MinCorr,
    #                            " and its dropoff is characterized by power", CorPower));
      print.flush(paste(spaces, "    Simulated colors:",
                         paste(SortedModColors[1:NMods], collapse = ", "), " and ", grey));
      print.flush(paste(spaces, "    Module sizes:", paste(no.in.modules, collapse = ", ")));
      print.flush(paste(spaces, "    near module sizes:", paste(no.near.modules, collapse = ", ")));
      print.flush(paste(spaces, "    Min correaltion:", paste(MinCorr, collapse = ", ")));
      if (!is.null(LeaveOut)) print.flush(paste(spaces, "    _leaving out_ modules",
                                                paste(SortedModColors[LeaveOut], collapse = ", ")));

    }

    truemodule=rep(grey, NGenes);
    TrueAllColors=rep(grey, NGenes);	# These have the colors for left-out modules as well.

    # This matrix contains the simulated expression values (rows are genes, columns samples)
    # Each simulated cluster has a distinct mean expression across the samples

    ExprData = matrix(rnorm(NGenes*NSamples), nrow = NSamples, ncol = NGenes)

    gene.index = 0;		# Where to put the current gene into ExprData

    for(mod in c(1:NMods))
    {
       NModGenes = no.in.modules[mod];
       NNearGenes = no.near.modules[mod];
       if (simulate.module[mod])
       {
         ME = SeedVectors[, mod];
         EffMaxCorr = MaxCorr[mod] * RandomMaxCorr(NModGenes, N0 = MaxCorrN0, AbsMin = MinMaxCorr,
                                AbsMax = MaxMaxCorr);
         EffMinCorr = MinCorr[mod] * EffMaxCorr;
         if (verbose>3)
           print.flush(paste(spaces, "    ..working on", SortedModColors[mod],
                             "module: EffMaxCor = ", EffMaxCorr, ", EffMinCor =", EffMinCorr));
         ExprData[, (gene.index+1):(gene.index+NModGenes+NNearGenes)] =
            SimulateModule(ME, NModGenes, NNearGenes, EffMinCorr, EffMaxCorr,
                           CorPower, verbose, print.level);
         truemodule[(gene.index+1):(gene.index+NModGenes)] = SortedModColors[mod];
       }
       TrueAllColors[(gene.index+1):(gene.index+NModGenes)] = SortedModColors[mod];
       gene.index = gene.index + NModGenes + NNearGenes;
    }

    if (NOrderedSmallLayers>0)
    {
      OrderVector = c(1:NGenes)
      for (layer in 1:NOrderedSmallLayers)
      {
        if (verbose>1) print.flush(paste(spaces, "Simulating ordereded extra layer", layer));
        ExprData = ExprData + SimulateSmallLayer(OrderVector, NSamples, MinCorr[1],
                                      MinMaxCorr, MaxMaxCorr,
                                      CorPower, AverageNGenesInSmallModule,
                                      AverageExprInSmallModule, InvDensityOfSmallModules,
                                      verbose-1, print.level+1);
        collect_garbage();
      }
    }

    if (NUnorderedSmallLayers>0) for (layer in 1:NUnorderedSmallLayers)
    {
      if (verbose>1) print.flush(paste(spaces, "Simulating unordereded extra layer", layer));
      OrderVector = sample(NGenes)
      ExprData = ExprData + SimulateSmallLayer(OrderVector, NSamples, MinCorr[1],
                                      MinMaxCorr, MaxMaxCorr, CorPower,
                                      AverageNGenesInSmallModule,
                                      AverageExprInSmallModule, InvDensityOfSmallModules,
                                      verbose = verbose-1, print.level = print.level+1);
      collect_garbage();
    }

    if (verbose>1) print.flush(paste(spaces, "  Adding background noise with amplitude", BackgrNoise));

    ExprData = ExprData + rnorm(n = NGenes*NSamples, sd = BackgrNoise);

    list(ExprData = ExprData, TrueColors = truemodule, TrueAllColors = TrueAllColors,
         ColorOrder = SortedModColors)

} # end of function

#--------------------------------------------------------------------------------------
#
# SimulateSets
#
#--------------------------------------------------------------------------------------
# Simulate several sets with some of the modules left out.
# SeedVectors are specified as one column for one module and all sets;
# LeaveOut must be a matrix of No.Modules x No.Sets of TRUE/FALSE values;
# MinCorr must be a single number here; ModProps are a single vector, since the proportions should be the
# same for all sets.
# No.Samples is a vector specifying the number of samples in each set; this must be compatible with the
# dimensions of the seed vectors.

SimulateSets = function(SeedVectors, NGenes, No.Samples, ModProps,
                          MinCorr = 0.5, MaxCorr = 1, MinMaxCorr = 0.5, MaxMaxCorr = 1, MaxCorrN0 = 100,
                          CorPower = 1, BackgrNoise = 0.1, LeaveOut = NULL,
                          ModuleColors = NULL,
			  NOrderedSmallLayers = 0, NUnorderedSmallLayers = 0,
                          AverageNGenesInSmallModule = 10, AverageExprInSmallModule = 0.2,
                          InvDensityOfSmallModules = 2,
                          verbose = 1, print.level = 0)
{
  d = dim(SeedVectors)
  if (length(d)!=2) stop("SimulateTwoSets: Input Error: SeedVectors must have 2 dimensions.");
  NSets = length(No.Samples);

  NMods = d[2];
  NSamples = sum(No.Samples);

  SetStarts = c(1, cumsum(No.Samples)+1);
  SetEnds = cumsum(No.Samples);

  if (d[1] != NSamples) stop(paste("SimulateSets: Input Error: Incompatible number of",
                                              "samples from SeedVectors and No.Samples."));
  d2 = length(ModProps)-1;
  if (d2 != d[2]) stop(paste("SimulateTwoSets: Input Error: Incompatible number of modules determined",
                                "from the dimensions of SeedVectors and ModProps"));
  d3 = dim(LeaveOut);
  if ( (d3[1] != NMods) | (d3[2] != NSets) )
    stop(paste("SimulateTwoSets: Input Error: incompatible",
               "dimensions of LeaveOut with NMods or NSets"));

  ExprData = vector(mode="list", length = NSets);
  TrueColors = NULL;
  TrueAllColors = NULL;
  ColorOrder = NULL;

  for (set in 1:NSets)
  {
    SetSeedVectors = scale(SeedVectors[c(SetStarts[set]:SetEnds[set]), ]);
    SetLeaveOut = LeaveOut[, set];
    # Convert SetLeaveOut from boolean to a list of indices where it's TRUE
    SetMinCorr = rep(MinCorr, NMods);
    SetMaxCorr = rep(MaxCorr, NMods);
    SetLO = c(1:NMods)[SetLeaveOut];
    SimulatedData = SimulateExprData(SetSeedVectors, NGenes, ModProps,
                          MinCorr = SetMinCorr, MaxCorr = SetMaxCorr,
                          MinMaxCorr = MinMaxCorr, MaxMaxCorr = MaxMaxCorr, MaxCorrN0 = MaxCorrN0,
                          CorPower = CorPower,
                          BackgrNoise = BackgrNoise, LeaveOut = SetLO,
                          ModuleColors = ModuleColors,
			  NOrderedSmallLayers = NOrderedSmallLayers,
                          NUnorderedSmallLayers  = NUnorderedSmallLayers ,
                          AverageNGenesInSmallModule = AverageNGenesInSmallModule,
                          AverageExprInSmallModule = AverageExprInSmallModule,
                          InvDensityOfSmallModules = InvDensityOfSmallModules,
                          verbose = verbose-1, print.level = print.level+1);
    ExprData[[set]] = list(data = SimulatedData$ExprData);
    TrueColors = cbind(TrueColors, SimulatedData$TrueColors);
    TrueAllColors = cbind(TrueAllColors, SimulatedData$TrueAllColors);
    ColorOrder = cbind(ColorOrder, SimulatedData$ColorOrder);
  }
  list(ExprData = ExprData, TrueColors = TrueColors, TrueAllColors = TrueAllColors, ColorOrder = ColorOrder);
}

#--------------------------------------------------------------------------------------
#
# SimulateOrderedPCs
#
#--------------------------------------------------------------------------------------
# For given seed vectors, number of genes, an array of fractions of the number in each module and the
# noise level, create a simulated expression data and feed it through the usual module-finding code.
#
# Returns a list containing the module eigengenes found through the standard calculation as well as the
# eigengenes found by performing a PC decomposition using the true module colors.

SimulateOrderedPCs = function(SeedVectors, No.Genes, ModProps,
                              MinCorr = 0.5, MinMaxCorr = 0.5, MaxMaxCorr = 1, MaxCorrN0 = 100,
                              CorPower = 1, LeaveOut = NULL,
                              DissimilarityLevel = 1, BranchHeightCutoff = 0.997,
                              ModuleMinSize = 80, verbose = 4, print.level = 0)
{
  No.Samples = dim(SeedVectors)[[1]];

  SimulatedData = SimulateExprData(SeedVectors, No.Genes, ModProps, MinCorr = MinCorr,
                                   MinMaxCorr = MinMaxCorr, MaxMaxCorr = MaxMaxCorr, MaxCorrN0 = MaxCorrN0,
                                   CorPower = CorPower, LeaveOut = LeaveOut,
                                   verbose = verbose, print.level = print.level);

  PCs = GetOrderedPCs(SimulatedData$ExprData, DissimilarityLevel, BranchHeightCutoff, ModuleMinSize,
                      verbose, print.level);

  SPCs = list(PCs = PCs$PCs, OrderedPCs = PCs$OrderedPCs,  SeedVectors = SeedVectors, Network = PCs$Network,
              SimulatedColors = SimulatedData$TrueColors, ExprData = SimulatedData$ExprData)

  SPCs;

}

#--------------------------------------------------------------------------------------
#
# GetOrderedPCs
#
#--------------------------------------------------------------------------------------
# For given expression data, perform the standard network and ME calculations

GetOrderedPCs = function(ExprData, DissimilarityLevel = 1,
                         BranchHeightCutoff = 0.997, ModuleMinSize = 80,
                         verbose = 4, print.level = 0)
{
  No.Samples = dim(ExprData)[[1]];

  Subsets = rep(1, times = No.Samples);

  Network = GetNetwork(ExprData = ExprData, Subsets = Subsets, DegreeCut = 0, ProbeSelection = "all",
                       DissimilarityLevel = DissimilarityLevel,
                       BranchHeightCutoff = BranchHeightCutoff, ModuleMinSize = ModuleMinSize,
                       verbose = verbose, print.level = print.level);

  # Dynamic pruning of the clustering tree

  DynamicColors = as.character(cutreeDynamic(hierclust = Network$ClusterTree[[1]]$data, deepSplit = FALSE,
                          maxTreeHeight = BranchHeightCutoff, minModuleSize = ModuleMinSize));

  Network$Colors[,1] = DynamicColors;

  # Module Eigengenes organized in a list: PCs[[i]]$data is the dataframe of principal
  # components in subset i.

  PCs = NetworkModulePCs(ExprData, Network,
                         verbose = verbose, print.level = print.level);

  OrderedPCs = OrderPCs(PCs, GreyLast=TRUE, GreyName = "MEgrey");

  SPCs = list(PCs = PCs, OrderedPCs = OrderedPCs, Network = Network);

  SPCs;
}
