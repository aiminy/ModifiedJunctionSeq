
#' JunctionSeq:::MatchStaisticsWithPvalue(Re.example.statistics, re.example.gene.based)

MatchStaisticsWithPvalue <- function(Re.example.statistics, re.example.gene.based)
  {
  
  #install.packages("acepack")
  #library(acepack)
  
  #TWOPI <- 8*atan(1)
  #x <- runif(200,0,TWOPI)
  #y <- exp(sin(x)+rnorm(200)/2)
  #a <- ace(x,y)
  
  #a <- ace(x,y)
  
  temp1=fData(Re.example.statistics)
  temp2=pData(re.example.gene.based)
  
  index<-paste0(temp2$geneID,temp2$mostSigID)
 
  temp4<-cbind(temp2,index)
  
  index<-paste0(temp1$geneID,temp1$countbinID)
  temp5<-cbind(temp1,index)
  
  #temp5[which(temp5$index %in% "ENSG00000223972+ENSG00000227232+ENSG00000243485+ENSG00000274890E004"),]
  
  temp6<-merge(temp4,temp5,by="index",sort = FALSE, suffixes = c(".p",".s"))
  
  par(mfrow=c(3,1))
  hist(temp6$pvalue,xlab="test statistics",main="Histogram of test statistics")
  hist(temp6$numExons,xlab="number of exons", main="Histogram of number of exons")
  
  plot(temp6$pvalue~temp6$numExons,xlab="#exon",ylab="Statistics")
  
  x=temp6$numExons
  y=temp6$pvalue
  a <- ace(x,y)
  
  
  par(mfrow=c(3,1))
  plot(a$y,a$ty)  # view the response transformation
  plot(a$x,a$tx)  # view the carrier transformation
  plot(a$tx,a$ty) # examine the linearity of the fitted model
  
  cor(a$tx,a$ty)
  
  
  cor(temp6$pvalue,temp6$numExons)
  #identify relationship between statistics and number of exon
  
  plot(log(temp6$pvalue)~log(temp6$numExons))
  cor(sqrt(temp6$pvalue),temp6$numExons,use = "complete.obs")
  cor(log10(abs(temp6$pvalue)),log10(temp6$numExons),use="complete.obs")
  cor(exp(temp6$pvalue),exp(temp6$numExons),use="complete.obs")
  
  
  Y=a$ty
  X=a$tx
  
  m<-lm(y~x)
  m<-loess(y~x,span=0.3)
  #m<-nls(y~x)
  
  
  
  #YY=mean(y)+residuals(fm)
  
  YY=mean(y)+residuals(m)
  
  y.re<-cbind(y,residuals(m))
 
  index.n<-which(y.re[,1]<0)
  index.p<-which(y.re[,1]>=0)
  
  ynm<-mean(y.re[index.n,1])
  ypm<-mean(y.re[index.p,1])
  
  y.re[index.n,2]<-ynm+y.re[index.n,2]
  y.re[index.p,2]<-ypm+y.re[index.p,2]
  
  
  YYY<-y.re[,2]
    
  hist(YYY)
  
  hist(YY)
  hist(residuals(m))
  hist(y)
  
  plot(Y~X)
  plot(YY~x)
  plot(YY~y)
  
  round(cor(Y,X),3)
  round(cor(YY,x),3)
  round(cor(y,x),3)
  
  p.new.2<-pchisq(YYY, df = 1, lower.tail = FALSE)
  
  p.new<-pchisq(YY, df = 1, lower.tail = FALSE)
  p.old<-pchisq(y, df = 1, lower.tail = FALSE)
  hist(p.new)
  hist(p.old)
  
  hist(p.new.2)
  
  temp7<-cbind(temp6,p.old,p.new)
  
  par(mfrow=c(3,1))
  plot(x,y)

  xplot=seq(2,690,1)
  
  for(degree in 1:4){
    fm= lm(y~poly(x,degree))
    assign(paste("y",degree,sep="."),fm) 
    # this isn't needed, but handy
    lines(xplot,predict(fm,data.frame(x=xplot)), col=degree)
  }
  
  fm= lm(y~poly(x,2.999999))
  YY=mean(y)+residuals(fm)
  p.new<-pchisq(YY, df = 1, lower.tail = FALSE)
  
  temp7<-cbind(temp6,p.old,p.new)
  temp8<-temp7[,c(2:7,44,9,44,11:18)]
   
  colnames(temp8)<-colnames(pData(re.example.gene.based))
  pData(re.example.gene.based)<-temp8
  rownames(pData(re.example.gene.based))<-pData(re.example.gene.based)[,1]
   
  re.example.gene.based.testable.reformat.2<-pData(re.example.gene.based)
  
  #UseLogistic2CKBias(re.example.gene.based.testable.reformat.2,0.05,type="exon")
  
  UseLogistic2CKBias(re.example.gene.based.testable.reformat.2)
  
}