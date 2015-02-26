#' A function that creates a heatmap with column indicators on top of each column.
#' @description This function creates a nice heatmap for publication. It puts on top of each column a color indicator of what sample type it is, does clustering, and adds a table in the right margin of the figure with logFoldChanges and respective significance indicators (stars).
#' @param mat numeric matrix to cluster. Usually th raw expression data.
#' @param coefs estimates to be drawn tin the tabl (eg lgFCHs)
#' @param fdrs p values or fdr indicating significance of the estimates; should be the same dim as coefs
#' @param colfc factor to add colored legend to the heatmap
#' @param facsepr factor to separate columns in the heatmaps
#' @param hmeth aglomeration strategy 
#' @param dmeth aglomeration distance to cluster
#' @param ColD logical TRUE if cluster by rows and columns
#' @param sn edited names for the rows (if different than the rownames) eg symbols for probesets 
#' @param gn gene name
#' @keywords heatmap , expression vizualisation
#' @examples 
#' doClusterTableFCH_fromMatrix(mat,ps=rownames(mat),colfac,facsepr=NULL,hmeth='average',dmeth=cor.dist,coefs,fdrs,main="", ColD=FALSE, ss=ps, gn=NULL, breaks=NULL,cexg=1,margins=c(5,20) )

doClusterTableFCH_fromMatrix<-function(mat,ps=rownames(mat),colfac,facsepr=NULL,hmeth='average',dmeth=cor.dist,coefs,fdrs,main="",
                                       ColD=FALSE, ss=ps, gn=NULL, breaks=NULL,cexg=1,margins=c(5,20)){
  #mat: numeric matrix to cluster
  #coefs: estimates to be drawn tin the table (eg lgFCHs)
  #fdrs: p values or fdr indicating significance of the estimates; should be the same dim as coefs
  #colfc: factor to add colored legend to the heatmap
  #facsepr: factor to separate columns in the heatmaps
  #hmeth, dmeth: aglomeration strategy and distance to cluster
  #ColD; logical TRUE if cluster by rows and columns
  #sn edited names for the rows (if different than the rownames) eg symbols for probesets
  # gn: gene name, 
  
  require(stringr)
  require(weights)
  
  mat<-mat[ps,]
  maxss<-max(sapply(ss,nchar))
  adspace<-function(x,n){paste(x,substr(" ---------------------",1,n-nchar(x)),sep="")}
  ss<-paste(sapply(ss,adspace,maxss+3),":  ",sep='')
  
  ADDzero<-function(x){
    if(grepl(pattern="\\.",x=x)==TRUE){
      # print(".")
      a1<-unlist(strsplit(x,"\\."))[1]
      a2<-unlist(strsplit(x,"\\."))[2]
      if (is.na(a2)==T){a2<-"0"}
      b1<-paste(str_dup(" ",4-nchar(a1)),a1,sep="")
      b1<-a1
      b2<-paste(a2,str_dup("0",2-nchar(a2)),sep="")
      out<-paste(b1,b2,sep=".")
    } else{
      out <- paste(paste(x,".00",sep=""),sep="")
    }
    return(out)
  }
  adzero<-function(xv){
    out<-sapply(xv,ADDzero)
    nmax<-max(sapply(out,nchar))
    sapply(out,function(x,nmx){paste(str_dup(" ",nmx-nchar(x)),x,sep="")}, nmax)
  }
  
  reformatps<-function(p){as.numeric(format(p,digit=3,drop0trailing = TRUE))}
  transformfch<-function(lgfch){fch<-sign(lgfch)*(2^abs(lgfch)); return(fch)}
  coef1<-apply(coefs[ps,],2,function(x){adzero(as.character(signif(transformfch(x),3)))})
  rownames(coef1)<-ps
  fdrs1<-apply(fdrs[ps,],2,function(x){ifelse(x<0.0001,signif(x,1),round(x,4))})
  fdrs2<-apply(fdrs1,2,starmaker,p.levels=c(.01, .05, .1), symbols=c("**", "*", "+"))
  rownames(fdrs2)<-ps
  adSpace <- function(x){
    out <- paste(x,str_dup(" ",4-nchar(x)),sep="")
  }
  fdrs2<-apply(fdrs2,2,adSpace)
  coef2<-coef1 #apply(coef1,2,function(x){paste(x,str_dup(" ",),sep="")})
  print(head(coef2))
  print(head(fdrs2))
  a<-DoHeatmap(mat, colfac=colfac, symb=ss, dmeth=dmeth, hmeth=hmeth, cex.genes=cexg, ColD=ColD, main=main, margins=c(5,10),breaks=breaks)
  
  Tab<-cbind(Symbol=ss[a$rowInd], Desc=substr(gn[a$rowInd],1,40))
  for ( i in c(1:ncol(coef2))){
    Tab<-cbind(Tab,paste(coef2[a$rowInd,i],fdrs2[a$rowInd,i],sep=''))
  }
  Tab<-print.table(Tab)
  ssTab<-apply(Tab[,-2],1,function(x){stackchar(x,sep="")})
  mat2<-mat[a$rowInd,];
  rownames(mat2)<-ssTab
  ssTab_bold<-do.call(expression,sapply(as.character(ssTab),function(.x){substitute(bold(.x),list(.x=.x))}))
  par(family='mono')
  a<-DoHeatmap(mat2, colfac=colfac, facsepr=facsepr, symb=ssTab_bold, dmeth=dmeth, hmeth=hmeth, cex.genes=cexg, ColD=ColD, main=main ,margins=margins, breaks=breaks)
  return(invisible(Tab))
}