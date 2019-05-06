snp.frequency <- function(data, maf=TRUE) {
  
  ct<-apply(gt_matrix, 2, sum, na.rm=T)
  norm<-apply(!is.na(gt_matrix), 2, sum)
  
  freq<-ct/norm;
  
  if (maf==TRUE) {
    f.inv<-1-freq;
    freq<-apply(cbind(freq, f.inv), 1, min);
  }
  return(freq)
}

selection.statistics<-function(data, derived=TRUE, rho=-1, fixed=FALSE, ehh.position=-1, f.cut=0.05, f.half=0.1) {
  
  # number of sequences
  ns <- nrow(gt_matrix)
  # number of segregating sites
  s <- ncol(gt_matrix)
  
  # ns <-data$n.seq;
  # s <- data$l.seq;
  d <- dist(gt_matrix, method="manhattan");
  
  pwd = 2*sum(d)/(ns*(ns-1));
  
  f<-apply(gt_matrix, 2, sum);
  n.sing<-sum(f==1);
  if (derived == FALSE) n.sing<-n.sing + sum(f==(data$n.seq-1));
  
  con=vector(length=10);
  con[1]=sum(1/c(1:(ns-1)));
  con[2]=sum((1/c(1:(ns-1)))^2);
  con[3]=(ns+1)/(3*ns-3);
  con[4]=2*(ns*ns+ns+3)/(9*ns*ns-9*ns);
  con[5]=con[3]-1/con[1];
  con[6]=con[4]-(ns+2)/(ns*con[1])+con[2]/(con[1]*con[1]);
  con[7]=con[5]/con[1];
  con[8]=con[6]/(con[1]*con[1]+con[2]);
  con[9]=2*(ns*con[1]-2*(ns-1))/((ns-1)*(ns-2));
  con[10]=con[9]+(ns-2)/((ns-1)*(ns-1))+((2/(ns-1))*(1.5-(2*(con[1]+1/ns)-3)/(ns-2)-1/ns));
  con[11]=(ns*ns*con[2]/((ns-1)*(ns-1))+con[1]*con[1]*con[10]-2*ns*con[1]*(con[1]+1)/((ns-1)*(ns-1)))/(con[1]*con[1]+con[2]);       
  con[12]=ns/(ns-1)*(con[1]-ns/(ns-1))-con[11];
  
  if (derived==TRUE) {
    f<-snp.frequency(gt_matrix, maf=FALSE)*ns;
    h<-hist(f[f>0 & f<ns], breaks=c(0:(ns-1)), plot=F)$counts;
    theta.der<-2*sum(h*c(1:(ns-1))^2)/(ns*(ns-1));
  }
  
  if (derived == FALSE) {
    return(list(n.seq=ns, s=s, pi=round(pwd,2), theta.Watterson=round(data$l.seq/con[1],2), 
                theta.singleton = n.sing,
                Tajima.D= round((pwd-s/con[1])/sqrt(con[7]*s+con[8]*s*(s-1)),2), 
                FuLi.D = round((ns/(ns-1)*s-con[1]*n.sing)/sqrt(con[12]*s+con[11]*s*s),2)));
  }
  
  if (ehh.position>0 & fixed==FALSE) {
    pos<-which(data$pos==ehh.position);
    iehh<-calculate.iehh(data=data, snp.pos=pos, rho=rho, f.cut=f.cut, f.half=f.half);
    
    return(list(n.seq=ns, s=s, pi=round(pwd,2), theta.Watterson=round(data$l.seq/con[1],2),
                theta.derived = round(theta.der,2), theta.singleton = n.sing,
                Tajima.D= round((pwd-s/con[1])/sqrt(con[7]*s+con[8]*s*(s-1)),2), 
                FuLi.D = round((ns/(ns-1)*s-con[1]*n.sing)/sqrt(con[12]*s+con[11]*s*s),2),
                FayWu.H = round(pwd,2)-round(theta.der, 2),
                ehh.position=ehh.position, ehh.freq = sum(data$seq[,pos])/data$n.seq, iehh=iehh));
  }
  
  if (ehh.position<1 | fixed==TRUE) {
    
    return(list(n.seq=ns, 
                s=s, 
                pi=round(pwd,2), 
                theta.Watterson=round(data$l.seq/con[1],2), 
                theta.derived = round(theta.der,2), 
                theta.singleton = n.sing,
                Tajima.D= round((pwd-s/con[1])/sqrt(con[7]*s+con[8]*s*(s-1)),2), 
                FuLi.D = round((ns/(ns-1)*s-con[1]*n.sing)/sqrt(con[12]*s+con[11]*s*s),2),
                FayWu.H = round(pwd,2)-round(theta.der, 2)));
    
    
  }
  
}


# gt_matrix is a genotype matrix data frame that has markers as columns and samples as rows
# do not include column of strain names
selection.statistics<-function(data, derived=TRUE, rho=-1, fixed=FALSE, ehh.position=-1, f.cut=0.05, f.half=0.1) {
  
  # number of sequences
  ns <- nrow(gt_matrix)
  # number of segregating sites
  s <- ncol(gt_matrix)
  
  # ns <-data$n.seq;
  # s <- data$l.seq;
  d <- dist(gt_matrix, method="manhattan");
  
  pwd = 2*sum(d)/(ns*(ns-1));
  
  f<-apply(gt_matrix, 2, sum);
  n.sing<-sum(f==1);
  if (derived == FALSE) n.sing<-n.sing + sum(f==(data$n.seq-1));
  
  con=vector(length=10);
  con[1]=sum(1/c(1:(ns-1)));
  con[2]=sum((1/c(1:(ns-1)))^2);
  con[3]=(ns+1)/(3*ns-3);
  con[4]=2*(ns*ns+ns+3)/(9*ns*ns-9*ns);
  con[5]=con[3]-1/con[1];
  con[6]=con[4]-(ns+2)/(ns*con[1])+con[2]/(con[1]*con[1]);
  con[7]=con[5]/con[1];
  con[8]=con[6]/(con[1]*con[1]+con[2]);
  con[9]=2*(ns*con[1]-2*(ns-1))/((ns-1)*(ns-2));
  con[10]=con[9]+(ns-2)/((ns-1)*(ns-1))+((2/(ns-1))*(1.5-(2*(con[1]+1/ns)-3)/(ns-2)-1/ns));
  con[11]=(ns*ns*con[2]/((ns-1)*(ns-1))+con[1]*con[1]*con[10]-2*ns*con[1]*(con[1]+1)/((ns-1)*(ns-1)))/(con[1]*con[1]+con[2]);       
  con[12]=ns/(ns-1)*(con[1]-ns/(ns-1))-con[11];
  
  if (derived==TRUE) {
    f<-snp.frequency(gt_matrix, maf=FALSE)*ns;
    h<-hist(f[f>0 & f<ns], breaks=c(0:(ns-1)), plot=F)$counts;
    theta.der<-2*sum(h*c(1:(ns-1))^2)/(ns*(ns-1));
  }
  
  # Wattersons_Theta <- round(s/con[1],2) # Watterson's Theta
  # Average_PWD <- round(pwd,2) # average pairwise differences
  # TajimasD <- round((pwd-s/con[1])/sqrt(con[7]*s+con[8]*s*(s-1)),2) # Tajima's D
  # 
  fq<-apply(gt_matrix, 2, sum, na.rm=T)
  n.na<-apply(is.na(gt_matrix), 2, sum)
  n.sing <- sum(fq==1)+sum(fq==(gt_matrix-n.na-1))
  
  # FuLiF <- round((ns/(ns-1)*s-con[1]*n.sing)/sqrt(con[12]*s+con[11]*s*s),2) # Fu and Li's D
  
  gene_summary <- data.frame(Wattersons_Theta = round(s/con[1],2),
                             Average_PWD = round(pwd,2) ,
                             TajimasD = round((pwd-s/con[1])/sqrt(con[7]*s+con[8]*s*(s-1)),2),
                             Singletons = n.sing,
                             FuLiF = round((ns/(ns-1)*s-con[1]*n.sing)/sqrt(con[12]*s+con[11]*s*s),2))
  return(gene_summary)
}
