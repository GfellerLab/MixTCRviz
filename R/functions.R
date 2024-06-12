
#####
# Define some function
#####

build_stat <- function(es, chain.list=c("TRA","TRB"), sp="HomoSapiens", comp.VJL=0){

  L <- list()
  countL <- list()
  countV <- list()
  countJ <- list()
  countV.L <- list()
  countJ.L <- list()
  countVJ <- list()
  countCDR1 <- list()
  countCDR2 <- list()
  countCDR3.L <- list()
  countVJ.L <- list()
  if(comp.VJL==1){
    countCDR3.VJL <- list()
    countL.VJ <- list()
  }

  for(chain in chain.list){

    seg <- c(paste(chain,"V", sep=""),paste(chain,"J", sep=""))
    cdr3 <- paste("cdr3_",chain,sep="")

    countV[[chain]] <- table(es[,seg[1]])
    countJ[[chain]] <- table(es[,seg[2]])

    countL[[chain]] <- table(nchar(es[,cdr3]))
    L[[chain]] <- as.numeric(names(countL[[chain]]))
    if(length(countL[[chain]])>0){
      names(countL[[chain]]) <- paste("L",names(countL[[chain]]),sep="_")
    }  
    countVJ[[chain]] <- table(es[,seg[1]], es[,seg[2]])

    countV.L[[chain]] <- list()
    countJ.L[[chain]] <- list()
    countVJ.L[[chain]] <- list()
    for(lg in L[[chain]]){
      ind <- which(nchar(es[,cdr3])==lg)
      lg.c <- paste("L",lg,sep="_")
      countV.L[[chain]][[lg.c]] <- table(es[ind,seg[1]])
      countJ.L[[chain]][[lg.c]] <- table(es[ind,seg[2]])
      countVJ.L[[chain]][[lg.c]] <- table(es[ind,seg[1]],es[ind,seg[2]])
    }
        
    if(comp.VJL==1){  
      for(V in names(countV[[chain]])){
        indv <- which(es[,seg[1]]==V)
        for(J in names(countJ[[chain]])){
          indj <- which(es[indv,seg[2]]==J)
          ind <- indv[indj]
          s <- paste(V,J, sep="_")
          if(length(ind)>0){
            countL.VJ[[chain]][[s]] <- table(nchar(es[ind,cdr3]))
            if(length(countL.VJ[[chain]][[s]])>0){
              names(countL.VJ[[chain]][[s]]) <- paste("L", names(countL.VJ[[chain]][[s]]),sep="_")
            }
            countCDR3.VJL[[chain]][[s]] <- count_aa(es[ind,cdr3], keep.gap=0)
          } else {
            countL.VJ[[chain]][[s]] <- table(NA)
            countCDR3.VJL[[chain]][[s]] <- table(NA)
          }
        }
      }
    }
    if(length(countV[[chain]])>0){
      countCDR1[[chain]] <- count_aa(cdr123[[sp]][[chain]][es[,seg[1]],"CDR1"], keep.gap=1)
      countCDR2[[chain]] <- count_aa(cdr123[[sp]][[chain]][es[,seg[1]],"CDR2"], keep.gap=1)
    }
    countCDR3.L[[chain]] <- count_aa(es[,cdr3], keep.gap=0)


  }
  if(comp.VJL==0){
    count <- list(L, countL, countV, countJ, countV.L, countJ.L, countCDR1, countCDR2, countCDR3.L, countVJ, countVJ.L)
    names(count) <- c("L", "countL", "countV", "countJ", "countV.L", "countJ.L", "countCDR1", "countCDR2", "countCDR3.L", "countVJ", "countVJ.L")
  } else {
    count <- list(L, countL, countV, countJ, countV.L, countJ.L, countL.VJ, countCDR1, countCDR2, countCDR3.L, countCDR3.VJL, countVJ, countVJ.L)
    names(count) <- c("L", "countL", "countV", "countJ", "countV.L", "countJ.L", "countL.VJ", "countCDR1", "countCDR2", "countCDR3.L", "countCDR3.VJL", "countVJ", "countVJ.L")
    
  }
  return(count)
}

#Compute the counts of each aa at each position.
#Would be better to have the option of treating gaps in CDR1/2 as separate amino acids, while missing data can be treated as 'unspecific'
count_aa <- function(cdr.seq, keep.gap=0){

  if(keep.gap == 0){    #All gaps, including "g" are treated as unspecific data. This can be useful for visualisation.
    tgap <- c(gap,"g")
    taa.list <- aa.list
  } else {   # "x" are treated as additonal aa, other 'gaps' are discarded. This is more correct for modelling CDR1/CDR2 loops.
    tgap <- c()
    taa.list <- c(aa.list,"g")
  }

  #First get the list of length
  l.seq <- nchar(cdr.seq)
  L <- sort(unique(l.seq))
  L <- L[L>0]  # This ensures that the cases of empty sequences are never counted.
  m.list <- list()
  for(lg in L){
    ind <- which(l.seq==lg)
    tcdr.seq <- cdr.seq[ind]
    m <- matrix(0, nrow=length(taa.list), ncol=lg)
    rownames(m) <- taa.list
    for(p in 1:lg){
      s <- substr(tcdr.seq,p,p)
      tb <- table(s)
      for(a in names(tb)){
        if(keep.gap==0){
          if(a %in% tgap){
            m[,p] <- m[,p] + tb[a]/20
          } else {
            m[a,p] <- m[a,p]+tb[a]
          }
        } else {  #Here "x" are treated as a separate amino acid and missing data (e.g., "*", "X",etc.) are not included
          if(a %in% gap==F){
            m[a,p] <- m[a,p]+tb[a]
          }
        }
      }
    }
    lc <- paste("L",lg,sep="_")
    m.list[[lc]] <- m
  }
  return(m.list)
}

#This assumes the matrix includes gaps (21 rows)
build_cdr12_motif <- function(cdr.seq, keep.gap=0){

  if(keep.gap==0){
    g <- cdr.seq["g",]
    for(a in aa.list){
      cdr.seq[a,] <- cdr.seq[a,]+g/length(aa.list)
    }
    cdr.seq <- cdr.seq[aa.list,]
  }
  cdr.seq <- scale(cdr.seq, center=F, scale=colSums(cdr.seq))
  return(cdr.seq)
}


# Compute the weighted average of motifs for different V/J usage
# This is important since the choice of V/J has huge influence on the CDR3 motif.
# By considering V/J usage in baseline, it helps identifying what is specific to the input TCR
weighted_countCDR3 <- function(countCDR3.VJL.baseline, countVJ.L.es){

  countCDR3 <- list()
  L <- names(countVJ.L.es)
  for(lg.c in L){
    
    lg <- as.numeric(unlist(strsplit(lg.c, split="_"))[2])

    countCDR3[[lg.c]] <- matrix(0, nrow=N.aa, ncol=lg)
    rownames(countCDR3[[lg.c]]) <- aa.list
    if(sum(countVJ.L.es[[lg.c]])>0){
      w <- countVJ.L.es[[lg.c]]/sum(countVJ.L.es[[lg.c]])
      for(v in rownames(countVJ.L.es[[lg.c]])){
        for(j in colnames(countVJ.L.es[[lg.c]])){
          if(w[v,j]>0){
            s <- paste(v,j,sep="_")
            if(lg.c %in% names(countCDR3.VJL.baseline[[s]])){
              if(length(countCDR3.VJL.baseline[[s]][[lg.c]]) > 0){
                countCDR3[[lg.c]] <- countCDR3[[lg.c]]+scale(countCDR3.VJL.baseline[[s]][[lg.c]], center=F, scale=colSums(countCDR3.VJL.baseline[[s]][[lg.c]]))*w[v,j]
              }
            }
          }
        }
      }
    }
    if(colSums(countCDR3[[lg.c]])[1]!=0){
      countCDR3[[lg.c]] <- scale(countCDR3[[lg.c]],center=F, scale=colSums(countCDR3[[lg.c]]))
    } else {
      countCDR3[[lg.c]][1:N.aa,1:lg] <- 0.05 #This is the case if the were no data for all V-J for this length in the repertoire. Ideally, we should indicate this in the title of the plot
    }
  }
  return(countCDR3)
}

# Compute the weighted average of CDR3 length for the observed V/J usage
# This is important since CDR3 length is primarily determined by the length of the V and J segments
weighted_countL <- function(countL.VJ.baseline, countVJ.es){


  countL <- rep(0,times=Lmax-Lmin+1)
  names(countL) <- paste("L", Lmin:Lmax, sep="_")

  countVJ.es <- countVJ.es/sum(countVJ.es)
  for(v in rownames(countVJ.es)){
    for(j in colnames(countVJ.es)){
      if(countVJ.es[v,j]>0){
        s <- paste(v,j,sep="_")
        if(s %in% names(countL.VJ.baseline)){
         if(length(countL.VJ.baseline[[s]])>0){
           countL.VJ.baseline[[s]] <- countL.VJ.baseline[[s]]/sum(countL.VJ.baseline[[s]])*countVJ.es[v,j]
           for(lc in names(countL.VJ.baseline[[s]])){
             countL[lc] <- countL[lc]+countL.VJ.baseline[[s]][lc]
           }          
          }
        }
      }
    }
  }
  return(countL)
}



#Infer the short name for the mhc from the long name (HLA-A*02:01 -> A0201)
#The input should be a matrix/dataframe with the MHC columns and a specices column, since it only works for human
#This is only needed for the MixTCRpred1.0 input.
find_mhc <- function(m){
  mhc <- unlist(lapply(1:dim(m)[1], function(x){
    if(m[x,"MHC"]=="HLA-A:01") {h <- "A0101"}
    else{
      if(m[x,"species"] == "HomoSapiens"){
        a <- unlist(strsplit(m[x,"MHC"], split="-", fixed=T))
        h <- gsub('*', '', a[2], fixed=T)
        h <- gsub(':', '', h, fixed=T)
        h <- gsub("/", "_", h, fixed=T)
      } else {
        h <- m[x,"MHC"]
      }
    }
    return(h)
  }))
}

# Plot the comparison between V or J usage.
# We should have different criteria for selecting the labels for the ES versus Pred
plotVJ <- function(count.es, count.rep, info, comp.baseline){

  #count.es is on the Y axis, count.rep on the X
  gene <- info[1]
  n <- sum(count.es)
  if(comp.baseline==0){n.rep <- sum(count.rep)}
  v <- c(names(count.rep), names(count.es))
  nm <- unique(v)
  type1 <- info[2]
  type2 <- info[3]
  cn <- c(type1, type2)
  count <- matrix(0, nrow=length(nm), ncol=2)
  rownames(count) <- nm
  colnames(count) <- cn
  for(i in 1:length(count.es)){ count[names(count.es[i]),type1] <- count.es[i] }
  for(i in 1:length(count.rep)){ count[names(count.rep[i]),type2] <- count.rep[i] }
  count <- scale(count, center=F, scale=colSums(count))
  count.df <- data.frame(count)
  colnames(count.df) <- c("Y","X")

  lim.y <- max(count[,c(type1)] )*1.3
  lim.x <- max(count[,c(type2)] )*1.3

  #lim.y <- round(min( max(count[,c(type1)] )+0.1, max( count[,c(type1)] )*1.4),1)
  #lim.x <- min(max(count[,c(type2)])+0.1, max(count[,c(type2)])*1.4)
  #if(lim.x < 0.1){lim.x <- round(lim.x+0.005,2)}
  #else { lim.x <- round(lim.x+0.05,1) }

  label <- nm
  label[count.df[,"Y"] < 0.05 & count.df[,"X"] < 0.05] <- NA
  ratio <- count.df[,"Y"]/(count.df[,"X"]+0.001)
  label[which(ratio < 1.5 & count.df[,"X"] < 0.3 & count.df[,"Y"] < 0.3) ] <- NA
  #If there are too many labels, show only those with FC > 2
  if(length(which(!is.na(label)))>8){
    label[which( (ratio < 2) & count.df[,"Y"]<0.1)] <- NA
  }

  ylab <- paste(type1," (",n,")", sep="")

  if(comp.baseline==1){  xlab <- type2 } else {xlab <- paste(type2," (",n.rep,")", sep="")}

  #Plot the comparison between input and repertoires
  count.plot <- list()
  if(length(count.es)>0){
    count.plot <- ggplot(count.df, aes(x=X, y=Y, label=label)) +
      geom_point() + geom_abline(col="orange",linetype="dashed",linewidth=1) +
      ggtitle(gene) +
      xlim(0, lim.x) + ylim(0,lim.y) +
      theme(plot.title = element_text(size = 14, hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=14)) +
      geom_label_repel(size = 3, nudge_y=0.02, box.padding = 0.15) +
      xlab(xlab) + ylab(ylab)
  } else {count.plot <- ggplot()}

  return(count.plot)

}


correct.VJnames <- function(es.all, name.list, segment.list){


  for(sp in species.list){
    for(s in segment.list){

      ind <- which(es.all[,s] %in% name.list[[sp]]==F & es.all[,"species"]==sp)
      if(length(ind)>0){
        nm <- strsplit(es.all[ind,s], split="*", fixed=T)
        gene <- unlist(lapply(nm, function(x){x[1]}))
        allele <- unlist(lapply(nm, function(x){x[2]}))

        ga <- unlist(lapply(1:length(gene), function(x){ clean.name.allele(gene[x],allele[x],sp)}))

        #Check the cases where the segment was not NA, but was put to NA (i.e., mapping of gene name failed)
        i <- which(es.all[ind,s] != "" & is.na(ga)==T)
        if(length(i)>0){
          print(paste(s, " gene names not in IMGT which could not be corrected:"))
          print(names(table(es.all[ind[i],s])))
        }
        es.all[ind,s] <- ga
      }
    }
  }
  return(es.all)
}


plotLD <- function(countL.es, countL.rep,info, plot.oneline){

 
  L.es <- as.numeric(lapply(names(countL.es), function(x){unlist(strsplit(x,split="_"))[2]}))
  L.rep <- as.numeric(lapply(names(countL.rep), function(x){unlist(strsplit(x,split="_"))[2]}))
  
  L.all <- sort(unique(L.es, L.rep))
 
  ct <- 1
  ld.es <- c()
  ld.rep <- c()

  for(lg in L.all){
    lc <- paste("L",lg,sep="_")
    if(!is.na(countL.es[lc])){ ld.es[ct] <- countL.es[lc]} else {ld.es[ct] <- 0}
    if(!is.na(countL.rep[lc])){ld.rep[ct] <- countL.rep[lc]} else {ld.rep[ct] <- 0}
    ct <- ct+1
  }
  if(sum(ld.es)>0){
    ld.es <- ld.es/sum(ld.es)
  }
  if(sum(ld.rep)>0){
    ld.rep <- ld.rep/sum(ld.rep)
  }
  
  ##########
  #Plot the comparison for length distribution
  ##########

  v1 <- c(L.all,L.all);
  v2 <- c(ld.es, ld.rep);
  v3 <- c( rep(info[2], length(L.all)), rep(info[3], length(L.all))) ;
  ld.df <- data.frame(v1,v2,v3)
  ld.df$v3 <- factor(ld.df$v3, levels=c(info[2], info[3]))

  legend.size <- 12
  if(plot.oneline==1){
    if(nchar(info[2])>23){legend.size=11}
    if(nchar(info[2])>25){legend.size=10}
  }
  ld.plot <-  ggplot(ld.df) + geom_point(aes(x=v1, y=v2, color=v3)) + #ggtitle(paste("N =",sum(countL.es))) +
    geom_line(aes(x=v1, y=v2, color=v3)) +
    theme(legend.key.size = unit(0.2, 'cm'), legend.position="top", legend.title=element_blank(),  legend.text=element_text(size=legend.size)) +
    xlab(paste("Length_CDR3",info[1],sep="")) + ylab("") +
    guides(color = guide_legend(nrow = 2)) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14), plot.title = element_text(size=15,hjust = 0.5))

  return(ld.plot)

}


plotCDR3 <- function(countL.es, countL.rep, countCDR3.es, countCDR3.rep, info, comp.baseline, plot.oneline=0, plot.logo.length=0, plot.cdr3.subtract.baseline=0){
  
  L.es <- as.numeric(lapply(names(countL.es), function(x){unlist(strsplit(x,split="_"))[2]}))
  L.rep <- as.numeric(lapply(names(countL.rep), function(x){unlist(strsplit(x,split="_"))[2]}))
  
  L.TR <- intersect(L.es,L.rep)
  
  
  
  pwm.rep <- list()
  pwm.es <- list()
  
  logo.CDR3.L.es <- list()
  logo.CDR3.L.rep <- list()
  
  if(length(L.TR)>0){
    
    tl <- countL.es[paste("L",L.TR,sep="_")]
    lmax.es <- as.numeric(unlist(strsplit(names(tl[which.max(tl)]), split="_"))[2])
    lmax <- lmax.es
  
  
  if(plot.oneline==1){
    if(lmax<15){
      axis.size.max <- 8
    } else {
      axis.size.max <- 7
    }
    title.size <- 11
  } else {
    axis.size.max <- 10
    title.size <- 12
  }
  ylab <- ""
  
  
  if(plot.logo.length==0){
    L.TR <- c(lmax)
  }
  
  for(l in L.TR){
    
    lc <- paste("L",l,sep="_")
    
    pwm.es[[lc]] <- scale(countCDR3.es[[lc]], center=F, scale=colSums(countCDR3.es[[lc]]))
    pwm.rep[[lc]] <- scale(countCDR3.rep[[lc]], center=F, scale=colSums(countCDR3.rep[[lc]]))
    
    if(plot.cdr3.subtract.baseline==1){
      
      #Compute the logo representing the difference in frequencies renormalized by information content
      #This is currently not supported by ggseqlogoMOD...
      #Compute the matrix representing the size of each letter
      size.es <- apply(pwm.es[[lc]], 2, function(x){ ind <- which(x!=0); IC <- log(N.aa)/log(2)+sum(x[ind]*log(x[ind])/log(2)); return(IC*x) })
      size.rep <- apply(pwm.rep[[lc]], 2, function(x){ ind <- which(x!=0); IC <- log(N.aa)/log(2)+sum(x[ind]*log(x[ind])/log(2)); return(IC*x) })
      x.norm <- size.es-size.rep
      y.inc <- 1
    } 
    if(plot.cdr3.subtract.baseline==2){
      
      #Compute the logo based on normalised fold-change
      rc <- max(5,0.1*countL.es[[lc]])
      #rc <- 5
      pseudo <- 1/N.aa*rc/countL.es[[lc]]
      x.baseline <- pwm.rep[[lc]]+pseudo
      x.es <- pwm.es[[lc]]+pseudo
      x.norm <- x.es/x.baseline
      x.norm <- scale(x.norm,center = F, scale=colSums(x.norm))
      y.inc <- 4
    }
    
    title <- paste("CDR3", info[1],"_",l," ",info[2], " (",countL.es[[lc]],")", sep="")
    logo.CDR3.L.es[[lc]] <- ggseqlogoMOD(data=pwm.es[[lc]], additionaAA=additionalAA,  axisTextSizeX = 12, axisTextSizeY = 8) +
      labs(title=title) + ylab(ylab) + theme(plot.title=element_text(size=15, hjust=0.5))
    
    if(plot.cdr3.subtract.baseline==0){
      title.baseline <- paste("CDR3", info[1],"_",l," ",info[3], sep="")
    } else if(plot.cdr3.subtract.baseline==1){
      title.baseline <- paste("CDR3", info[1],"_",l," subtract ",info[3], sep="")
    } else if(plot.cdr3.subtract.baseline==2){
      title.baseline <- paste("CDR3", info[1],"_",l," renorm ",info[3], sep="")
    }
    
    
    if(comp.baseline==0){title.baseline <- paste(title.baseline, " (",countL.rep[[lc]],")", sep="")} else {title.baseline <- title.baseline}
    
    if(plot.cdr3.subtract.baseline==0){
      logo.CDR3.L.rep[[lc]] <- ggseqlogoMOD(data=pwm.rep[[lc]], additionaAA=additionalAA,  axisTextSizeX = 12, axisTextSizeY = 8) +
        labs(title=title.baseline) + ylab(ylab) + theme(plot.title=element_text(size=15, hjust=0.5))
    } else if(plot.cdr3.subtract.baseline==1){
      y.min <- min(apply(x.norm, 2, function(x){ sum(x[x<0]) }))
      y.max <- max(apply(x.norm, 2, function(x){ sum(x[x>0]) }))
      y.min <- max(-log(N.aa)/log(2), y.inc*y.min)
      y.max <- log(N.aa)/log(2) # min(log(N.aa)/log(2), y.inc*y.max)
      logo.CDR3.L.rep[[lc]] <- ggseqlogo(data=x.norm, method='custom') +
        labs(title=title.baseline) + ylim(y.min,y.max) + ylab(ylab) +
        theme(plot.title=element_text(size=title.size, hjust=0.5)) + theme(legend.position = 'none')
    } else if(plot.cdr3.subtract.baseline==2){
      IC.max <- max(unlist(apply(x.norm, 2, function(x){ ind <- which(x!=0); IC <- log(N.aa)/log(2)+sum(x[ind]*log(x[ind])/log(2)); return(IC) })))
      y.max <- min(IC.max*y.inc, log(N.aa)/log(2))
      logo.CDR3.L.rep[[lc]] <- ggseqlogoMOD(data=x.norm, additionaAA=additionalAA,  axisTextSizeX = 12, axisTextSizeY = 8, ylim=c(0, y.max)) +
        labs(title=title.baseline) + ylab(ylab) + theme(plot.title=element_text(size=15, hjust=0.5)) 
    }
    #For the special case where l==lmax, build the logo with different graphical parameters, depending on the plot.online
    #So far, we redo everything, since the graphical outline has to be a little bit different,
    # but this is not optimal since any change has to be performed multiple times
    if(l==lmax){
      title <- paste("CDR3", info[1],"_",l," ",info[2], " (",countL.es[[lc]],")", sep="")
      if(plot.cdr3.subtract.baseline==0){
        title.baseline <- paste("CDR3", info[1],"_",l," ",info[3], sep="")
      } else if(plot.cdr3.subtract.baseline==1){ 
        title.baseline <- paste("CDR3", info[1],"_",l," subtract ",info[3], sep="")
      } else if(plot.cdr3.subtract.baseline==2){ 
        title.baseline <- paste("CDR3", info[1],"_",l," renorm ",info[3], sep="")
      }
      if(comp.baseline==0){title.baseline <- paste(title.baseline, " (",countL.rep[[lc]],")", sep="")}
      
      if(plot.oneline==1 & (nchar(title)>26 | nchar(title.baseline)>26)){
        title <- paste("CDR3", info[1],"_",l," ",info[2], "\n(",countL.es[[lc]],")", sep="")
        if(plot.cdr3.subtract.baseline==0){
          title.baseline <- paste("CDR3", info[1],"_",l," ",info[3],"\n", sep="")
        } else if(plot.cdr3.subtract.baseline==1){ 
          title.baseline <- paste("CDR3", info[1],"_",l," subtract \n",info[3], sep="")
        } else if(plot.cdr3.subtract.baseline==2){ 
          title.baseline <- paste("CDR3", info[1],"_",l," renorm \n",info[3], sep="")
        }
        if(comp.baseline==0){title.baseline <- paste(title.baseline, "(",countL.rep[[lc]],")", sep="")}
      }
      
      logo.CDR3.L.es.max <- ggseqlogoMOD(data=pwm.es[[lc]], additionaAA=additionalAA,  axisTextSizeX = axis.size.max, axisTextSizeY = 8) +
        labs(title=title) + ylab(ylab) + theme(plot.title=element_text(size=title.size, hjust=0.5))
      
      if(plot.cdr3.subtract.baseline==0){
        logo.CDR3.L.rep.max <- ggseqlogoMOD(data=pwm.rep[[lc]], additionaAA=additionalAA,  axisTextSizeX = axis.size.max, axisTextSizeY = 8) +
          labs(title=title.baseline) + ylab(ylab) + theme(plot.title=element_text(size=title.size, hjust=0.5))
      } else if(plot.cdr3.subtract.baseline==1){
        y.min <- min(apply(x.norm, 2, function(x){ sum(x[x<0]) }))
        y.max <- max(apply(x.norm, 2, function(x){ sum(x[x>0]) }))
        y.min <- max(-log(N.aa)/log(2), 2*y.min)
        y.max <- log(N.aa)/log(2) # min(log(N.aa)/log(2), 2*y.max)
        logo.CDR3.L.rep.max <- ggseqlogo(data=x.norm, method='custom') +
            ggtitle(title.baseline) + ylim(y.min,y.max) + ylab(ylab) +
            theme(plot.title=element_text(size=title.size, hjust=0.5)) + theme(legend.position = 'none')
      } else if (plot.cdr3.subtract.baseline==2){ 
        IC.max <- max(unlist(apply(x.norm, 2, function(x){ ind <- which(x!=0); IC <- log(N.aa)/log(2)+sum(x[ind]*log(x[ind])/log(2)); return(IC) })))
        y.max <- min(IC.max*y.inc, log(N.aa)/log(2))
        logo.CDR3.L.rep.max <- ggseqlogoMOD(data=x.norm, additionaAA=additionalAA,  axisTextSizeX = axis.size.max, axisTextSizeY = 8, ylim=c(0, y.max)) +
          labs(title=title.baseline) + ylab(ylab) + theme(plot.title=element_text(size=title.size, hjust=0.5))
      }
    }
  }
  ls <- list(logo.CDR3.L.es, logo.CDR3.L.rep, L.TR, lmax, logo.CDR3.L.es.max, logo.CDR3.L.rep.max)
  names(ls) <- c("ES", "Baseline", "length", "lmax", "ES_max", "Baseline_max")
   
  }  else { 
    ls <- list(list(), list(), c(), 0, ggplot(), ggplot())
    names(ls) <- c("ES", "Baseline", "length", "lmax", "ES_max", "Baseline_max")
    
  }
  
  return(ls)
}


check_input <- function(es.all, chain.list.output="AB",name="input1"){

  if(chain.list.output=="AB"){
    col <- c("TRAV","TRAJ","cdr3_TRA","TRBV","TRBJ","cdr3_TRB")
  }
  if(chain.list.output=="A"){
    col <- c("TRAV","TRAJ","cdr3_TRA")
  }
  if(chain.list.output=="B"){
    col <- c("TRBV","TRBJ","cdr3_TRB")
  }
  
  #Check missing input
  for(cl in col){
    if(cl %in% colnames(es.all) == F){
      cn <- colnames(es.all)
      es.all <- cbind(es.all,"")
      colnames(es.all) <- c(cn, cl)
      print(paste("Missing",cl,"information in",name))
    }
  }
  if("species" %in% colnames(es.all) == F){
    cn <- colnames(es.all)
    es.all <- cbind(es.all,"HomoSapiens")
    colnames(es.all) <- c(cn, "species")
    print(paste("Missing host species in",name,", using HomoSapiens as default"))
  }
  if("model" %in% colnames(es.all) == F){
    cn <- colnames(es.all)
    es.all <- cbind(es.all,"Model_default")
    colnames(es.all) <- c(cn, "model")
    print(paste("Missing model information in ",name,", using the Model_default as default"))
  }
  return(es.all) 

}

clean_input <- function(es.all, use.allele=0, correct.gene.names=1, use.mouse.strain=0, chain.list.output="AB"){
  
  if(chain.list.output=="AB"){
    col <- c("TRAV","TRAJ","cdr3_TRA","TRBV","TRBJ","cdr3_TRB")
    segment.list <- c("TRAV","TRAJ","TRBV","TRBJ")
    cdr3.list <- c("cdr3_TRA","cdr3_TRB")
  }
  if(chain.list.output=="A"){
    col <- c("TRAV","TRAJ","cdr3_TRA")
    segment.list <- c("TRAV","TRAJ")
    cdr3.list <- c("cdr3_TRA")
  }
  if(chain.list.output=="B"){
    col <- c("TRBV","TRBJ","cdr3_TRB")
    segment.list <- c("TRBV","TRBJ")
    cdr3.list <- c("cdr3_TRB")
  }
  
  if(use.allele==1){
    name.list <- gene.allele.list
  } else {
    name.list <- gene.list
  }
  
  #Set to NA CDR3 sequences with incompatible lengths or weird characters
  
  for(cdr3 in cdr3.list){
    ind <- which(nchar(es.all[,cdr3]) < Lmin | nchar(es.all[,cdr3]) > Lmax |
                   grepl('X|x|Z|z|-|_|\\.|\\*', es.all[,cdr3]) == T)
    
    es.all[ind,cdr3] <- NA
  }
  
  #Remove alleles
  if(use.allele==0){
    for(s in segment.list){
      es.all[,s] <- unlist(lapply(es.all[,s], function(x){unlist(strsplit(x,split="*", fixed=T))[1]})) # This is needed if we did not correct gene name
    }
  }
  
  ###################
  # Correct gene names
  # If alleles, it will correct the gene name, and keep the allele. If the allele cannot be found, it will remove it
  # If genes, it will correct the gene name
  # If gene name cannot be corrected, it gives NA
  ###################
  
  if(correct.gene.names==1){
    es.all <- correct.VJnames(es.all, name.list, segment.list)
  }
  if(correct.gene.names==0){
    for(s in segment.list){
      ind <- which(es.all[,s] %in% name.list[[sp]] == F & !is.na(es.all[,s]) )
      if(length(ind)>=1){
        print(c(paste("WARNING: ",s," names in Input TCRs absent from IMGT: ",sep=""), sort(unique(es.all[ind,s])) ))
      }
      es.all[ind,s] <- NA
    }
  }
  
  #Replace empty values by NA
  for(i in col){
    es.all[which(es.all[,i] == ''),i] <- NA
  }
  
  ####
  # Here we should add the check between CDR3 sequences and VJ genes 
  ####
  
  
  ################
  # Do an extra correction for mouse entries, where only gene level analyses are allowed
  # and TRAV genes can be merged
  ################
  
  ind <- which(es.all[,"species"]=="MusMusculus")
  
  if(length(ind)>0){  
  
    if(use.allele==1){
      #Remove the alleles (if(use.allele==0), this was done before)
      for(s in segment.list){
        es.all[ind,s] <- unlist(lapply(es.all[ind,s], function(x){unlist(strsplit(x,split="*", fixed=T))[1]}))
      }
    }
   
    if(use.mouse.strain==0 & chain.list.output != "B"){
      es.all[ind,] <- merge_mouse_TRAV(es.all[ind,])  #WARNING: This only works if alleles have been removed (so far always the case in mouse)
    }
  }
  
  #Remove empty lines
  ind <- apply(es.all,1,function(x){ s <- length(which(is.na(x[col])==F)); return(s)})
  es.all <- es.all[which(ind>0),]
  
  return(es.all)
  
}

merge_mouse_TRAV <- function(es){

  #This has to be run after alleles have been removed and genes have been corrected

  if("TRAV" %in% colnames(es)){
    ind <- which(es[,"species"]=="MusMusculus")
    v.cor <- as.character(unlist(lapply(es[["TRAV"]][ind], function(y){
      if (y %in% names(merge.mouse.TRAV)){
        y <- merge.mouse.TRAV[y]
      }
      return(y)
    })))
    es[ind,"TRAV"] <- v.cor
  }

  return(es)
}




