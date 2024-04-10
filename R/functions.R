
#####
# Define some function
#####

build_stat <- function(es){
  
  
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
  
  for(chain in chain.list){
    
    seg <- c(paste(chain,"V", sep=""),paste(chain,"J", sep=""))
    cdr3 <- paste("cdr3_",chain,sep="")
    
    countV[[chain]] <- table(es[,seg[1]])
    countJ[[chain]] <- table(es[,seg[2]])

    countL[[chain]] <- table(nchar(es[,cdr3]))
    
    L[[chain]] <- as.numeric(names(countL[[chain]]))

    countVJ[[chain]] <- table(es[,seg[1]], es[,seg[2]])
        
    countV.L[[chain]] <- list()
    countJ.L[[chain]] <- list()
    countVJ.L[[chain]] <- list()
    for(lg in L[[chain]]){
      ind <- which(nchar(es[,cdr3])==lg)
      countV.L[[chain]][[lg]] <- table(es[ind,paste(chain,"V",sep="")])
      countJ.L[[chain]][[lg]] <- table(es[ind,paste(chain,"J",sep="")])
      countVJ.L[[chain]][[lg]] <- table(es[ind,paste(chain,"V",sep="")],es[ind,paste(chain,"J",sep="")])
    }
    
    if(length(countV[[chain]])>0){
      countCDR1[[chain]] <- count_aa(cdr123[[sp]][[chain]][es[,seg[1]],"CDR1"], keep.gap=1)
      countCDR2[[chain]] <- count_aa(cdr123[[sp]][[chain]][es[,seg[1]],"CDR2"], keep.gap=1)
    }
    countCDR3.L[[chain]] <- count_aa(es[,cdr3], keep.gap=0)
    
    
  }
  count <- list(L, countL, countV, countJ, countV.L, countJ.L, countCDR1, countCDR2, countCDR3.L, countVJ, countVJ.L)
  names(count) <- c("L", "countL", "countV", "countJ", "countV.L", "countJ.L", "countCDR1", "countCDR2", "countCDR3.L", "countVJ", "countVJ.L")
  
  return(count)
}

#Compute the counts of each aa at each position.
#Would be better to have the option of treating gaps in CDR1/2 as separate amino acids, while missing data can be treated as 'unspecific'
count_aa <- function(x, keep.gap=0){
  
  if(keep.gap == 0){    #All gaps, including "g" are treated as unspecific data. This can be useful for visualisation.
    tgap <- c(gap,"g")
    taa.list <- aa.list
  } else {   # "x" are treated as additonal aa, other 'gaps' are discarded. This is more correct for modelling CDR1/CDR2 loops.
    tgap <- c()
    taa.list <- c(aa.list,"g")
  }
  
  #First get the list of length
  lc <- nchar(x)
  lg <- sort(unique(lc))
  lg <- lg[lg>0]  # This ensures that the cases of empty sequences are never counted.
  m.list <- list()
  for(l in lg){
    ind <- which(lc==l)
    tx <- x[ind]
    m <- matrix(0, nrow=length(taa.list), ncol=l)
    rownames(m) <- taa.list
    for(p in 1:l){
      s <- substr(tx,p,p)
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
    m.list[[l]] <- m
  }
  return(m.list)
}

#This assumes the matrix includes gaps (21 rows)
build_cdr12_motif <- function(x, keep.gap=0){
  
  if(keep.gap==0){
    g <- x["g",]
    for(a in aa.list){
      x[a,] <- x[a,]+g/length(aa.list)
    }
    x <- x[aa.list,]
  }
  x <- scale(x, center=F, scale=colSums(x))
  return(x)
}

#Compute the weighted average of VJ usage for different CDR3 lengths
#This is no longer useful.
weighted_VJcount <- function(x,y){
  
  #Build the list of all segment names
  nm <- c()
  for(i in names(y)){
    nm <- c(nm, names(x[[as.numeric(i)]]) )
  }
  nm <- unique(nm)
  
  ct <- rep(0, length(nm))
  names(ct) <- nm
  for(i in names(y)){
    v <- x[[as.numeric(i)]]
    v <- v/sum(v)  #This is the normalized distribution of V or J for a given length
    
    #print(v[1:10])
    
    for(j in 1:length(v)){
      ct[names(v)[j]] <- ct[names(v)[j]] + v[j]*y[i]
    }
  }
  return(ct)
}

# Compute the weighted average of motifs for different V/J usage
# This is important since the choice of V/J has huge influence on the CDR3 motif.
# By considering V/J usage in baseline, it helps identifying what is epitope-specific
weighted_countCDR3 <- function(x,y){
  
  countCDR3 <- list()
  L <- which(y != "NULL")
  for(lg in L){
    countCDR3[[lg]] <- matrix(0, nrow=N.aa, ncol=lg)
    rownames(countCDR3[[lg]]) <- aa.list
    if(sum(y[[lg]])>0){
      w <- y[[lg]]/sum(y[[lg]])
      for(v in rownames(y[[lg]])){
        for(j in colnames(y[[lg]])){
          if(w[v,j]>0){
            s <- paste(v,j,sep="_")
            if(length(x[[s]])>=lg){
              if(length(x[[s]][[lg]]) > 0){
                countCDR3[[lg]] <- countCDR3[[lg]]+scale(x[[s]][[lg]], center=F, scale=colSums(x[[s]][[lg]]))*w[v,j]
              }
            }
          }
        }
      }
    }
    if(colSums(countCDR3[[lg]])[1]!=0){
      countCDR3[[lg]] <- scale(countCDR3[[lg]],center=F, scale=colSums(countCDR3[[lg]]))
    } else {
      countCDR3[[lg]][1:N.aa,1:lg] <- 0.05 #This is the case if the were no data for all V-J for this lenght in the repertoire
      
    }
  }
  return(countCDR3)
}

# Compute the weighted average of CDR3 length for the observed V/J usage
# This is important since CDR3 length is primarily determined by the length of the V and J segments
weighted_countL <- function(x,y){
  
  countL <- rep(0,times=Lmax-Lmin+1)
  names(countL) <- as.character(Lmin:Lmax)
  
  y <- y/sum(y)
  for(v in rownames(y)){
    for(j in colnames(y)){
      if(y[v,j]>0){
        s <- paste(v,j,sep="_")
        if(length(x[[s]])>0){
          x[[s]] <- x[[s]]/sum(x[[s]])*y[v,j]
          for(l in names(x[[s]])){
            countL[l] <- countL[l]+x[[s]][l]
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
plotVJ <- function(count.es, count.rep, info){
  
  #count.es is on the Y axis, count.rep on the X
  gene <- info[1]
  
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
  
  lim.y <- max(count[,c(type1)] )*1.2
  lim.x <- max(count[,c(type2)] )*1.2
  
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
  
  if(gene=="TRAV" | gene=="TRBV"){
    ylab <- type1
  } else{ylab=""}
  
  #Plot the comparison between epitope specific and repertoires
  count.plot <- list()
  if(length(count.es)>0){
    count.plot <- ggplot(count.df, aes(x=X, y=Y, label=label)) + 
      geom_point() + geom_abline(col="orange",linetype="dashed",linewidth=1) + 
      ggtitle(gene) + xlim(0, lim.x) + ylim(0,lim.y) + 
      theme(plot.title = element_text(size = 15, hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=15)) +
      geom_label_repel(size = 3, nudge_y=0.02, box.padding = 0.15) +
      xlab(type2) + ylab(ylab)
  } else {count.plot <- ggplot()}
  
  return(count.plot)
  
}


correct.VJnames <- function(es.all, name.list){
  
  segment.list <- c("TRAV", "TRAJ", "TRBV", "TRBJ")
  
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


plotLD <- function(lc.es,lc.rep,info){
  
  l.all <- sort(as.numeric(unique(c(names(lc.es), names(lc.rep)))))
  
  ct <- 1
  ld.es <- c()
  ld.rep <- c()
  
  for(l in as.character(l.all)){
    if(!is.na(lc.es[l])){ld.es[ct] <- lc.es[l]} else {ld.es[ct] <- 0}
    if(!is.na(lc.rep[l])){ld.rep[ct] <- lc.rep[l]} else {ld.rep[ct] <- 0}
    ct <- ct+1
  }
  ld.es <- ld.es/sum(ld.es)
  ld.rep <- ld.rep/sum(ld.rep)
  
  ##########
  #Plot the comparison for length distribution
  ##########
  
  v1 <- c(l.all,l.all); 
  v2 <- c(ld.es, ld.rep);
  v3 <- c( rep(info[2], length(l.all)), rep(info[3], length(l.all))) ; 
  ld.df <- data.frame(v1,v2,v3)
 # ld.df$v3 <- factor(ld.df$v3, levels=c("Epitope specific", "Baseline"))
  ld.df$v3 <- factor(ld.df$v3, levels=c(info[2], info[3]))
  
  ld.plot <-  ggplot(ld.df) + geom_point(aes(x=v1, y=v2, color=v3)) + ggtitle(paste("N =",sum(lc.es))) +
    geom_line(aes(x=v1, y=v2, color=v3)) +
    theme(legend.key.size = unit(0.2, 'cm'), legend.position="top", legend.title=element_blank(),  legend.text=element_text(size=12)) +
    xlab(paste("Length_CDR3",info[1],sep="")) + ylab("Distribution") + 
    guides(color = guide_legend(nrow = 2)) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=15), plot.title = element_text(size=15,hjust = 0.5))
  
  return(ld.plot)
  
  
}


plotCDR3 <- function(lc.es, lc.rep, countCDR3.es, countCDR3.rep, info, comp.baseline){
  
  L.TR <- as.numeric(intersect(names(lc.es), names(lc.rep)))
  
  pwm.rep <- list()
  pwm.es <- list()
  
  logo.CDR3.L.es <- list()
  logo.CDR3.L.rep <- list()
  
  tl <- lc.es[as.character(L.TR)]
  lmax.es <- as.numeric(names(tl[which.max(tl)]))
  lmax <- lmax.es
  
  
  for(l in L.TR){
    
    lc <- as.character(l)
    
    pwm.es[[l]] <- scale(countCDR3.es[[l]], center=F, scale=colSums(countCDR3.es[[l]]))
    pwm.rep[[l]] <- scale(countCDR3.rep[[l]], center=F, scale=colSums(countCDR3.rep[[l]]))
    
    title <- paste("CDR3", info[1]," ",info[2],"\n(",lc.es[[lc]],")", sep="")
    ylab <- ""
    
    logo.CDR3.L.es[[l]] <- ggseqlogoMOD(data=pwm.es[[l]], additionaAA=additionalAA,  axisTextSizeX = 10, axisTextSizeY = 10) + 
      labs(title=title) + ylab(ylab) + theme(plot.title=element_text(size=12, hjust=0.5))
 
    
    title <- paste("CDR3", info[1]," ",info[3], sep="")
    ylab <- ""
    if(comp.baseline==0){title <- paste(title, "\n(",lc.rep[[lc]],")", sep="")}
    
    logo.CDR3.L.rep[[l]] <- ggseqlogoMOD(data=pwm.rep[[l]], additionaAA=additionalAA,  axisTextSizeX = 10, axisTextSizeY = 10) +
      labs(title=title) + ylab(ylab) + theme(plot.title=element_text(size=12, hjust=0.5))
    
  }
  
  ls <- list(logo.CDR3.L.es, logo.CDR3.L.rep, L.TR, lmax)
  names(ls) <- c("ES", "Baseline", "length", "lmax")
  return(ls)
}

