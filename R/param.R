library(ggplot2)

#This initiates a few values and load the mapping of CDR1/2 sequences from the gene*allele names.

Lmin <- 7
Lmax <- 22

aa <- "ACDEFGHIKLMNPQRSTVWY"
aa.list <- unlist(strsplit(aa, split=""))
N.aa <- length(aa.list)
#map.aa <- 1:20
#names(map.aa) <- aa.list

gap <- c("-", "X", ".", "*","_")  # "x" is kept for 'real' gaps in CDR1/2, while all the other symbols are for missing data.
th <- theme(plot.title = element_text(size = 8, hjust=0.5), axis.title=element_text(size=8))
yl <- ylim(0,2)

chain.list <- c("TRA", "TRB")
chain.small <- c("a","b")
names(chain.small) <- chain.list

species.list <- c("HomoSapiens", "MusMusculus")

### Load the aligned CDR1 and CDR2 sequences

cdr123 <- list()
Jseq <- list()

for(sp in species.list){

  tcdr123 <- list()
  tJseq <- list()
  
  for(chain in chain.list){
    tcdr123[[chain]] <- read.csv(file=paste(imgt.path,"CDR123/",sp,"/",chain,"V_allele.csv", sep=""), row.names = 1)
    tJseq[[chain]] <- read.csv(file=paste(imgt.path,"CDR123/",sp,"/",chain,"J_allele.csv", sep=""), row.names = 1)
    
    #Take the data without allele information
    tcdr123[[chain]] <- rbind(tcdr123[[chain]],read.csv(file=paste(imgt.path,"CDR123/",sp,"/",chain,"V.csv", sep=""), row.names = 1))
    tJseq[[chain]] <- rbind(tJseq[[chain]],read.csv(file=paste(imgt.path,"CDR123/",sp,"/",chain,"J.csv", sep=""), row.names = 1))
    
    #Change the gap ("-") int "x"
    tcdr123[[chain]][] <- data.frame(lapply(tcdr123[[chain]], function(x){gsub(pattern="-", replacement = "g", x)}))
    tJseq[[chain]][] <- data.frame(lapply(tJseq[[chain]], function(x){gsub(pattern="-", replacement = "g", x)}))
    
  }

  cdr123[[sp]] <- list(tcdr123[["TRA"]],tcdr123[["TRB"]])
  names(cdr123[[sp]]) <- c("TRA","TRB")
  Jseq[[sp]] <- list(tJseq[["TRA"]],tJseq[["TRB"]])
  names(Jseq[[sp]]) <- c("TRA","TRB")
  
}

