# Description ---------------------------------------------------------------
# In this file we create and export some variables used by the package

if (FALSE){
  # Don't run this code when building the package.

# Mapping of the gene names -----------------------------------------------
  segment.list <- c("TRAV", "TRBV", "TRAJ", "TRBJ")
  species.list <- c("HomoSapiens", "MusMusculus")

  # Take the list of alleles
  # This is only done for human
  # We need to do it for mouse
  # Alternatively, we could take the fasta files.


  gene.allele.list <- list()
  gene.list <- list()
  allele.default <- list()
  for(tsp in species.list){
    gene.allele.list[[tsp]] <- c()
    gene.list[[tsp]] <- c()
    for(s in segment.list){
      gene.allele.list[[tsp]] <- c(gene.allele.list[[tsp]], rownames(read.csv(file=paste("data_raw/CDR123/",tsp,"/",s,"_allele.csv", sep=""), row.names = 1)))
      gene.list[[tsp]] <- c(gene.list[[tsp]], rownames(read.csv(file=paste("data_raw/CDR123/",tsp,"/",s,".csv", sep=""), row.names = 1)))
    }

    # Build the mapping to the most likely allele (needed when allele names are not given, and we infer them)
    allele.default[[tsp]] <- rep("01", length(gene.list[[tsp]]))
    names(allele.default[[tsp]]) <- gene.list[[tsp]]
    # Do some manual correction
    if(tsp=="HomoSapiens"){
      allele.default[[tsp]]["TRAV14/DV4"] <- "02"
      allele.default[[tsp]]["TRAV36/DV7"] <- "02"
    }
  }

  #For mouse TRAV, take the merging between different strains
  m <- read.csv(paste("data_raw/CDR123/MusMusculus/TRAV_merge.csv", sep=""))
  merge.mouse.TRAV <- m[,2]
  names(merge.mouse.TRAV) <- m[,1]

  mp.TRAV <- read.csv("data_raw/TidyVJ/mapping_TRAV.csv", header=F, skip=1)
  mp.TRBV <- read.csv("data_raw/TidyVJ/mapping_TRBV.csv", header=F, skip=1)
  mp.TRAJ <- read.csv("data_raw/TidyVJ/mapping_TRAJ.csv", header=F, skip=1)
  mp.TRBJ <- read.csv("data_raw/TidyVJ/mapping_TRBJ.csv", header=F, skip=1)
  mp <- rbind(mp.TRAV,mp.TRBV,mp.TRAJ,mp.TRBJ)
  # Keep only the entries that could be mapped
  mp <- mp[mp[,5]=="yes",]

  #Build the mapping separately for human and mouse
  map <- list()
  pos.h <- which(mp[,4]=="Human")
  map[["HomoSapiens"]] <- mp[pos.h,2]
  names(map[["HomoSapiens"]]) <- mp[pos.h,1]
  pos.m <- which(mp[,4]=="Mouse")
  map[["MusMusculus"]] <- mp[pos.m,2]
  names(map[["MusMusculus"]]) <- mp[pos.m,1]

  # usethis::use_data(gene.allele.list, gene.list, allele.default,
  #   merge.mouse.TRAV, map, overwrite=T, internal=F)
  # # Uncomment if you want to have some of these variable available outside
  # # of MixTCRviz (used through MixTCRviz::gene.allele.list for example or
  # # when library("MixTCRviz") is used).



# Defining other parameters -----------------------------------------------
  # This initiates a few values and load the mapping of CDR1/2 sequences from
  # the gene*allele names.
  library(ggplot2)

  Lmin <- 7
  Lmax <- 22

  aa <- "ACDEFGHIKLMNPQRSTVWY"
  aa.list <- unlist(strsplit(aa, split=""))
  N.aa <- length(aa.list)
  #map.aa <- 1:20
  #names(map.aa) <- aa.list

  gap <- c("-", "X", ".", "*","_")  # "x" is kept for 'real' gaps in CDR1/2, while all the other symbols are for missing data.
  th <- theme(plot.title = element_text(size = 8, hjust=0.5),
    axis.title=element_text(size=8))
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
      tcdr123[[chain]] <- read.csv(file=paste("data_raw/CDR123/",sp,"/",chain,"V_allele.csv", sep=""), row.names = 1)
      tJseq[[chain]] <- read.csv(file=paste("data_raw/CDR123/",sp,"/",chain,"J_allele.csv", sep=""), row.names = 1)

      #Take the data without allele information
      tcdr123[[chain]] <- rbind(tcdr123[[chain]],read.csv(file=paste("data_raw/CDR123/",sp,"/",chain,"V.csv", sep=""), row.names = 1))
      tJseq[[chain]] <- rbind(tJseq[[chain]],read.csv(file=paste("data_raw/CDR123/",sp,"/",chain,"J.csv", sep=""), row.names = 1))

      #Change the gap ("-") int "x"
      tcdr123[[chain]][] <- data.frame(lapply(tcdr123[[chain]], function(x){gsub(pattern="-", replacement = "g", x)}))
      tJseq[[chain]][] <- data.frame(lapply(tJseq[[chain]], function(x){gsub(pattern="-", replacement = "g", x)}))

    }

    cdr123[[sp]] <- list(tcdr123[["TRA"]],tcdr123[["TRB"]])
    names(cdr123[[sp]]) <- c("TRA","TRB")
    Jseq[[sp]] <- list(tJseq[["TRA"]],tJseq[["TRB"]])
    names(Jseq[[sp]]) <- c("TRA","TRB")

  }

  # usethis::use_data(cdr123, Jseq, th, yl, aa, aa.list, N.aa, chain.small, gap,
  #   Lmin, Lmax, species.list, overwrite=T, internal=F)


# Saving all these variable for internal use within the package ------------------
# functions

  usethis::use_data(gene.allele.list, gene.list, allele.default,
    merge.mouse.TRAV, map, cdr123, Jseq, th, yl, aa, aa.list, N.aa,
    chain.small, gap, Lmin, Lmax, species.list, overwrite=F, internal=T)



  # EOF ---------------------------------------------------------------------
}

