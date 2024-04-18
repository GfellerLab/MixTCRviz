
#This loads the mapping of the gene names, and the function to handle cases where the mapping does not work

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
    gene.allele.list[[tsp]] <- c(gene.allele.list[[tsp]], rownames(read.csv(file=paste(imgt.path,"CDR123/",tsp,"/",s,"_allele.csv", sep=""), row.names = 1)))
    gene.list[[tsp]] <- c(gene.list[[tsp]], rownames(read.csv(file=paste(imgt.path,"CDR123/",tsp,"/",s,".csv", sep=""), row.names = 1)))
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
m <- read.csv(paste(imgt.path,"/CDR123/MusMusculus/TRAV_merge.csv", sep=""))
merge.mouse.TRAV <- m[,2]
names(merge.mouse.TRAV) <- m[,1]

mp.TRAV <- read.csv(paste(imgt.path,"TidyVJ/mapping_TRAV.csv", sep=""), header=F, skip=1)
mp.TRBV <- read.csv(paste(imgt.path,"TidyVJ/mapping_TRBV.csv", sep=""), header=F, skip=1)
mp.TRAJ <- read.csv(paste(imgt.path,"TidyVJ/mapping_TRAJ.csv", sep=""), header=F, skip=1)
mp.TRBJ <- read.csv(paste(imgt.path,"TidyVJ/mapping_TRBJ.csv", sep=""), header=F, skip=1)
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


# Take the gene + allele.
# If allele, try to correct the gene*allele and return it. 
# If not possible, try to correct only the gene, and return only the gene
# If no allele, try to correct the gene if needed
# Cases that could not be corrected are returned as NA
clean.name.allele <- function(g, a, sp){
  
  if(sp != "HomoSapiens" & sp != "MusMusculus"){
    print("Undefined species: ",sp)
  }
  
  if(is.na(g) | g==""){
    ga <- NA
  } else {
    #An allele is given
    if(!is.na(a) & a != ""){
      ga <- paste(g,a,sep="*")
      if(ga %in% gene.allele.list[[sp]] == F){
        #Try correcting the gene name
        if(g %in% gene.list[[sp]] == F ){
          if(!is.na(map[[sp]][g])){
            v <- paste(map[[sp]][g],a,sep="*")
            if(v %in% gene.allele.list[[sp]] == T){
              ga <- v #After correcting the gene, the gene*allele is ok
            } else {
              ga <- map[[sp]][g] #After correcting the gene, the gene*allele is not ok -> give back only the gene
            }
          } else {
            ga <- NA #The gene is wrong and could not be corrected
          }
        } else {
          ga <- g #The gene was ok, but the allele was problematic -> give back only the gene
        }
      }
    } else {  #Allele is not given
      
      if(g %in% gene.list[[sp]] == F ){
        if(!is.na(map[[sp]][g])){
          ga <- map[[sp]][g] #The gene was wrong, but can be corrected
        } else {
          ga <- NA #The gene was wrong and could not be corrected
        }
      } else {
        ga <- g
      }
    }
  }
  return(ga)  
}



#We check the first aa is compatible with the V gene and the last two with the J gene.
#Cases where one aa is missing due to a different definition of CDR3 are corrected when the first/last 3 aa are compatible with the V/J gene
#Cases where the V/J genes are not given are only kept if starting with C and ending with F/W
#Cases of CDR3 incompatible with V/J are put to "".
clean.cdr3 <- function(v,j,cdr3,sp,chain){
  
  cdr3.cor <- cdr3
  v.cor <- v
  j.cor <- j
  
  if(cdr3 != ""){
    
    ######
    # Check with the beginning versus the V gene
    ######
    first <- str_sub(cdr3,1,1)
    if(v != "" & !is.na(cdr123[[sp]][[chain]][v,3])){
      ref.first <- str_sub(cdr123[[sp]][[chain]][v,3],1,1)
      if(ref.first != ""){  #This is needed because some V alleles have no CDR3 (rare alleles)
        if(first != ref.first){
          first2 <- str_sub(cdr3,1,3)
          ref.next <- str_sub(cdr123[[sp]][[chain]][v,3],2,4)
          if(first2 == ref.next){ # Missing first amino acid
            cdr3.cor <- paste(ref.first,cdr3,sep="")
            #print(c(cdr3,cdr3.cor,v, sp))
          } else { #First amino acid incompatible AND not due to different CDR3 definition => put both the CDR3 and the V gene to ""
            cdr3.cor <- ""
            v.cor <- ""
          }
        } 
      } else {
        if(first != "C"){
          cdr3.cor <- ""
        }
      }
    } else {
      if(first != "C"){
        cdr3.cor <- ""
      }
      v.cor <- ""
    }
    
    ######
    # Check with the end versus the J gene
    ######
    if(j != "" & !is.na(Jseq[[sp]][[chain]][j,1])){
      last2 <- str_sub(cdr3,-2,-1)
      ref.last2 <- str_sub(Jseq[[sp]][[chain]][j,1],-2,-1)
      
      #print(c(last2,ref.last2, sp, j))
      if(last2 != ref.last2){ #Check if the last two aa correspond to the J gene
        ref.next <- str_sub(Jseq[[sp]][[chain]][j,1],-4,-2)
        last3 <- str_sub(cdr3,-3,-1)
        if(last3==ref.next){ #If this is due to a missing residue, add the missing aa
          cdr3.cor <- paste(cdr3,str_sub(Jseq[[sp]][[chain]][j,1],4,4), sep="")
          #print(c(cdr3,cdr3.cor,j, sp))
        } else {
          cdr3.cor <- ""
          j.cor <- ""
        }
      }
    } else { #J gene is not given or not valid (mixing of chains)
      last <- str_sub(cdr3,-1,-1)
      if(last != "F" & last != "W"){
        cdr3.cor <- ""
      }
      j.cor <- ""
    }
    
  } 
  return(c(v.cor, j.cor, cdr3.cor))
}
