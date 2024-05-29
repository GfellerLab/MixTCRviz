# Defining functions to handle cases where the mapping does not work.

#' Correcting gene / allele names.
#'
#' Take the gene + allele.
#' If allele, try to correct the gene*allele and return it.
#' If not possible, try to correct only the gene, and return only the gene
#' If no allele, try to correct the gene if needed
#' Cases that could not be corrected are returned as NA
#' @export
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



#' We check the first aa is compatible with the V gene and the last two with the J gene.
#' Cases where one aa is missing due to a different definition of CDR3 are corrected when the first/last 3 aa are compatible with the V/J gene
#' Cases where the V/J genes are not given are only kept if starting with C and ending with F/W
#' Cases of CDR3 incompatible with V/J are put to "".
#' @export
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
          } else { #First amino acid incompatible AND not due to different CDR3 definition => put both the CDR3 to ""
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
