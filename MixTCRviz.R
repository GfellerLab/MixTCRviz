
if (!require("ggseqlogo")){
  devtools::install_github("GfellerLab/ggseqlogo")
}

library(ggplot2)
library(ggpubr)
library(ggseqlogo)  ## WARNING: This is the modified version of ggseqlogo, with the small aa. Load if from our github
library(seqinr)
library(ggrepel)
library(mgcv)

### CURRENTLY, THIS HAS TO BE GIVEN BY THE USER ###
MixTCRviz.path <- "/Users/davidgfeller/Research/LICR/Research/TCR/Benchmarking/MixTCRviz/"


imgt.path <- paste(MixTCRviz.path,"/data/", sep="") #This is the path to the reference data. Always declare it
source(paste(MixTCRviz.path,"/R/mapping_function.R", sep=""))
source(paste(MixTCRviz.path,"/R/param.R", sep=""))
source(paste(MixTCRviz.path,"/R/functions.R", sep=""))

MixTCRviz <- function(input1, output.path, 
                      input2="", baseline.file="",
                      use.allele=0, correct.gene.names=1, use.mouse.strain=0,
                      renormVJ=1, N.min=10,
                      plot=1, plot.cdr12.motif=0, plot.oneline=0, plot.logo.length=0,
                      chain.list.output="AB", es.name="Epitope specific"){
  
  # Mandatory Parameters:
  # - input1: Required: Input TCRs. Can be a filename or a data.frame
  # - output.path: Required: path for output
  
  # Optional parameters
  # - input2: Reference TCRs. Can be a filename or a data.frame
  # - correct.gene.names <- 1 # Attempt to correct V/J names not in IMGT
  # - plot <- 0  # Build the plots.
  # - plot.cdr12.motif <- 0 #Build the plots for the CDR12 motifs 
  # - plot.logo.length <- 1 #Plot the logos of different lengths
  # - plot.oneline <- 0
  # - N.min <- 10
  # - renormVJ <- 1
  # - use.allele <- 0 #0: The allele are removed or absent in the input file; 1: the allele are present or added (not recommanded)
  # - use.mouse.strain <- 0 #0: Merge the TRAV segments in mouse corresponding to different strains (TRAV10/TRAV10D/TRAAV10N)
  # - chain.list.output="AB": Decide which chains to plot 
  # - es.name="Epitope specific": Decide the name for the input TCRs (e.g., "Epitope specific") 
  
  #######
  # Choose some parameters
  #######
  
  keep.gap.pwm <- 0 # 1: means that 'g' (gaps in CDR1/2) are treated as an additional aa in the logos. 0: means that 'g' are treated as unspecific (i.e., 0.05) in the logos
  if(keep.gap.pwm==1){taa.list <- c(aa.list,"g"); additionalAA <- "g"} else {taa.list <- aa.list; additionalAA <- ""}
  
  min.logo <- 10 #Minimum number of sequences to plot the logos (when plotting different lengths)
  
  col.TCR <- c("TRAV", "TRAJ", "cdr3_TRA", "TRBV", "TRBJ", "cdr3_TRB")
  
  if(plot.cdr12.motif==1){plot.oneline=0}
  if(plot.oneline==1){
    th <- theme(plot.title = element_text(size = 8, hjust=0.5), axis.title=element_text(size=4))
  }
  
  if(N.min < 1){N.min <- 1}
  
  if(!file.exists(output.path)){
    dir.create(output.path, recursive = TRUE);
  }
  
  if(is.character(input2)==T){
    if(input2 == ""){
      comp.baseline <- 1
      baseline.name <- "Baseline"
    } else {
      comp.baseline <- 0
      baseline.name <- "Reference"
    }
  } else if(is.data.frame(input2)==T){
    comp.baseline <- 0
    baseline.name <- "Reference"
  }
  
  if(use.allele==1){
    name.list <- gene.allele.list
  } else {
    name.list <- gene.list
  }
  
  ############################
  # Load all the input data
  ############################
  
  if(is.character(input1)==T){
    es.all <- read.csv(input1)
  } else if(is.data.frame(input1)==T){
    es.all <- input1
  }
  
  #Check if some columns are missing
  es.all <- check_input(es.all)
  
  #Remove alleles
  if(use.allele==0){
    for(s in segment.list){
      es.all[,s] <- unlist(lapply(es.all[,s], function(x){unlist(strsplit(x,split="*", fixed=T))[1]})) # This is needed if we did not correct gene name
    }
  }  
  
  ###################
  # Correct gene names
  # If use.allele==1, it will correct the gene name, and keep the allele. If the allele cannot be found, it will remove it
  # If use.allele==0, it will correct the gene name
  # If gene name cannot be corrected, it gives NA
  ###################
  
  if(correct.gene.names==1){
    es.all <- correct.VJnames(es.all, name.list)
  }
  
  #Replace empty values by NA
  for(i in col.TCR){
    es.all[which(es.all[,i] == ''),i] <- NA
  }
  
  
  ########################
  #Select the model which should be considered
  ########################
  
  #Take samples where there is at least one chain with enough data
  
  md <- unique(es.all[,"model"])
  st <- lapply(md, function(x){ 
    i <- which(es.all[,"model"]==x);
    nA <- length(which(!is.na(es.all[i,"TRAV"]) & !is.na(es.all[i,"TRAJ"]) & !is.na(es.all[i,"cdr3_TRA"]) )) 
    nB <- length(which(!is.na(es.all[i,"TRBV"]) & !is.na(es.all[i,"TRBJ"]) & !is.na(es.all[i,"cdr3_TRB"]) )) 
    if(nA>=N.min | nB>=N.min){
      return(1)
    } else {
      return(0)
      print(paste("Model",x,"will not be considered (less than N.min data)"))
    }
  })
  model.list <- md[st==1]
  
  if(length(model.list)==0){
    stop("No model has enough data to use MixTCRviz. Check your input or lower the N.min parameter (default 10)")
  }
  
  #######################################
  # Analyse the input
  #######################################
  
  for(model in model.list){
    
    use.allele.es <- use.allele
    
    print (model)
    
    pg.all <- list()
    pg.length <- list()
    tl.logo <- list() #Keep track of the number of logos of different length
    
    ######
    # Select the input models with enough data
    ######
    
    es <- es.all[which(es.all[,"model"]==model),]
    sp <- unique(es[1,"species"])
    if(length(sp)>1){
      stop("Multiple species provided for the same model. Use a different model name for each species")
    } else {
      if(sp %in% species.list==F){
        stop("Unknown species. MixTCRviz supports only \"HomoSapiens\" or \"MusMusculus\"")
      }
    }
    if(sp=="MusMusculus" & use.allele.es==1){
      print("Alleles currently not supported in mouse. The data will be treated at the gene level")
      use.allele.es <- 0
      for(s in segment.list){
        es[,s] <- unlist(lapply(es[,s], function(x){unlist(strsplit(x,split="*", fixed=T))[1]}))
      }  
    }
    #This is to handle cases where use.allele was forced to 0 (mouse)
    if(use.allele.es==1){
      name.list <- gene.allele.list
    } else {
      name.list <- gene.list
    }
    
    #Merge the TRAV genes in mouse across different strains
    if(sp=="MusMusculus" & use.mouse.strain==0 & use.allele.es==0){
      es <- merge_mouse_TRAV(es)  #WARNING: This only works if use.allele.es==0 (so far always the case in mouse)
    }
    
    #Do a check of the gene names missing from IMGT
    if(correct.gene.names==0){
      for(s in segment.list){
        ind <- which(es[,s] %in% name.list[[sp]] == F & !is.na(es[,s]) )
        if(length(ind)>=1){
          print(c(paste("WARNING: ",s," names in Input TCRs absent from IMGT: ",sep=""), sort(unique(es[ind,s])) ))
        }
        es[ind,s] <- NA
      }
    }
    es <- es[,col.TCR]
    for(i in col.TCR){
      es[which(es[,i] == ''),i] <- NA
    }
    
    count <- build_stat(es)
    L.es <- count$L
    countL.es <- count$countL
    countV.es <- count$countV
    countJ.es  <- count$countJ
    countV.L.es <- count$countV.L
    countJ.L.es  <- count$countJ.L
    countVJ.es <- count$countVJ
    countVJ.L.es <- count$countVJ.L
    countCDR1.es  <- count$countCDR1
    countCDR2.es <- count$countCDR2
    countCDR3.L.es <- count$countCDR3.L
    
    summary <- list(L.es, countL.es, countV.es, countJ.es, countV.L.es, countJ.L.es, countCDR1.es, countCDR2.es, countCDR3.L.es, countVJ.es, countVJ.L.es)
    names(summary) <- c("L", "countL", "countV", "countJ", "countV.L", "countJ.L", "countCDR1", "countCDR2", "countCDR3.L", "countVJ", "countVJ.L")
    dir <- paste(output.path,"/stats/", sep="")
    if(!file.exists(dir)){
      dir.create(dir);
    }
    saveRDS(summary, file=paste(output.path,"/stats/",model,".rds", sep=""))
    
    
    if(plot==1){
      
      if(comp.baseline==1){
        
        #####
        # Load the baseline repertoire corresponding to the species
        #####
        
        if(baseline.file==""){
          if(sp=="HomoSapiens"){
            if(use.allele.es==1){
              all.baseline <- readRDS(file=paste(MixTCRviz.path,"/Rdata/summary_",sp,".rds", sep=""))
            } else {
              all.baseline <- readRDS(file=paste(MixTCRviz.path,"/Rdata/summary_",sp,"_noallele.rds", sep=""))
            }
          }
          if(sp=="MusMusculus"){
            if(use.allele.es==0){  #N.B. Currently always the case in mouse
              if(use.mouse.strain==1){
                all.baseline <- readRDS(file=paste(MixTCRviz.path,"/Rdata/summary_",sp,"_noallele_SEQTR.rds", sep=""))
              } else {
                all.baseline <- readRDS(file=paste(MixTCRviz.path,"/Rdata/summary_",sp,"_noallele_noStrain_SEQTR.rds", sep=""))
              }
            }
          }
        } else {
          all.baseline <- readRDS(file=baseline.file)
        }
        
        L.baseline <- all.baseline$L
        countL.baseline <- all.baseline$countL
        countV.baseline <- all.baseline$countV
        countJ.baseline  <- all.baseline$countJ
        countV.L.baseline <- all.baseline$countV.L
        countJ.L.baseline  <- all.baseline$countJ.L
        countL.VJ.baseline  <- all.baseline$countL.VJ
        countVJ.baseline  <- all.baseline$countVJ
        countCDR1.baseline  <- all.baseline$countCDR1
        countCDR2.baseline <- all.baseline$countCDR2
        countCDR3.L.baseline <- all.baseline$countCDR3.L
        countCDR3.VJL.baseline <- all.baseline$countCDR3.VJL
        
        
      } else {
        
        #####
        # Load the other repertoire to compare with
        # Still need to implement a few check that the data are available.
        #####
        
        if(is.character(input2)==T){
          es2.all <- read.csv(input2)
        } else if (is.data.frame(input2)==T){
          es2.all <- input2
        }
        es2.all <- check_input(es2.all)
        
        for(i in col.TCR){
          es2.all[which(es2.all[,i] == ''),i] <- NA
        }
        
        
        ind <- which(es2.all[,"model"]==model)
        if(length(ind)==0){
          stop(c("Missing model in Reference: ", model))
        } else {
          es2 <- es2.all[ind,]
        }    
        
        #Remove alleles
        if(use.allele.es==0){
          for(s in segment.list){
            es2[,s] <- unlist(lapply(es2[,s], function(x){unlist(strsplit(x,split="*", fixed=T))[1]}))
          }
        }
        if(correct.gene.names==1){
          es2 <- correct.VJnames(es2, name.list)
        }
        
        if(sp=="MusMusculus" & use.mouse.strain==0 & use.allele.es==0){
          es2 <- merge_mouse_TRAV(es2)
        }
        
        for(s in segment.list){
          #Put to NA the V/J segments for which we do not have sequences (this should not be the case if the repertoire is good.
          ind <- which(es2[,s] %in% name.list[[sp]] == F & !is.na(es2[,s]))
          if(length(ind)>=1){
            print(c(paste("WARNING: ",s," names in Reference absent from IMGT: ",sep=""), sort(unique(es2[ind,s])) ))
          }
          es2[ind,s] <- NA
        }
        
        count <- build_stat(es2)
        L.baseline <- count$L
        countL.baseline <- count$countL
        countV.baseline <- count$countV
        countJ.baseline  <- count$countJ
        countV.L.baseline <- count$countV.L
        countJ.L.baseline  <- count$countJ.L
        countCDR1.baseline  <- count$countCDR1
        countCDR2.baseline <- count$countCDR2
        countCDR3.L.baseline <- count$countCDR3.L
      }
      
      for(chain in chain.list){
        
        if(comp.baseline==1){
          #Check segments that were in the ES, but not in baseline
          miss.V.baseline <- setdiff(names(countV.es[[chain]]), names(countV.baseline[[chain]]))
          miss.J.baseline <- setdiff(names(countJ.es[[chain]]), names(countJ.baseline[[chain]]))
          if(length(miss.V.baseline)>=1){print(c(paste("WARNING: ",chain,"V in Input TCRs, but absent from baseline: ", sep=""), miss.V.baseline))}
          if(length(miss.J.baseline)>=1){print(c(paste("WARNING: ",chain,"J in Input TCRs, but absent from baseline: ", sep=""), miss.J.baseline))}
        }
        
        #Make sure there are CDR3 sequences
        if(length(countL.es[[chain]])>0){
          
          if(length(countV.es[[chain]])==0){
            print(paste("WARNING: No ", chain,"V segment", sep=""))
          }
          if(length(countJ.es[[chain]])==0){
            print(paste("WARNING: No ",chain,"J segment", sep=""))
          }
          
          
          
          info <- c(chain.small[chain], es.name, baseline.name)
          
          if(comp.baseline==1){
            #ld.plot <- plotLD.weighted(countL.es, countL.baseline, info)
            if(length(countV.es[[chain]])>0 & length(countJ.es[[chain]])>0 & renormVJ==1){
              info[3] <- paste(info[3], "VJ",sep=" | ")
              ld.plot <- plotLD(countL.es[[chain]], weighted_countL(countL.VJ.baseline[[chain]], countVJ.es[[chain]]), info)
            } else{
              ld.plot <- plotLD(countL.es[[chain]], countL.baseline[[chain]], info)
            }
          } else {
            ld.plot <- plotLD(countL.es[[chain]], countL.baseline[[chain]], info)
          }
          
          #######
          # Plot comparison of V/J usage
          #######
          
          infoV <- c(paste(chain,"V", sep=""), es.name, baseline.name)
          infoJ <- c(paste(chain,"J", sep=""), es.name, baseline.name)
          
          if(length(countV.es[[chain]])>0){
            countV.plot <- plotVJ(countV.es[[chain]], countV.baseline[[chain]], infoV)
          } else {
            countV.plot <- ggplot()
          }
          if(length(countJ.es[[chain]])>0){
            countJ.plot <- plotVJ(countJ.es[[chain]], countJ.baseline[[chain]], infoJ)
          } else {
            countJ.plot <- ggplot()
          }
          #######
          # Plot comparison of motifs for CDR1 and CDR2, but this is redundant with V/J plots
          # No correction based on VJ usage
          #######
          
          logo <- list()
          
          if(plot.cdr12.motif==1){
            
            if(length(countV.es[[chain]])>0){
              
              for(cdr in c("CDR1", "CDR2")){
                lg <- nchar(cdr123[[sp]][[chain]][1,cdr])  #Check the length of cdr1 and cdr2
                
                if(cdr=="CDR1"){ct <- countCDR1.es[[chain]][[lg]]; ct2 <- countCDR1.baseline[[chain]][[lg]]; ylab <- es.name; ylab.baseline <- baseline.name}
                if(cdr=="CDR2"){ct <- countCDR2.es[[chain]][[lg]]; ct2 <- countCDR2.baseline[[chain]][[lg]]; ylab <- ""; ylab.baseline <- ""}
                
                pwm <- build_cdr12_motif(ct, keep.gap=keep.gap.pwm)  #Useful if we keep the gaps
                logo1 <- ggseqlogoMOD(data=pwm, additionaAA=additionalAA, axisTextSizeX = 10, axisTextSizeY = 10) + 
                  ggtitle(paste(cdr,chain.small[chain],sep="")) + ylab(ylab) + th + theme(plot.title=element_text(size=12))
                pwm2 <- build_cdr12_motif(ct2, keep.gap=keep.gap.pwm) 
                logo2 <- ggseqlogoMOD(data=pwm2, additionaAA=additionalAA,  axisTextSizeX = 10, axisTextSizeY = 10) +
                  ggtitle(paste(cdr,chain.small[chain], sep="")) + ylab(ylab.baseline) + th + theme(plot.title=element_text(size=12))
                
                logo[[cdr]] <- ggarrange(logo1, logo2, nrow=2)
                
              }
              
            } else {
              for(cdr in c("CDR1", "CDR2")){
                logo[[cdr]] <- ggplot()
              }
            }
          }
          
          #######
          # Plot comparison of motifs of CDR3 for specific length
          #######
          
          #Take the length with max Input TCRs
          if(plot.cdr12.motif==0){
            info <- c(chain.small[chain], es.name, baseline.name)
          } else {
            info <- c(chain.small[chain], "","")
          }
          if(comp.baseline==1){
            #Here we include a correction based on VJ usage for each length.
            if(length(countV.es[[chain]])>0 & length(countJ.es[[chain]])>0 & renormVJ==1){
              info[3] <- paste(info[3], "VJ",sep=" | ")
              wt <- weighted_countCDR3(countCDR3.VJL.baseline[[chain]], countVJ.L.es[[chain]])
              CDR3 <- plotCDR3(countL.es[[chain]], countL.baseline[[chain]], countCDR3.L.es[[chain]], wt, info, comp.baseline)
            } else {
              CDR3 <- plotCDR3(countL.es[[chain]], countL.baseline[[chain]], countCDR3.L.es[[chain]], countCDR3.L.baseline[[chain]], info, comp.baseline)
            }
            
          } else {
            CDR3 <- plotCDR3(countL.es[[chain]], countL.baseline[[chain]], countCDR3.L.es[[chain]], countCDR3.L.baseline[[chain]], info, comp.baseline)
          }

          logo.CDR3.L.es <- CDR3$ES
          
          logo.CDR3.L.baseline <- CDR3$Baseline
          L.inter <- CDR3$length
          lmax <- CDR3$lmax
          
          logo[["CDR3"]] <- ggarrange(logo.CDR3.L.es[[lmax]], logo.CDR3.L.baseline[[lmax]], nrow=2)
          
          #############
          #Build the full Figure
          #############
          
          if(plot.cdr12.motif==1){
            g <- ggarrange(ld.plot, countV.plot, countJ.plot, ncol=3)
            pg.cdr12 <- ggarrange(logo[["CDR1"]], logo[["CDR2"]], ncol=2)
            pg.all[[chain]] <- ggarrange(g, ggarrange(pg.cdr12, logo[["CDR3"]], ncol=2, widths=c(1.2,1)), heights=c(1,1.5), nrow=2)
          } else {
            if(plot.oneline==1){
              pg.all[[chain]] <- ggarrange(ld.plot, countV.plot, countJ.plot,  logo[["CDR3"]], ncol=4, nrow=1)
            } else {
              pg.all[[chain]] <- ggarrange(countV.plot, countJ.plot, ld.plot, logo[["CDR3"]], ncol=2, nrow=2)
            }
          }
          
          
          
          
          if(plot.logo.length==1){
            
            tl.logo[[chain]] <- intersect(as.numeric(names(countL.es[[chain]][countL.es[[chain]]>=min.logo])), L.inter) #Currently the min.logo limitation does not apply to L.inter
            
            if(length(tl.logo[[chain]])>0){
              logo.sub <- list()
              logo.sub.baseline <- list()
              
              plotV_J.L <- list()
              ct <- 1
              for(t in tl.logo[[chain]]){
                #print(t)
                logo.sub[[ct]] <- logo.CDR3.L.es[[t]]
                logo.sub.baseline[[ct]] <- logo.CDR3.L.baseline[[t]]
                
                #Add the comparison of V/J usage
                plotV.L <- plotVJ(countV.L.es[[chain]][[t]],countV.L.baseline[[chain]][[t]], c(paste(chain,"V", sep=""), es.name, baseline.name))
                plotJ.L <- plotVJ(countJ.L.es[[chain]][[t]],countJ.L.baseline[[chain]][[t]], c(paste(chain,"J", sep=""), es.name, baseline.name))
                
                plotV_J.L[[ct]] <- ggarrange(plotV.L, plotJ.L, ncol=2, nrow=1)
                
                ct <- ct+1
                
              }
              g1 <- ggarrange(plotlist=plotV_J.L, nrow=length(tl.logo[[chain]]), ncol=1)
              g2 <- ggarrange(plotlist=logo.sub, nrow=length(tl.logo[[chain]]), ncol=1)
              g3 <- ggarrange(plotlist=logo.sub.baseline, nrow=length(tl.logo[[chain]]), ncol=1)
              
              pg.length[[chain]] <- ggarrange(plotlist=list(g1,g2,g3), nrow=1, ncol=3)

              
            } else {
              tl.logo[[chain]] <- 1
              pg.length[[chain]] <- ggplot()
            }
            
          }
          #End of the part specific for plotting data
          
        } 
        else{
          #This is the case where no CDR3 is given
          pg.all[[chain]] <- ggplot()
          pg.length[[chain]] <- ggplot()
          
          print(paste("WARNING: No CDR3",chain.small[chain]," data", sep=""))
        }
      } #End of the loop over both chains
      
      
      
      dir <- paste(output.path,"/plots/", sep="")
      if(!file.exists(dir)){
        dir.create(dir);
      }
      
      if(chain.list.output=="A"){   pg.both <- pg.all[[chain.list[1]]]; div=2   }
      if(chain.list.output=="B"){   pg.both <- pg.all[[chain.list[2]]]; div=2   }
      if(chain.list.output=="AB"){  pg.both <- ggarrange(pg.all[[chain.list[1]]], pg.all[[chain.list[2]]], ncol=2); div=1  }
      
      fig <- annotate_figure(pg.both, top = text_grob(model, face = "bold", size = 12))
      
      if(plot.cdr12.motif==1){
        width <- 18
        height <- 9 
      } else {
        if(plot.oneline==1){
          width <- 20
          height <- 3.5
        } else {
          width <- 15
          height <- 8
        }
      }
      ggsave(fig, filename=paste(output.path,"/plots/",model,".pdf", sep=""), width = width/div, height = height)
      
      if(plot.logo.length==1){
        dir <- paste(output.path,"/plots/CDR3_length/", sep="")
        if(!file.exists(dir)){
          dir.create(dir);
        }
        
        for(chain in chain.list){
          g.final <- pg.length[[chain]]; 
          mx <- length(tl.logo[[chain]])
          if(chain == "TRA" & chain.list.output!="B"){
            ggsave(g.final, filename=paste(dir, model,"_",chain,".pdf", sep=""), width = 20/div, height = 2.5*mx)
          }
          if(chain == "TRB" & chain.list.output!="A"){
            ggsave(g.final, filename=paste(dir, model,"_",chain,".pdf", sep=""), width = 20/div, height = 2.5*mx)
          }
        }
      }
    }
  }
}


