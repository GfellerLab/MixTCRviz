#' MixTCRviz: Plot TCR motifs
#'
#' MixTCRviz is an R package to display TCR motifs for a set of TCRs provided by the user.
#' Typically, the input TCRs correspond to TCRs binding to a specific epitope, or isolated in a specific experiment.
#' MixTCRviz compares V usage, J usage, CDR3 length distribution and CDR3 sequence motifs to those expected from baseline TCR repertoire.
#'
#' @param input1 csv files or data.frame with the input TCRs. Columns should
#'   ideally include: TRAV, TRAJ, cdr3_TRA, TRBV, TRBJ, cdr3_TRB, model, species
#'   * The "TRAV", "TRAJ", "TRBV", "TRBJ" should follow the IMGT
#'   nomenclature, with or without allele (see below for potential name correction).
#'   If a column is missing, empty values will be used.
#'   * The "cdr3_TRA" and "cdr3_TRB" columns should provide CDR3A/CDR3B sequences, following the standard definition (e.g., CAVNSDGQKLLF).
#'   Cases with non-amino acid characters, or length < 7 or > 22 will be not be considered (i.e., put to NA).
#'   If a column is missing, empty values will be used.
#'   * The "model" column typically describes the epitopes/experiments/classes/... Each model will be treated independently in MixTCRviz.
#'   If missing, all TCRs will be considered as coming from the same model, which can be specified in the MixTCRviz function (model.default, default value is Model_default).
#'   * The "species" column indicates the species of the TCR (e.g., NOT the source organism of the epitope). It should be "HomoSapiens" or "MusMusculus"
#'   If missing, all TCRs will be considered as coming from the same species, which can be specified in the MixTCRviz function (species.default, default value is HomoSapiens).
#' @param output.path name of the output directory (if not already existing, it
#'   will be created). If existing the files with the same name will be overwritten.
#' @param input2 (default=NULL): csv file or data.frame containing the second set of
#'   TCRs to be used in comparisons. Same format as input1. 
#'   In particular, the set of models in the "model" field need to be the same as in input1
#'   so that comparisons are performed for each model separately.
#'   If no "model" is given, all TCRs will be considered as coming from the same model ("Model_default").
#'   If input2 is provided, the comparisons is performed with this second input, and not the baseline
#'   repertoire.
#' @param baseline.file (default=NULL) .rds or .rda file containing all data about the
#'   baseline repertoire to be used, in the same format as the default repertoires. 
#'   If .rda file, the object containing the data about the repertoires must be called 'baseline'.
#'   If empty, the default baseline repertoires are used. 
#' @param use.allele (default=0)
#'    * 0: All V/J alleles are merged at the gene level (recommended).
#'    * 1: Alleles are kept. If some entries do not include alleles, the most frequent one is added.
#'   Currently, use.allele is fixed to 0 for mouse TCRs.
#' @param use.mouse.strain (default=0)
#'    * 0: Merge the different TRAV genes corresponding to different mouse strains (e.g., TRAV10D, TRAVN10 and TRAV10 become all TRAV10).
#'   This is recommended since differences between input TCRs and baseline TCRs
#'   may result from the use of different mouse strains, and not be related to
#'   specificity in input TCRs. In addition different TCR reconstruction tools
#'   use different approaches to call these genes, which can create artificial
#'   enrichment in one of them versus the baseline.
#'    * 1: The different TRAV segments (e.g., TRAV10D, TRAVN10 and TRAV10) are treated separately.
#' @param renormVJ (default=1)
#'    * 0: Compare CDR3 length distribution and motif with those from the baseline repertoire or input2
#'    * 1: Compare CDR3 length distribution and motif with those from the baseline repertoire or input2 with the V-J
#'   usage observed in the input TCRs. In the plots the mark " | VJ" is used to
#'   indicate that CDR3 length distributions and motifs correspond to those
#'   expected with the V-J usage in the input TCRs.
#' @param species.default (default="HomoSapiens"). Option to provide the species for all the TCR in input. 
#'   This is useful if your input does not contain a "species" column.
#'   In case the input contains the "species" column, species.default is not considered.
#'   Should be either "HomoSapiens" or "MusMusculus"
#' @param model.default (default="Model_default") Option to provide the model for all the TCR in input, which is also used as the name of the output files. 
#'   This is useful if your input does not contain a "model" column.
#'   In case the input contains the "model" column, model.default is not considered.
#' @param set.cdr3a.length (default=NA) Length for the CDR3a motif to be shown in the main plot.
#'   By default the value corresponding to the most frequent CDR3a length in input1 (and also present in input2 if input2 is given) is chosen.
#' @param set.cdr3b.length (default=NA) Length for the CDR3b motif to be shown in the main plot.
#'   By default the value corresponding to the most frequent CDR3b length in input1 (and also present in input2 if input2 is given) is chosen.
#' @param N.min (default=10) Minimum number of TCR (i.e., V-J-CDR3) for at least
#'   one chain. This number is computed after cleaning the data.
#' @param output.stat (default=1) Create a stat/ folder with .rds objects summarizing the raw statistics for each model.
#'   This includes countL, countV, countJ, countCDR3.L, etc. for each chain used in input.
#'   Create also a processed_data/ folder with the data for each model after the different processing (e.g., removing alleles, correcting V/J names, removing inconsistent CDR3)
#' @param correct.gene.names (default=1)
#'    * 0: Do not attempt to correct V/J gene names. Put to NA genes not in IMGT.
#'    * 1: Attempt to correct V/J gene names not in IMGT based on our internal dictionary. Put to NA genes that could not be corrected
#' @param check.cdr3.mode (default=1)
#'    * 0: Keep all CDR3 without any correction.
#'    * 1: Remove V and CDR3 when the first CDR3 amino acid is incompatible with the V segment; 
#'        Remove J and CDR3 when the last two amino acids are not compatible with the J segments.
#'    * 2: remove V and CDR3 when the M first CDR3 amino acids are incompatible with the V segment, with M depending on the V gene and the CDR3 length; 
#'        Remove J and CDR3 when the last two amino acids are not compatible with the J segments, with M depending on the J gene and CDR3 length.
#'    * Note: 0: all PCR / sequencing / TCR reconstruction errors are kept. 
#'            1: several PCR / sequencing / TCR reconstruction errors are removed, but only based on the first and last amino acids in CDR3.
#'            2: most PCR / sequencing / TCR reconstruction errors are removed.
#' @param verbose (default=1)
#'    * 0: Do not write any QC in the output
#'    * 1: Write max 10 examples of putative issues with the data (V/J names, CDR3 sequences, etc.) in the terminal
#'    * 2: Write all the putative issues with the data (V/J names, CDR3 sequences, etc.) in the terminal
#' @param plot (default=1)
#'    * 0: only create .rds object with the statistics
#'    * 1: Plot the data in output.path/plots/ and create .rds object with the statistics for each model in output.path/stats/.
#' @param plot.cdr12.motif (default=0)
#'    * 0: Only show motifs for CDR3.
#'    * 1: Include sequence motifs for CDR1 and CDR2. With this choice, plot.oneline is set to 0.
#' @param plot.oneline (default=0)
#'    * 0: Show the data on two lines (better for clarity).
#'    * 1: Show all plots in a single line (can be useful to compare different models).
#'    * 2: Show only V/J usage and length (i.e., do not show CDR3 motifs)
#' @param plot.logo.length (default=0)
#'    * 0: Show only the CDR3 motifs for the most frequent CDR3 length.
#'    * 1: Show in a separate plot the V usage, J usage and CDR3 motifs for all CDR3 length, for both alpha and beta chains.
#' @param plot.cdr3.norm (default=0)
#'    * 0: Show the CDR3 motifs of the baseline repertoire.
#'    * 1: Show the CDR3 motifs of the input TCRs after subtracting the baseline repertoire.
#'    * 2: Show the CDR3 motifs of the input TCRs after normalising by the baseline repertoire (motif of normalised fold-change).
#' @param chain.list.output (default="AB")
#'    * A: Only the alpha chain is plotted in output;
#'    * B: only the beta chain is plotted in output;
#'    * AB both chains are plotted in output
#' @param output.format (default="pdf"): Choose the format for the output
#'      plots (can be "pdf", "png" or "jpg").
#' @param input1.name (default="Input"): Provide a generic name for
#'   the input TCRs in the plots (e.g., Epitope Specific). Avoid names with more
#'   than 20 characters
#' @param input2.name (default=NULL): If a second set of TCRs is provided
#'   (i.e., input2 != NULL), Provide a generic name for the input2 TCRs in the
#'   plots. Avoid names with more than 20 characters.
#'
#' @returns Nothing.
#' @export
MixTCRviz <- function(input1, output.path,
                      input2=NULL, baseline.file=NULL,
                      use.allele=0, correct.gene.names=1, use.mouse.strain=0, check.cdr3.mode=1,
                      renormVJ=1, N.min=10, output.stat=1, set.cdr3a.length=NA, set.cdr3b.length=NA,
                      species.default="HomoSapiens", model.default="Model_default", verbose=1,
                      plot=1, plot.cdr12.motif=0, plot.oneline=0, plot.logo.length=0, plot.cdr3.norm=0,
                      chain.list.output="AB", input1.name="Input", input2.name=NULL, output.format="pdf"){
  
  
  #######
  # Choose some parameters
  #######
  
  if(!is.na(set.cdr3a.length)){
    if(is.numeric(set.cdr3a.length)==F | set.cdr3a.length < Lmin | set.cdr3a.length > Lmax | set.cdr3a.length%%1 != 0){
      print(paste("Invalid value for set.cdr3a.length",". Default value will be used.", sep=""))
      set.cdr3a.length=NA
    }
  }
  if(!is.na(set.cdr3b.length)){
    if(is.numeric(set.cdr3b.length)==F | set.cdr3b.length < Lmin | set.cdr3b.length > Lmax | set.cdr3b.length%%1 != 0){
      print(paste("Invalid value for set.cdr3b.length=",set.cdr3b.length,". Default value will be used.", sep=""))
      set.cdr3b.length=NA
    }
  }
  
  if(chain.list.output=="A"){
    chain.list <- c("TRA"); 
    set.cdr3.length <- c(set.cdr3a.length)
  } else if(chain.list.output=="B"){
    chain.list <- c("TRB"); 
    set.cdr3.length <- c(set.cdr3b.length)
  } else if(chain.list.output=="AB"){
    chain.list <- c("TRA","TRB"); 
    set.cdr3.length <- c(set.cdr3a.length, set.cdr3b.length)
  } else { 
    stop(paste("chain.list.output ", chain.list.output, " not supported by MixTCRviz", sep=""))
  }

  names(set.cdr3.length) <- chain.list
  
  col.TCR <- c(); segment.list <- c()
  for(chain in chain.list){
    col.TCR <- c(col.TCR, paste(chain,"V",sep=""), paste(chain,"J",sep=""), paste("cdr3_",chain,sep=""))
    segment.list <- c(segment.list, paste(chain,"V",sep=""), paste(chain,"J",sep=""))
  }
  
  if(output.format != "pdf" & output.format != "png" & output.format != "jpg"){
    stop("'output.format' should be pdf, png or jpg!")
  }
  
  if(species.default %in% species.list == F){
    stop("Wrong choice for species.default. Should be either \"HomoSapiens\" or \"MusMusculus\"")
  }
  
  if(check.cdr3.mode %in% c(0,1,2) == F){
    check.cdr3.mode <- 1
    print("Using the standard mode to clean CDR3 sequences")
  }
  if(plot.oneline %in% c(0,1,2)==F){
    print(paste("plot.oneline=",plot.oneline," not supported, using default plot.oneline=0", sep=""))
    plot.oneline <- 0
  }
  
  
  keep.gap.pwm <- 0 # 1: means that 'g' (gaps in CDR1/2) are treated as an additional aa in the logos. 0: means that 'g' are treated as unspecific (i.e., 0.05) in the logos
  if(keep.gap.pwm==1){taa.list <- c(aa.list,"g"); additionalAA <- "g"} else {taa.list <- aa.list; additionalAA <- ""}
  
  min.logo <- 5 #Minimum number of sequences to plot the logos (when plotting different lengths)
  
  if(plot.cdr12.motif==1){ plot.oneline <- 0 }
  if(plot.oneline != 0){
    th <- theme(plot.title = element_text(size = 8, hjust=0.5), axis.title=element_text(size=4))
  }
  
  if(plot.cdr3.norm != 0 & plot.cdr3.norm != 1 & plot.cdr3.norm != 2){
    plot.cdr3.norm <- 0
  }
  
  if(N.min < 1 | is.numeric(N.min)==FALSE){N.min <- 1}
  
  if(!dir.exists(output.path)){
    dir.create(output.path, recursive = TRUE);
  }
  es.name <- input1.name
  
  
  if(is.null(input2)){
    comp.baseline <- 1
    if(is.null(input2.name)){
      baseline.name <- "Baseline"
    } else if(is.character(input2.name)){ 
      baseline.name <- input2.name
    } else {
      stop("Invalid value for input2.name. Should be a character or NULL")
    }
  } else if(is.character(input2) | is.data.frame(input2)) {
    comp.baseline <- 0
    if(is.null(input2.name)){ 
      baseline.name <- "Input2"
    } else if(is.character(input2.name)){ 
      baseline.name <- input2.name
    } else {
      stop("Invalid value for input2.name. Should be a character or NULL")
    }
  } else {
    stop("Invalid format for input2. Should be a filename or a data.frame, or left empty")
  }
  
  ############################
  # Load all the input data
  ############################
  
  if(is.character(input1)==T){
    if(file.exists(input1)){
      es.all <- read.csv(input1)
    } else {
      stop("Missing file for input1")
    }
  } else if(is.data.frame(input1)==T){
    es.all <- input1
  } else {
    stop("Invalid value for input1. Should be a .csv filename or a data.frame")
  }
  

  # Check the input
  print("Check input1")
  es.all <- check_input(es.all, chain.list.output, "input1", species.default, model.default)
  es.all <- clean_input(es.all, use.allele, correct.gene.names, use.mouse.strain, chain.list.output, species.default, check.cdr3.mode, verbose)
  
  #############
  # Load input2
  #############
  
  if(comp.baseline==0){
    if(is.character(input2)==T){
      if(file.exists(input2)){
        es2.all <- read.csv(input2)
      } else {
        stop("Missing file for input2")
      }
    } else if (is.data.frame(input2)==T){
      es2.all <- input2
    }
    print("Check input2")
    es2.all <- check_input(es2.all, chain.list.output, "input2", species.default, model.default)
    es2.all <- clean_input(es2.all, use.allele, correct.gene.names, use.mouse.strain, chain.list.output, species.default, check.cdr3.mode, verbose)
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
    if( (nA>=N.min & chain.list.output=="A") |(nB>=N.min & chain.list.output=="B") | ((nA>=N.min | nB>=N.min) & chain.list.output=="AB")){
      return(1)
    } else {
      warning("Model ",x," will not be considered (less than N.min=",
              N.min," data)")
      return(0)
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
    
    print(paste("Model:",model))
    
    pg.all <- list()
    pg.length <- list()
    tl.logo <- list() #Keep track of the number of logos of different length
    
    ######
    # Select the input models with enough data
    ######
    
    es <- es.all[which(es.all[,"model"]==model),]
    sp <- unique(es[,"species"])
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
    }
   
    
    count <- build_stat(es=es, chain.list=chain.list, sp=sp, comp.VJL=0)
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
    
    s <- unlist(strsplit(model, split="_"))
    MHC <- s[1]
    epitope <- s[2]
    
    info <- c(model,sp,MHC,epitope)
    names(info) <- c("model", "species", "MHC", "epitope")
    
    summary <- list(info, L.es, countL.es, countV.es, countJ.es, countV.L.es, countJ.L.es, countCDR1.es, countCDR2.es, countCDR3.L.es, countVJ.es, countVJ.L.es)
    names(summary) <- c("info", "L", "countL", "countV", "countJ", "countV.L", "countJ.L", "countCDR1", "countCDR2", "countCDR3.L", "countVJ", "countVJ.L")
    if(output.stat==1){
      dir <- paste(output.path,"/stats/", sep="")
      if(!dir.exists(dir)){
        dir.create(dir);
      }
      saveRDS(summary, file=paste(output.path,"/stats/",model,".rds", sep=""))
      dir <- paste(output.path,"/processed_data/", sep="")
      if(!dir.exists(dir)){
        dir.create(dir);
      }
      write.csv(es, file=paste(output.path,"/processed_data/",model,".csv", sep=""), quote=F, row.names = F, na = "")
    }
    
    if(plot==1){
      
      if(comp.baseline==1){
        
        #####
        # Load the baseline repertoire corresponding to the species
        #####
        
        if(is.null(baseline.file)){
          
          #These are the default repertoires.
          if(sp=="HomoSapiens"){
            if(use.allele.es==1){
              baseline <-  MixTCRviz::summary_HomoSapiens
            } else {
              baseline <- MixTCRviz::summary_HomoSapiens_noallele
            }
          }
          if(sp=="MusMusculus"){
            if(use.allele.es==0){  #N.B. Currently always the case in mouse
              if(use.mouse.strain==1){
                baseline <- MixTCRviz::summary_MusMusculus_noallele_SEQTR
              } else {
                baseline <- MixTCRviz::summary_MusMusculus_noallele_noStrain_SEQTR
              }
            }
          }
        } else {
          if(file.exists(baseline.file)){
            st <- unlist(strsplit(baseline.file, split=".", fixed = T))
            if(st[length(st)]=="rds"){
              baseline <- readRDS(file=baseline.file)
            } else if(st[length(st)]=="rda" | st[length(st)]=="rdata"){
              load(baseline.file)
            } else{
              stop("Unsupported format for baseline.file. Should be .rds, .rda or .rdata")
            }
          } else {
            stop("Missing baseline file")
          }
        }
        
        L.baseline <- baseline$L
        countL.baseline <- baseline$countL
        countV.baseline <- baseline$countV
        countJ.baseline  <- baseline$countJ
        countV.L.baseline <- baseline$countV.L
        countJ.L.baseline  <- baseline$countJ.L
        countL.VJ.baseline  <- baseline$countL.VJ
        countVJ.baseline  <- baseline$countVJ
        countCDR1.baseline  <- baseline$countCDR1  
        countCDR2.baseline <- baseline$countCDR2  
        countCDR3.L.baseline <- baseline$countCDR3.L
        countCDR3.VJL.baseline <- baseline$countCDR3.VJL
        
      } else {
        
        #####
        # Load the other repertoire to compare with
        #####
        
        ind <- which(es2.all[,"model"]==model)
        if(length(ind)==0){
          stop(c("Missing model in Input2: ", model))
        } else {
          es2 <- es2.all[ind,]
        }
        
        count <- build_stat(es=es2, chain.list=chain.list, sp=sp, comp.VJL=renormVJ)
        
        L.baseline <- count$L
        countL.baseline <- count$countL
        countV.baseline <- count$countV
        countJ.baseline  <- count$countJ
        countV.L.baseline <- count$countV.L
        countJ.L.baseline  <- count$countJ.L
        countCDR1.baseline  <- count$countCDR1
        countCDR2.baseline <- count$countCDR2
        countCDR3.L.baseline <- count$countCDR3.L
        
        if(renormVJ==1){
          countL.VJ.baseline <- count$countL.VJ
          countCDR3.VJL.baseline <- count$countCDR3.VJL
        }
        
      }
      
      for(chain in chain.list){
        
        #print(chain)
        
        if(comp.baseline==1){
          #Check segments that were in the ES, but not in baseline
          miss.V.baseline <- setdiff(names(countV.es[[chain]]), names(countV.baseline[[chain]]))
          miss.J.baseline <- setdiff(names(countJ.es[[chain]]), names(countJ.baseline[[chain]]))
          if(verbose>0){
            if(length(miss.V.baseline)>=1){
              print(paste("WARNING: ",chain,"V in Input TCRs, but absent from baseline: ", sep=""))
              print(miss.V.baseline)  
            }
            if(length(miss.J.baseline)>=1){
              print(paste("WARNING: ",chain,"J in Input TCRs, but absent from baseline: ", sep=""))
              print(miss.J.baseline)
            }
          }
        }
        
        #Make sure there are CDR3 sequences
        if(length(countL.es[[chain]])>0){
          
          if(verbose>0){
            if(length(countV.es[[chain]])==0){
              print(paste("WARNING: No ", chain,"V segment in input1", sep=""))
            }
            if(length(countJ.es[[chain]])==0){
              print(paste("WARNING: No ",chain,"J segment in input1", sep=""))
            }
            
            if(comp.baseline==0){
              if(length(countV.baseline[[chain]])==0){
                print(paste("WARNING: No ", chain,"V segment in input2", sep=""))
              }
              if(length(countJ.baseline[[chain]])==0){
                print(paste("WARNING: No ",chain,"J segment in input2", sep=""))
              }
            }
          }
          
          info <- c(chain.small[chain], paste(es.name, " (", sum(countL.es[[chain]]),")",sep=""), baseline.name)
          if(comp.baseline==0){
            info[3] <- paste(info[3], " (", sum(countL.baseline[[chain]]),")",sep="")
          }
          
          if(renormVJ==1){
            if(length(countVJ.es[[chain]])>0){
              info[3] <- paste(info[3], "VJ",sep=" | ")
              bs <- weighted_countL(countL.VJ.baseline[[chain]], countVJ.es[[chain]])
            } else {
              print(paste("No V-J information to compute baseline CDR3",chain.small[chain]," length distribution | VJ. Check your data or use renormVJ=0", sep=""))
              info[3] <- paste(info[3], "NA",sep=" | ")
              bs <- countL.baseline[[chain]]-countL.baseline[[chain]]
            }
          } else{
            bs <- countL.baseline[[chain]]
          }
          
          ld.plot <- plotLD(countL.es[[chain]], bs, info, plot.oneline)
          
          #######
          # Plot comparison of V/J usage
          #######
          
          infoV <- c(paste(chain,"V", sep=""), es.name, baseline.name)
          infoJ <- c(paste(chain,"J", sep=""), es.name, baseline.name)
          
          if(length(countV.es[[chain]])>0){
            countV.plot <- plotVJ(countV.es[[chain]], countV.baseline[[chain]], infoV, comp.baseline)
          } else {
            countV.plot <- ggplot()
          }
          if(length(countJ.es[[chain]])>0){
            countJ.plot <- plotVJ(countJ.es[[chain]], countJ.baseline[[chain]], infoJ, comp.baseline)
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
                lc <- paste("L",nchar(cdr123[[sp]][[chain]][1,cdr]),sep="_")  #Check the length of cdr1 and cdr2
                if(cdr=="CDR1"){
                  ct <- countCDR1.es[[chain]][[lc]];
                  ylab <- ""; ylab.baseline <- ""
                }
                if(cdr=="CDR2"){
                  ct <- countCDR2.es[[chain]][[lc]];
                  ylab <- ""; ylab.baseline <- ""
                }
                pwm <- build_cdr12_motif(ct, keep.gap=keep.gap.pwm)  #Useful if we keep the gaps
                logo1 <- ggseqlogoMOD(data=pwm, additionaAA=additionalAA, axisTextSizeX = 10, axisTextSizeY = 10) +
                  ggtitle(paste(cdr,chain.small[chain]," ",es.name," (", sum(countV.es[[chain]]),")",sep="")) + ylab(ylab) + th + theme(plot.title=element_text(size=12))
                #pwm2 <- build_cdr12_motif(ct2, keep.gap=keep.gap.pwm)
                #logo2 <- ggseqlogoMOD(data=pwm2, additionaAA=additionalAA,  axisTextSizeX = 10, axisTextSizeY = 10) +
                #ggtitle(paste(cdr,chain.small[chain]," ",baseline.name," ", sep="")) + ylab(ylab.baseline) + th + theme(plot.title=element_text(size=12))
                #logo[[cdr]] <- ggarrange(logo1, logo2, nrow=2)
                logo[[cdr]] <- ggarrange(logo1, nrow=1)
                
              }
              
            } else {
              for(cdr in c("CDR1", "CDR2")){
                logo[[cdr]] <- ggplot()
              }
            }
          }
          
          #######
          # Plot comparison of motifs of CDR3 for specific lengths
          #######
          
          #Take the length with max Input TCRs
          info <- c(chain.small[chain], es.name, baseline.name)
         
          #Here we include a correction based on VJ usage for each length.
          if(renormVJ==1){
            if(max(sapply(countVJ.L.es[[chain]],length))>0){  #Make sure we have V-J pairs for at least one length
              info[3] <- paste(info[3], "VJ",sep=" | ")
              bs <- weighted_countCDR3(countCDR3.VJL.baseline[[chain]], countVJ.L.es[[chain]])
            } else {
              print(paste("No V-J information to compute baseline CDR3",chain.small[chain]," motif | VJ. Check your data or use renormVJ=0",sep=""))
              info[3] <- paste(info[3], "NA",sep=" | ")
              bs <- lapply(countCDR3.L.baseline[[chain]], function(x){y <- x-x+0.05; return(y)})
            }
          } else {
            bs <- countCDR3.L.baseline[[chain]]
          }
          CDR3 <- plotCDR3(countL.es[[chain]], countL.baseline[[chain]], countCDR3.L.es[[chain]], 
                           bs, info, comp.baseline, plot.oneline, plot.logo.length, plot.cdr3.norm,
                           set.cdr3.length[[chain]])
          
          logo.CDR3.L.es <- CDR3$ES
          logo.CDR3.L.baseline <- CDR3$Baseline
          
          if(length(CDR3$length)>0){
            L.inter <- paste("L",CDR3$length, sep="_")
          } else {L.inter <- c()}
          
          lmax <- CDR3$lmax
          
          logo[["CDR3"]] <- ggarrange(CDR3$ES_max, CDR3$Baseline_max, nrow=2)
          
          #############
          #Build the full Figure
          #############
          
          if(plot.cdr12.motif==1){
            g <- ggarrange(countV.plot, countJ.plot, ld.plot, ncol=3)
            pg.cdr12 <- ggarrange(logo[["CDR1"]], logo[["CDR2"]], ncol=2)
            pg.all[[chain]] <- ggarrange(g, ggarrange(pg.cdr12, CDR3$ES_max, ncol=2, widths=c(1.2,1)), heights=c(1.5,1), nrow=2)
          } else {
            if(plot.oneline==0){
              pg.all[[chain]] <- ggarrange(countV.plot, countJ.plot, ld.plot, logo[["CDR3"]], ncol=2, nrow=2)
            } else if(plot.oneline==1){
              pg.all[[chain]] <- ggarrange(countV.plot, countJ.plot, ld.plot, logo[["CDR3"]], ncol=4, nrow=1)
            } else if(plot.oneline==2){
              pg.all[[chain]] <- ggarrange(countV.plot, countJ.plot, ld.plot, ncol=3, nrow=1)
            }
          }
          
          
          if(plot.logo.length==1){
            
            tl.logo[[chain]] <- intersect(names(countL.es[[chain]][countL.es[[chain]]>=min.logo]), L.inter) #Currently the min.logo limitation does not apply to L.inter
            
            if(length(tl.logo[[chain]])>0){
              logo.sub <- list()
              logo.sub.baseline <- list()
              
              plotVJ.L <- list()
              ct <- 1
              for(t in tl.logo[[chain]]){
                logo.sub[[ct]] <- logo.CDR3.L.es[[t]]
                logo.sub.baseline[[ct]] <- logo.CDR3.L.baseline[[t]]
                
                #Add the comparison of V/J usage
                plotV.L <- plotVJ(countV.L.es[[chain]][[t]],countV.L.baseline[[chain]][[t]], c(paste(chain,"V", sep=""), es.name, baseline.name), comp.baseline)
                plotJ.L <- plotVJ(countJ.L.es[[chain]][[t]],countJ.L.baseline[[chain]][[t]], c(paste(chain,"J", sep=""), es.name, baseline.name), comp.baseline)
                
                plotVJ.L[[ct]] <- ggarrange(plotV.L, plotJ.L, ncol=2, nrow=1)
                
                ct <- ct+1
                
              }
              g1 <- ggarrange(plotlist=plotVJ.L, nrow=length(tl.logo[[chain]]), ncol=1)
              g2 <- ggarrange(plotlist=logo.sub, nrow=length(tl.logo[[chain]]), ncol=1)
              g3 <- ggarrange(plotlist=logo.sub.baseline, nrow=length(tl.logo[[chain]]), ncol=1)
              
              pg.length[[chain]] <- ggarrange(plotlist=list(g1,g2,g3), nrow=1, ncol=3)
              
              
            } else {
              tl.logo[[chain]] <- c()
              pg.length[[chain]] <- ggplot()
            }
            
          }
          #End of the part specific for plotting data
          
        } else {
          #This is the case where no CDR3 is given
          pg.all[[chain]] <- ggplot()
          pg.length[[chain]] <- ggplot()
          tl.logo[[chain]] <- c()
          
          print(paste("WARNING: No CDR3",chain.small[chain]," data in input1", sep=""))
        }
      } #End of the loop over both chains
      
      
      
      dir <- paste(output.path,"/plots/", sep="")
      if(!dir.exists(dir)){
        dir.create(dir);
      }
      
      if(chain.list.output=="A" | chain.list.output=="B"){   pg.both <- pg.all[[chain.list[1]]]; div=2   }
      if(chain.list.output=="AB"){  pg.both <- ggarrange(pg.all[[chain.list[1]]], pg.all[[chain.list[2]]], ncol=2); div=1  }
      
      #fig <- annotate_figure(pg.both, top = text_grob(model, face = "bold", size = 12))
      fig <- pg.both
      
      if(plot.cdr12.motif==1){
        width <- 20
        height <- 6
      } else {
        if(plot.oneline==0){
          width <- 15
          height <- 8
        } else if(plot.oneline==1){
          width <- 20
          height <- 3
        } else if(plot.oneline==2){
          width <- 15
          height <- 3
        }
      }
      ggsave(fig, filename=paste(output.path,"/plots/",model,".", output.format, sep=""), device=output.format, width = width/div, height = height)
      
      if(plot.logo.length==1){
        dir <- paste(output.path,"/plots/CDR3_length/", sep="")
        if(!dir.exists(dir)){
          dir.create(dir);
        }
        
        for(chain in chain.list){
          g.final <- pg.length[[chain]];
          mx <- length(tl.logo[[chain]])
          if(mx>0){
            ggsave(g.final, filename=paste(dir, model,"_",chain,".", output.format, sep=""), device=output.format, width = 20, height = 2.5*mx)
          }
        }
      }
    }
  }
}

