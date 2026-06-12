
#Run it from the MixTCRviz/ folder

#Do this if you do not want to install the package
#devtools::load_all(".")

#Do this if you have installed the package
library(MixTCRviz)

#Compare input TCRs (specific for A0201_LLWNGPMAV) with baseline repertoire
cat("Building a standard TSP\n")
MixTCRviz(input1="test/test.csv", output.path="test/out/", output.stat = T)

list_1 <- readRDS("test/out/stats/A0201_LLWNGPMAV.rds")
list_2 <- readRDS("test/out_compare/stats/A0201_LLWNGPMAV.rds")

#Test MixTCRviz including V/J allele (NOT recommended unless the baseline comes from the same donor as the data in input1).
use.allele <- T
if(use.allele){
  cat("Building TSP with allele information (not recommanded without a baseline from the same donor)\n")
  MixTCRviz(input1="test/test_allele.csv", output.path="test/out/", use.allele=T, model.default = "A0201_LLWNGPMAV_allele")
}

#Test the MixTCRviz with a specific baseline, including V/J allele information.
#If input2 is to be used as a baseline, remember to set renormVJ to TRUE.
use.own.baseline <- T
if(use.own.baseline){
  cat("\nBuilding TSP with allele information and a specific baseline (will take roughly 20 secondes)\n")
  MixTCRviz(input1="test/test_allele.csv", input2="test/baseline_allele.csv", output.path="test/out/", use.allele=T, 
            renormVJ = T, model.default = "A0201_LLWNGPMAV_allele_baseline")
}

comp <- identical(list_1, list_2)
if(comp){
  print("No problem detected")
} else {
  stop("There were some issues...")
}

