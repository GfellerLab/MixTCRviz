
#Run it from the MixTCRviz/ folder

#Do this if you do not want to install the package
devtools::load_all(".")

#Do this if you want to install the package
#devtools::build()
#install.packages("../MixTCRviz_0.0.1.tar.gz", repos=NULL)
#library(MixTCRviz)

#Compare input TCRs (specific for A0201_LLWNGPMAV) with baseline repertoire
MixTCRviz(input1="test/test1.csv", output.path="test/out/test1")

list_1 <- readRDS("test/out/test1/stats/A0201_LLWNGPMAV.rds")
list_2 <- readRDS("test/out_compare/test1/stats/A0201_LLWNGPMAV.rds")

comp <- identical(list_1, list_2)
if(comp){
  print("No problem detected")
} else {
  stop("There were some issues...")
}
