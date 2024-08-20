
#Run it from the MixTCRviz folder

#Do this if you do not want to install the package
devtools::load_all(".")

#Do this if you have installed the package
library(MixTCRviz)

#Compare input TCRs (specific for A0201_LLWNGPMAV) with baseline repertoire
MixTCRviz::MixTCRviz(input1="test/test1.csv", output.path="test/out/test1")

