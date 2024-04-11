

source("MixTCRviz.R")

#Compare input TCRs (specific for A0201_LLWNGPMAV) wiht baseline repertoire
MixTCRviz(input1="test/test1.csv", output.path="test/out_compare/test1")

#Compare input TCRs (specific for A0201_LLWNGPMAV) with another set of TCRs recognizing the same epitope
MixTCRviz(input1="test/test1.csv",input2="test/test2.csv", output.path="test/out_compare/test2")


