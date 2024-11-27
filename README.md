
# MixTCRviz
R package to visualize TCR binding motifs.

MixTCRviz can be used freely by academic groups for non-commercial purposes (see license).
The product is provided free of charge, and, therefore, on an "as is"
basis, without warranty of any kind.

**FOR-PROFIT USERS**
If you plan to use MixTCRviz or any data provided with the script in any for-profit
application, you are required to obtain a separate license.
To do so, please contact Nadette Bulgin (nbulgin@lcr.org) at the Ludwig Institute for Cancer Research Ltd.

For scientific questions, please contact David Gfeller (David.Gfeller@unil.ch)

Copyright (2024) David Gfeller

*STILL CONFIDENTIAL - DO NOT SHARE OUTSIDE OF THE GFELLER LAB*



# INTRODUCTION

MixTCRviz is a tool to visualize important properties of a set of TCRs (e.g., epitope-specific TCRs).
It focuses on V usage, J usage, CDR3 length distribution and CDR3 motifs, for both the alpha and the beta chains.
These define the so-called **TCR binding motifs**.

By default, properties of the input TCRs are compared to those of a "baseline" TCR repertoire (from HomoSapiens or MusMusculus). Alternatively, users can choose to compare with another set of TCRs ("input2").

TCR binding motifs enable users to rapidly visualize and understand what are the key determinants of specificity within a set of TCRs, such as those binding to a given epitope.


# INSTALLATION

There are different ways of "installing" MixTCRviz:
  1) It can be directly installed from the GitHub page (you need to have
    credentials correctly set so that Rstudio can access the private repository):<br />
    ` r` <br />
    `# install.packages("devtools") # Only needed if you don't already have devtools` <br />
    `devtools::install_github("GfellerLab/MixTCRviz")` <br />

  2) You can download the MixTCRviz directory from the GitHub page and then open
    Rstudio setting its working directory as MixTCRviz folder. Then you can
    compile it and install it:<br />
    ` r` <br />
    `devtools::build()` <br />
    `install.packages("../MixTCRviz_0.0.2.tar.gz", repos=NULL)` <br />


  3) If you don't need to use it through Rscript and are in process of developing
    the package, it would be tiring to always build it, install it, reload it,
    etc. Instead, you can set the working directory to the MixTCRviz folder and
    then simply:<br />
    ` r` <br />
    `devtools::load_all()` <br />
    You can reuse this command everytime changes are made to the code so that
    the updates are directly available in the R environment. (You need to reuse
    this load_all in every new R session to use this in dev code instead of the
    installed package). And you should of course not do a "library(MixTCRviz)"
    if using the load_all, this command makes as if you had done such one (the
    library... would probably load the built+installed package not the
    development version).

You may be prompted to install several packages.

### Testing the INSTALLATION:

In the MixTCRviz directory run:

`Rscript test_MixTCRviz.R`

Alternatively, you can run the code test_MixTCRviz.R in Rstudio (or any R interface)

The output in test/out should be the same as in test/out_compare

# RUNNING


MixTCRviz should be primarily run in R, by calling the MixTCRviz function. To
this end, in R, you need to call the MixTCRviz function with the mandatory
parameters (e.g., MixTCRviz::MixTCRviz(input1=input, output.path=output)).
(note that the "MixTCRviz::" isn't needed if using the "load_all"
version from above or if you used library(MixTCRviz)).


## Mandatory parameters:

 - input1: Input1 can one of the different values:
   1) A .csv or .txt or .tsv file with the input TCRs
   2) A data.frame with the input TCRs
   3) A list in the MixTCRviz format.

   If using a filename or a data.frame:
    * Columns should ideally consist of "TRAV","TRAJ","cdr3_TRA","TRBV","TRBJ","cdr3_TRB"
    * If TCRs from multiple experiments/epitopes/classses/... are provided in the same file, the "model" column should indicate the models.
    ("Model_default" is used by default if no "model" column is provided and model.default=NULL, see below). The use of multiple model is only possible in input1.
    * If TCRs from multiple species are provided in the same file, the "species" should indicate the species of the TCRs
    ("HomoSapiens" or "MusMusculus", with "HomoSapiens" being the default if no "species" column and species.default=NULL)
    * The "TRAV", "TRAJ", "TRBV", "TRBJ" entries should follow the IMGT
   nomenclature, with or without allele (see below for potential name correction).
    * The "cdr3_TRA" and "cdr3_TRB" columns should provide CDR3A/CDR3B sequences, following the standard definition (e.g., CAVNSDGQKLLF).
   Cases with non-amino acid characters, or length < 7 or > 22 will be not be considered (i.e., put to NA).
   * Other formats are supported (see below)

 - output.path: name of the output directory (if not already existing, it
   will be created). If existing the files with the same name will be overwritten.
   It can be left empty ONLY with the option return.object=2.

## Optional parameters:


 - input2 (default=NULL): .csv file or data.frame containing a second set of
   TCRs to be used in comparisons. Same format as input1, except that all TCRs are assumed to come from a single model so any data in the "model" field is omitted.
   If input2 is provided, the comparisons is performed with this second input, and not the baseline
   repertoire. In this case, all models provided in input1 will be compared with the repertoire in input2, and multiple species should not be mixed in input1.
   Other than .csv filename or data.frame, a third alternative consists of a list or .rds with precomputed statistics, such as those generated by MixTCRviz in output.path/stats/.rds files.

 - baseline (default=NULL) R object, or .rds/.rda/.rdata file containing all data about the
   baseline repertoire to be used, in the same format as the default repertoires.
   If .rda/.rdata file, the object containing the data about the repertoires must be called 'baseline'.
   Alternatively one of the strings: "Default" or "SEQTR" for data in input1 generated by publicly available TCR-seq approaches (Default) or the SEQTR protocol (SEQTR)
   If empty, the default baseline repertoires are used.

 - use.allele (default=FALSE)
    * FALSE: All V/J alleles are merged at the gene level (recommended).
    * TRUE: Alleles are kept. If some entries do not include alleles, the most frequent one is added.
   Currently, use.allele is fixed to FALSE for mouse TCRs.

 - correct.gene.names (default=TRUE)
    * F: Do not attempt to correct V/J gene names. Put to NA genes not in IMGT.
    * T: Attempt to correct V/J gene names not in IMGT based on our internal dictionary. Put to NA genes that could not be corrected

 - use.mouse.strain (default=FALSE)
    * FALSE: Merge the different TRAV genes corresponding to different mouse strains (e.g., TRAV10D, TRAVN10 and TRAV10 become all TRAV10).
   This is recommended since differences between input TCRs and baseline TCRs
   may result from the use of different mouse strains, and not be related to
   specificity in input TCRs. In addition different TCR reconstruction tools
   use different approaches to call these genes, which can create artificial
   enrichment in one of them versus the baseline.
    * TRUE: The different TRAV segments (e.g., TRAV10D, TRAVN10 and TRAV10) are treated separately.

 - check.cdr3.mode (default=1)
    * 0: Keep all CDR3 without any correction.
    * 1: Remove V and CDR3 when the first start.lg CDR3 amino acid are incompatible with the V segment;
        Remove J and CDR3 when the last end.lg amino acids are not compatible with the J segments.
        Allow compatibility check with any allele.

 - start.lg (default=1): Number of amino acids to be checked for compatibility between V segment and beginning of CDR3. Needs to be an integer between 0 and 3

 - end.lg (default=2): Number of amino acids to be checked for compatibility between J segment and end of CDR3. Needs to be an integer between 0 and 5

 - renormVJ (default=NULL)
    * NULL: If left empty, 1 is used if comparison is performed with a baseline (input2==NULL)
      and 0 is used if comparison is performed with input2 (i.e., input2 != NULL)
    * 0: Compare CDR3 length distribution and motif with those from the baseline repertoire or input2
    * 1: Compare CDR3 length distribution and motif with those from the baseline repertoire with the V-J
   usage observed in the input TCRs. In the plots the mark " | P(VJ)" is used to
   indicate that CDR3 length distributions and motifs correspond to those
   expected with the V-J usage in the TCRs used for training.

 - N.min (default=10) Minimum number of TCR (i.e., V-J-CDR3) for at least
   one chain. This number is computed after cleaning the data.

 - output.stat (default=FALSE)
    * FALSE: Do not print the stat to a file
    * TRUE: create a output.path/stat/ folder with .rds files summarizing the raw statistics for each model. Only supported if return.object != 2


 - output.processed.data (default=FALSE)
    * FALSE: Do not print the processed/clean data to a file
    * TRUE: create a output.path/processed_data/ folder with the data for each model after the different processing steps (e.g., removing alleles, correcting V/J names, removing inconsistent CDR3, etc. depending on the cleaning option selected). Only supported if return.object != 2

 - set.cdr3a.length (default=NA) Length for the CDR3a motif to be shown in the main plot.
   By default the value corresponding to the most frequent CDR3a length in input1 (and also present in input2 if input2 is given) is chosen.

 - set.cdr3b.length (default=NA) Length for the CDR3b motif to be shown in the main plot.
   By default the value corresponding to the most frequent CDR3b length in input1 (and also present in input2 if input2 is given) is chosen.

 - species.default (default="HomoSapiens"). Option to provide the species for all the TCR in input.
   This is useful if your input does not contain a "species" column.
   In case the input contains the "species" column, species.default is not considered.
   Should be either "HomoSapiens" or "MusMusculus"

 - model.default (default="Model_default") Option to provide the model for all the TCR in input, which is also used as the name of the output files.
   This is useful if your input does not contain a "model" column.
   In case the input contains the "model" column, model.default is not considered.

 - verbose (default=1)
    * 0: Do not write any QC in the output
    * 1: Write max 10 examples of putative issues with the data (V/J names, CDR3 sequences, etc.) in the terminal
    * 2: Write all the putative issues with the data (V/J names, CDR3 sequences, etc.) in the terminal
    * 3: Write all putatitive issues + cases that were corrected in the terminal

 - plot (default=TRUE)
    * F: Do not build the motifs. This overrides plot.all.length=T and interactive.plots=T
    * T: Build the motifs

 - plot.cdr12.motif (default=FALSE)
    * F: Only show motifs for CDR3 of the most frequent length.
    * T: Include sequence motifs for CDR1 and CDR2. With this choice, plot.oneline is set to 0.

 - plot.oneline (default=0)
    * 0: Show the data on two lines (better for clarity).
    * 1: Show all plots in a single line (can be useful to compare different models).
    * 2: Show only V/J usage and length (i.e., do not show CDR3 motifs for the most frequent CDR3 length)

 - plot.all.length (default=FALSE)
    * F: Show only the CDR3 motifs for the most frequent CDR3 length.
    * T: Write in output.path/CDR3_length separate plot the V usage, J usage and CDR3 motifs
     for all CDR3 length. Only applicable if return.object != 2 and plot=T.

 - plot.cdr3.norm (default=0)
    * 0: Show the CDR3 motifs of the input and of the baseline repertoire or input2, possibly | P(VJ).
    * 1: Show the CDR3 motifs of the input TCRs after subtracting the baseline repertoire (not recommanded).
    * 2: Show the CDR3 motifs of the input TCRs after normalising by the baseline repertoire (motif of normalised fold-change, not recommanded)

 - plot.sd (default=T)
    * T: Show standard deviation for P(V), P(J) and P(L), if such data are provided for baseline and/or input1
    * F: Do not show standard deviation

 - plot.VJ.switch (default=1)
    * 1: Show the VJ usage as a scatter plot with inner colors, outer colors and shapes of the points based
       on the V/J gene names (see lookup table in MixTCRviz/figures/, scheme1).
    * 1.2: Show the VJ usage as a scatter plot with colors of the points based
       on the V/J gene names (see lookup table in MixTCRviz/figures/, scheme2).
    * 1.3: Show the VJ usage as a scatter plot with black points and V/J gene label in colors (see lookup table in MixTCRviz/figures/, scheme2).
    * 2: Show the VJ usage as a bar plot (see lookup table in MixTCRviz/figures/, scheme2).

 - plot.modelsCombined (default=FALSE)
    * FALSE or empty string: Show the data for each model separately.
    * TRUE or a string: Show the data for all models combined in a single
        figure (based on the format plot.VJ.switch=2). When given as a
        logical, the default name "modCombined" is used. Otherwise, the string
        is used for the resulting figure filename.

 - label.neg (default=FALSE): If TRUE, show also the labels of the genes most depleted in input1

 - label.diag (default=0.3): Decide to keep some label along the diagonal inthe upper corner, based on label.diag value. This can be useful when comparing two epitope-specific TCRs

 - label.min.fr (default=c(0.05,0.05)): Region (i.e., X - Y rectangle) of the left corner of V/J plots with no gene label. If only one number is provided, it will be used on both the X and Y axes.

 - keep.incomplete.chain (default=T): If False, incomplete chains are discarded.
     Even if input only consists of complete chains, incomplete chain can occur when one V/J gene cannot be corrected,
     or when there is some incompatibilities between V/J names and CDR3 sequences

 - chain.list.output (default="AB")
    * A: Only the alpha chain is plotted in output;
    * B: only the beta chain is plotted in output;
    * AB both chains are plotted in output

 - input1.name (default="Input"): Provide a generic name for
   the input TCRs in the plots (e.g., Epitope Specific). Avoid names with more
   than 20 characters

 - input2.name (default=NULL): If a second set of TCRs is provided
   (i.e., input2 != NULL), Provide a generic name for the input2 TCRs in the
   plots. Avoid names with more than 20 characters.

 - output.format (default="pdf"): Choose the format for the output
      plots (can be "pdf", "png" or "jpg").

 - interactive.plots (default=F):
    * F: Do not create an html file with interactive plots.
    * T: Create an html file with interactive plots. Only applicable if return.object != 2 and plot=T

 - print.size (default=TRUE): If TRUE, print the number of TCRs in input1 in the plots.

 - plot.title (default=TRUE): If TRUE, print the model name as title to the plots.

 - set.title (default=NULL): Set the title of the plots. If empty, model names in input1 are used as title.

 - build.clones (default=F): If TRUE and if the data are provided in format with clone.id, reconstruct the actual clones in output

 - keep.colnames.origin (default=F): If TRUE, keep the name of the columns dedicated the TCRs in the input. The option is not supported for format with clone.id

 - return.object (default=0):
    * 0: do not return any object.
    * 1: return an object with plots ($plot, if plot=T), statistics ($stat) and processed.data ($processed.data, if input1 is not a precomputed list with the stat) for each model
    * 2: Only return an object and do not write anything to file. In this case output.path can be omitted.
    Plots can be displayed in R studio, for instance. Quality will vary depending on your R studio settings, so this option is not recommended for high-quality figures.

# OUTPUT


MixTCRviz creates a directory (output.path). The output.path/ directory contains the motifs (e.g. pdf files) for each model.

- If output.stat==T, the output.path/stats/ contains .rds files with all the stats for each model.
- If output.processed.data==T, the output.path/processed_data/ contains .csv files with the actual data used to build the motifs.
- If plot.logo.length==1, the output.path/CDR3_length/ directory shows the V/J usage and CDR3 motifs for multiple lengths for both chains.
- If return.object==1 or 2, the motifs, stats and processed data are returns as lists

# Data format

By default, MixTCRviz uses column names c("TRAV","TRAJ","cdr3_TRA","TRBV","TRBJ","cdr3_TRB") to define a TCR,
"species" to indicate the species and "model" to define groups of TCRs (e.g., binding to the same epitope).
- For single-chain data, only one chain can be provided. In those cases, it is recommanded to define the chain in chain.output.list="A" or "B".
- "model" can be skipped, in which case all TCRs will be analyzed together and the output file will take the name given in model.default (default="Model_default").
- "species" can be skipped, in which case all TCRs are assumed to come from the same species given in species.default (default="HomoSapiens"). If you have data coming from multiple species, you need to have the "species" column.
- Other column names are supported, including "Va", "V_alpha", "CDR3a", "CDR3A", "CDR3_alpha", etc. for data with both chain.
 Or "V", "v_gene","V-region","aaSeqCDR3","CDR3",etc for single chain data, see list in data_raw/TidyVJ/mapping_colnames.csv for a full description

Other supported formats treating each chain in a different row include:

 - VDJdb with the columns: c("V", "J", "CDR3")
 - 10X Genomics format with columns: c("v_gene", "j_gene", "cdr3")
 - Qiagen with the columns: c("V-region", "J-region", "CDR3 amino acid seq")
 - Adaptive Biotech with the columns: c("vGeneName", "jGeneName", "aminoAcid")
 - Adaptive Biotech v4 with the columns: c("v_resolved", "j_resolved", "amino_acid")
 - AIRR with the columns: c("v_call", "j_call", "junction_aa")
 - MiXCR with the columns: c("allVHitsWithScore", "allJHitsWithScore", "aaSeqCDR3")

By default, both chains are treated independently, without reconstructing clones.
To reconstruct alpha-beta clones, you can use the build.clones=T.
However, you need to have *exactly* one colum indicating the clone_id labelled as "clone_id", "cell_id", "cloneId", "barcode" or "complex.id"





# OTHER INFORMATION

* If working with epitope-specific TCRs, we encourage users to define model names which capture both the epitope sequence and the MHC restriction. Using only epitope sequences as "model" is possible, but can lead to issues when the same epitope is restricted to different MHC.

* V/J genes are key to the TCR binding motifs in MixTCRviz and only V/J names compatible with the IMGT nomenclature can be considered. Even if correct.gene.names==1 allows to correct several wrong V/J names, we strongly encourage the users to use only V/J gene names compatible with IMGT

* In the default setting, error bars on the baseline distributions of V/J segments represent the variability observed across multiple studies, encompassing different sequencing protocols, different centers and different donors.

* PCR / Sequencing / TCR reconstruction errors frequently occur in TCR-Seq data. Although the option check.cdr3.mode = 1, can detect some of these errors, we encourage users to carefully check the quality of their CDR3 sequences.

* As with all motif visualisation tools, limited numbers of TCRs have a big impact on the interpretation of the results.
Therefore, it is hilghly recommanded to use sets of TCRs with enough sequences to be able to interpret frequencies plotted in MixTCRviz.

* When comparing to baseline TCR repertoire, we encourage to use renormVJ=1,
so that the comparisons of CDR3 length distributions and motifs is not confounded by the specific V/J usage in input1 (i.e., baseline shows the expected distributions and motifs knowing P(VJ) in input1).
This option is often less relevant when comparing two TCR datasets (input1 and input2), this is why renormVJ is by default put to 0 unless renormVJ = 1.


* The information about the species is important to choose an appropriate baseline.
In case you are using the same model (e.g.,MHC_peptide) with both human and mouse TCRs, we recommend distinguishing both cases in the model column
(e.g., A0201_LLWNGPMAV_HomoSapiens and A0201_LLWNGPMAV_MusMusculus) and indicating the species in the appropriate column.

* When using use.allele=T, it is important to realize that presence/absence of specific alleles in a given sample also reflects the genetic background of the donor.
For these reasons, some alleles may appear to be enriched in the input TCRs versus default baseline repertoires, but this enrichment may not be linked to any signal of specificity (e.g., epitope specificity).
In addition, determining the correct allele from TCR-Seq data can be challenging, and sequencing errors can easily result in wrong allele calls.
We therefore recommend analyzing data at the gene level (use.allele=F) to avoid confounding factors related to genetic background/TCR reconstruction issues.
Alternatively, the baseline TCR repertoire can also be sequenced in each patient, and used as input for MixTCRviz (baseline parameter).
