
###############
# MixTCRviz is a tool to visualize important properties of a set of TCRs (e.g., epitope-specific TCRs).
#
# MixTCRviz can be used freely by academic groups for non-commercial purposes (see license).
# The product is provided free of charge, and, therefore, on an "as is"
#  basis, without warranty of any kind.
#
# FOR-PROFIT USERS
# If you plan to use MixTCRviz or any data provided with the script in any for-profit
# application, you are required to obtain a separate  license.
# To do so, please contact Nadette Bulgin (nbulgin@lcr.org) at the Ludwig Institute for Cancer Research Ltd.
#
# If you use MixTCRviz in a publication, please cite:
# xxx
#
# For scientific questions, please contact David Gfeller (David.Gfeller@unil.ch)
#
# Copyright (2024) David Gfeller
###############


############
INTRODUCTION
############

MixTCRviz is a tool to visualize important properties of a set of TCRs (e.g., epitope-specific TCRs).
It focuses on V usage, J usage, CDR3 length distribution and CDR3 motifs, for both alpha and beta chain.

By default, properties of the input TCRs are compared to those of a "baseline" TCR repertoire (mouse or human).
Alternatively, users can choose to compare with another set of TCRs ("reference").

############
INSTALLATION
############

There are different ways of "installing" MixTCRviz:
  1) It can be directly installed from the GitHub page (you need to have
    credentials correctly set so that Rstudio can access the private repository):
    ``` r
    # install.packages("devtools") # Only needed if you don't already have devtools
    devtools::install_github("GfellerLab/MixTCRviz")
    ```
  2) You can download the MixTCRviz directory from the GitHub page and then open
    Rstudio setting its working directory as MixTCRviz folder. Then you can
    compile it and install it:
    ``` r
    devtools::build()
    install.packages("../MixTCRviz_0.0.1.tar.gz", repos=NULL)
    ```

  3) If you don't need to use it through Rscript and are in process of developing
    the package, it would be tiring to always build it, install it, reload it,
    etc. Instead, you can set the working directory to the MixTCRviz folder and
    then simply:
    ``` r
    devtools::load_all()
    # Could also give the path to MixTCRviz folder if not in the working directory
    ```
    You can reuse this command everytime changes are made to the code so that
    the updates are directly available in the R environment. (You need to reuse
    this load_all in every new R session to use this in dev code instead of the
    installed package). And you should of course not do a "library(MixTCRviz)"
    if using the load_all, this command makes as if you had done such one (the
    library... would probably load the built+installed package not the
    development version).

You may be prompted to install several packages.

Testing the INSTALLATION:

In the MixTCRviz directory run:

Rscript test_MixTCRviz.R

Alternatively, you can run the code test_MixTCRviz.R in Rstudio (or any R interface)

The output in test/out should be the same as in test/out_compare

############
RUNNING
############

MixTCRviz should be primarily run in R, by calling the MixTCRviz function. To
this end, in R, you need to call the MixTCRviz function with the mandatory
parameters (e.g., MixTCRviz::MixTCRviz(input1=input, output.path=output)).
(note that the "MixTCRviz::" isn't needed if using the "load_all"
version from above or if you used library(MixTCRviz)).


Input parameters of MixTCRviz:

Mandatory parameters:

- input1: csv files or data.frame with the input TCRs.
    Columns should ideally include: TRAV, TRAJ, cdr3_TRA, TRBV, TRBJ, cdr3_TRB, model, species
        - The "TRAV", "TRAJ", "TRBV", "TRBJ" should follow the IMGT nomenclature, with or without allele (see below for potential name correction).
            If a column is missing, empty values will be used.
        - The "cdr3_TRA" and "cdr3_TRB" columns should provide CDR3A/CDR3B sequences, following the standard definition (e.g., CAVNSDGQKLLF).
            Cases with non-amino acid characters, or length <7 or >22 will be not be considered (i.e., put to NA).
            If a column is missing, empty values will be used.
        - The "model" column typically describes the epitopes/experiments/classes/... Each model will be treated independently in MixTCRviz.
             If missing, all TCRs will be considered as coming from the same model (Model_default).
        - The "species" column indicates the species of the TCR (e.g., NOT the source organism of the epitope). It should be "HomoSapiens" or "MusMusculus"
             If not provided, "HomoSapiens" is taken as default.

- output.path: name of the output directory (if not already existing, it will be created).

Optional parameters:

- input2 (default=""): csv file or data.frame containing the "reference" TCRs to be used in comparisons. Same format as input1
    In particular, the set of models in the "model" field needs to be the same as in input1
    so that comparisons are performed for each model separately.
    If no "model" is given, all TCRs will be considered as coming from the same model ("Model_default").
    If input2 is provided, the comparisons is performed with this reference, and not the baseline
    repertoire.

- baseline.file (default=""): .rds file containing information about baseline repertoire.
    If empty, the default baseline repertoires are used.

- use.allele (default=0): 0: All V/J alleles are merged at the gene level (recommended).
    1: Alleles are kept, including mouse TRAV genes from different strains.
    Currently, use.allele is fixed to 0 for mouse TCRs.

- use.mouse.strain (default=0): 0: Merge the different TRAV genes corresponding to different mouse strains (e.g., TRAV10D, TRAVN10 and TRAV10 become all TRAV10).
    This is recommended since differences between input TCRs and baseline TCRs may result from the use of different mouse strains, and not be related to specificity in input TCRs.
    In addition different TCR reconstruction tools use different approaches to call these genes, which can create artificial enrichment in one of them versus the baseline.
    1: The different TRAV segments (e.g., TRAV10D, TRAVN10 and TRAV10) are treated separately.

- renormVJ (default=1): 1: When comparing to the baseline TCR repertoire,
    build CDR3 length distributions and CDR3 motifs corresponding to the V-J usage observed in the input TCRs.
    In the plots the mark " | VJ" is used to indicate that CDR3 length distributions and motifs correspond to those expected with the V-J usage in the input TCRs.
    0: Use CDR3 length distribution and motif from the full repertoire (this is the only option when input2 is provided for comparison with a user-given reference)

- N.min (default=10): Minimum number of TCR (i.e., V-J-CDR3) for at least one chain.

- correct.gene.names (default=1): 1: Attempt to correct V/J gene names not in IMGT based on internal map. Put to NA genes that could not be corrected
    0: Do not attempt to correct V/J gene names. Put to NA genes not in IMGT.

- plot (default=1): 1: Plot the data in output.path/plots/ and create .rds object with the statistics for each model in output.path/stats/.
    0: only create .rds object with the statistics

- plot.cdr12.motif (default=0): 0: Only show motifs for CDR3. 1: Include sequence motifs for CDR1 and CDR2.

- plot.oneline (default=0): 0: Show the data on two lines (better for clarity). 1: Show all plots in a single line (can be useful to compare different models).

- plot.logo.length (default=0): 0: Show only the CDR3 motifs for the most frequent CDR3 length.
    1: Show in a separate plot the V usage, J usage and CDR3 motifs for all CDR3 length, for both alpha and beta chains.

- plot.cdr3.norm (default=0):
    0: Show the CDR3 motifs of the baseline repertoire.
    1: Show the CDR3 motifs of the input TCRs after subtracting the baseline repertoire.
    2: Show the CDR3 motifs of the input TCRs after normalising by the baseline repertoire (motif of normalised fold-change).

- chain.list.output (default="AB"); A: Only the alpha chain is plotted in output; B: only the beta chain is plotted in output; AB both chains are plotted in output

- output.format (default="pdf"); Choose the format for the output plots.

- input1.name (default="Epitope Specific"); Provide a generic name for the input TCRs in the plots (e.g., Epitope Specific).
     Avoid names with more than 20 characters
 - input2.name (default="Reference"); If a second set of TCRs is provided (i.e., input2 != ""), Provide a generic name for the input2 TCRs in the plots.
     Avoid names with more than 20 characters

#############
OUTPUT
#############

MixTCRviz creates a directory (output.path). The output.path/plots/ directory contains the plots for each model.
The output.path/stats/ contains .rds objects summarizing the data for each model.
If plot.logo.length==1, the output.path/plots/CDR3_length/ directory shows the V/J usage and CDR3 motifs for multiple lengths for both chains.


#############
OTHER INFORMATION
#############

** When comparing to baseline TCR repertoire, we encourage to use renormVJ=1,
so that the comparisons of CDR3 length distributions and motifs are not biased by the V/J usage in the input TCRs.
With renormVJ=1, the baseline CDR3 length distributions and motifs will be different for different models.
With renormVJ=0, some of the differences between the model-specific and the baseline CDR3 length distributions and motifs
will be redundant with the differences observed in V/J usage.

** For a given model, MixTCRviz requires at least 10 TCR (V-J-CDR3) for at least one chain. Missing data (e.g., only one chain) will result in empty plots.

** By default, the minimum number of sequences (i.e. V-J-CDR3) for at least one chain is 10.
Using less data is possible by changing the N.min parameter. However the plots of MixTCRviz will no longer be much meaningful with too little data to perform meanigful statistics.

** Although we have compiled several frequent issues with V/J gene names, as well as a way to correct some of them,
we strongly recommend using V/J genes from IMGT. If other V/J gene names are used, it is likely that several entries will not be considered by MixTCRviz.

** The information about the species is important to choose an appropriate baseline.
In case you are using the same model (e.g.,MHC_peptide) with both human and mouse TCRs, we recommend distinguishing both cases in the model column
(e.g., A0201_LLWNGPMAV_HomoSapiens and A0201_LLWNGPMAV_MusMusculus) and indicating the species in the appropriate column.

** When chosing plot.logo.length=1, the baseline/reference V and J usage (x-axis) correspond to those of TCRs with the given lengths.
Hence the exact values change for different lengths.

** When usinge use.allele=1, it is important to realize that presence/absence of specific alleles in a given sample also reflects the genetic background of the donor.
For these reasons, some alleles may appear to be enriched in the input TCRs versus default baseline repertoires, but this enrichment may not be linked to any signal of specificity (e.g., epitope specificity).
In addition, determining the correct allele from TCR-Seq data can be challenging, and sequencing errors can easily result in wrong allele calls.
We therefore recommend analyzing data at the gene level (use.allele=0) to avoid confounding factors related to genetic background/TCR reconstruction issues.
Alternatively, the baseline TCR repertoire can also be sequenced in each patient, and used as input for MixTCRviz (baseline.file parameter).

** The .rds objects in the output.path/stats/ directory contain a summary of the statistics for each model for each chain.
This is provided as a list with the following entries:
- "L": list of CDR3 length
- "countL": Number of TCRs for each CDR3 length
- "countV": Number of TCRs with each V segment
- "countJ": Number of TCRs with each J segment
- "countV.L": Number of TCRs with each V segment for each CDR3 length
- "countJ.L": Number of TCRs with each J segment for each CDR3 length
- "countCDR1": Amino acid count at each positions of the CDR1 loops
- "countCDR2": Amino acid count at each positions of the CDR2 loop
- "countCDR3.L": Amino acid count at each positions of the CDR3 loop, for each CDR3 length
- "countVJ": Number of TCRs with each V-J combination
- "countVJ.L": Number of TCRs with each V-J combination, for each CDR3 length
