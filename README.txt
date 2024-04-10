
###############
# MixTCRviz is a tool to visualise and compare the properties of repertoires of epitope specific TCRs.
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

MixTCRviz is a tool to visualize important properties of epitope-specific TCRs from human or mouse.
It focuses on V usage, J usage, CDR3 length distribution and CDR3 motif of the most frequent CDR3 length, for both alpha and beta chain.

By default, properties of epitope specific TCRs are compared to those of a "baseline" TCR repertoire (mouse or human).
Alternatively, users can choose to compare with another set of TCRs ("reference").


############
INSTALLATION
############

- Download and open MixTCRviz.zip
- In MixTCRpred.R, provide the path to MixTCRviz: MixTCRviz.path <- "PATH_TO_MixTCRviz"


Testing the INSTALLATION:

In the MixTCRviz directory run:

Rscript test_MixTCRviz.R

The output in test/out_compare should be the same as in test/out

############
RUNNING
############

MixTCRviz should be primarily run in R, by calling the MixTCRviz function. To this end, in R, you need to:

1) Source the MixTCRviz.R (e.g., source("PATH_TO_MixTCRviz/MixTCRviz.R"))
2) Call the MixTCRviz function with the mandatory parameters (e.g., MixTCRviz(input.file1=input, output.path=output))

You may be prompted to install several packages.

Input parameters of MixTCRviz:

Mandatory parameters:

- input.file1: csv files with the epitope-specific TCRs.
    Columns MUST include: TRAV, TRAJ, cdr3_TRA, TRBV, TRBJ, cdr3_TRB, model, species
        - The "model" column describes the epitope. We recommand using MHC_peptide (A0201_LLWNGPMAV), although other names are allowed.
        - The "species" column should be "HomoSapiens" or "MusMusculus".
        - The "TRAV", "TRAJ", "TRBV", "TRBJ" should follow the IMGT nomenclature, with or without allele (see below for potential name correction)
        - The "cdr3_TRA" and "cdr3_TRB" columns should provide CDR3A/CDR3B sequences, following the standard definition (e.g., CAVNSDGQKLLF).

- output.path: name of the output directory (if not existing, will be created).

Optional parameters:

- input.file2 (default=""): File containing the "reference" TCRs to be used in comparisons. Same format as input.file1
    If provided, the comparisons is performed with this reference, and not the baseline repertoire.

- baseline.file (default=""): .rds file containing information about baseline repertoire.
    If empty, the default baseline repertoires are used.

- use.allele (default=0): 0: All V/J alleles are merged at the gene level (recommanded). 1: Alleles are kept.

- renormVJ (default=1): 1: When comparing to the baseline TCR repertoire,
    build CDR3 length distributions and CDR3 motifs corresponding to the V-J usage observed in epitope specific TCRs.
    In the plots the mark " | VJ" is used to indicate that CDR3 length distributions and motifs correspond to those expected wiht the V-J usage in the epitope-specific TCRs.
    0: Use CDR3 length distribution and motif from the full repertoire (this is the only option when input.file2 is provided for comparison with a user-given reference)

- N.min (default=10): Minimum number of TCR (i.e., V-J-CDR3) for at least one chain.

- correct.gene.names (default=1): 1: Attempt to correct V/J gene names not in IMGT based on internal map. Put to NA genes that could not be corrected
    0: Do not attempt to correct V/J gene names. Put to NA genes not in IMGT.

- plot (default=1): 1: Plot the data in output.path/plots/ and create .rds object with the statistics for each model in output.path/stats/.
    0: only create .rds object with the statistics

- plot.cdr12.motif (default=0): 0: Only show motifs for CDR3. 1: Include sequence motifs for CDR1 and CDR2.

- plot.oneline (default=0): 0: Show the data on two lines (better for clarity). 1: Show all plots in a single line (can be useful to compare different epitopes).

- plot.logo.length (default=0): 0: Show only the CDR3 motifs for the most frequent CDR3 length.
    1: Show in a separate plot the V usage, J usage and CDR3 motifs for all CDR3 length, for both alpha and beta chains.

#############
OUTPUT
#############

MixTCRviz creates a directory (output.path). The output.path/plots/ directory contains the plots for each model (e.g. epitope).
The output.path/stats/ contains .rds objects summarizing the data for each model (e.g., epitope, see below).
If plot.logo.length==1, the output.path/plots/CDR3_length/ directory shows the V/J usage and CDR3 motifs for multiple lengths.


#############
OTHER INFORMATION
#############

** When comparing to baseline TCR repertoire, we encourage to use renormVJ=1,
so that the comparisons of CDR3 length distributions and motifs are not biased by the V/J usage in epitope-specific TCRs.
With renormVJ=1, the baseline CDR3 length distributions and motifs will be different for different epitopes.
With renormVJ=0, some of the differences between the eptiope-specific and the baseline CDR3 length distributions and motifs
will be redundant with the differences observed in V/J usage.

** For a given epitope, MixTCRviz requires at least 10 TCR (V-J-CDR3) for at least one chain. Missing data (e.g., only one chain) will result in empty plots.

** By default, the minimum number of sequences (i.e. V-J-CDR3) for at least one chain is 10.
Using less data is possible by changing the N.min parameter. However the plots of MixTCRviz will not longer be much meaningful with too little data to perform the statistics.

** Although we have compiled several frequent issues with V/J gene names, as well as a way to correct some of them,
we strongly recommand using V/J genes from IMGT. If other V/J gene names are used, it is likely that several entries will not be considered by MixTCRviz.

** The information about the species is important to choose an appropriate baseline.
In case you are using the same model (e.g.,MHC_peptide) with both human and mouse TCRs, we recommand distinguishing both cases in the model name
(e.g., A0201_LLWNGPMAV_HomoSapiens and A0201_LLWNGPMAV_MusMusculus) and indicating the species in the appropriate column.

** When chosing plot.logo.length=1, the baseline/reference V and J usage (x-axis) correspond to those of TCRs with the given lengths.
Hence the exact values change for different lengths.

** When usinge use.allele=1, it is important to realize that presence/absence of specific alleles in a given patient also reflects the genetic background of this patient.
For these reasons, some alleles may appear to be enriched in epitope specific TCRs vesus default baseline repertoires, but this enrichment is not linked to any signal of epitope specificity.
In addition, determining the correct alleles from TCR-Seq data can be challeging, and sequencing errors can easily result in wrong allele calls.
We therefore recommand analysing data at the gene level to avoid confounding factors related to genetic background/TCR reconstruction issues.
Alternatively, the baseline TCR repertoire can also be sequenced in each patient, and used as input for MixTCRviz (baseline.file parameter).

** The .rds objects in the output.path/stats/ directory contain a summary of the statistics for each model for each chain.
This is provided as a list with the following entries:
- "L": list of CDR3 lenght
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
