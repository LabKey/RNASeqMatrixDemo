RNASeqMatrixDemo
================

Demo of LabKey Server's module file-based pipeline using RNASeq data.  This simple pipeline reformats an existing matrix into one appropriate for importing into LabKey, performs further analysis and generates images, and imports the analyis results.

Prerequisites:
---------------

1. Install LabKey Server (14.3 is required) and ensure the Microarray module is included.

2. Install this module by zipping up the root RNASeqMatrixDemo directory and copying to the server's "modules" directory.

3. Install R and the "edgeR" bioconductor package.
  - source("http://www.bioconductor.org/biocLite.R")
  - biocLite("edgeR")

4. Use the "data/GSE56845_counts.txt" RNASeq matrix or obtain the original "GSE56845_gene_counts_Rhesus_Ensembl.txt" from Geo:
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56845

> Barrenas F, Palermo RE, Agricola B, Agy MB, Aicher L, Carter V, Flanary L,
Green RR, McLain R, Li Q, Lu W, Murnane R, Peng X, Thomas MJ, Weiss JM, Anderson 
DM, Katze MG. Deep transcriptional sequencing of mucosal challenge compartment
from rhesus macaques acutely infected with simian immunodeficiency virus
implicates loss of cell adhesion preceding immune activation. J Virol. 2014 Jul
15;88(14):7962-72. doi: 10.1128/JVI.00543-14. Epub 2014 May 7. PubMed PMID:
24807713.



Setup
-----

1. Start the server.
2. Create a new assay folder, enable the RNASeqMatrixDemo module.
3. Upload the "GSE56845_counts.txt" matrix file to the folder.
3. Create a new ExpressionMatrix assay with the default fields.
4. Add the "Feature Annotation Sets" webpart to the portal.
5. Import the "illumina-feature-set.tsv" found in the "data" directory of this module.
4. Create a new "General" assay named "DifferentialExpression" and with the following fields:
  * Run Fields:
    - multiDimensionalScalingPlot (label=MDS, type=File)
    - meanVariancePlot (label=Variance, type=File)
  * Data Fields:
    - FeatureId (Text)
    - logConc (Number)
    - logFC (Number)
    - pvalue (Number)
    - padj (Number)


Running the Pipeline
--------------------

From the file-browser, select the matrix file then click "Import Data".

<img src="https://raw.githubusercontent.com/LabKey/RNASeqMatrixDemo/master/docs/img/import-data.png"/>

Next, enter a name for the run and select the target ExpressionMatrix assay and the feature annotation set you created during the setup.  By default, the "Import Values" checkbox is unchecked and the matrix will not be imported into the assay -- only pointers to the original data file and the output files will be captured as part of the run.

<img src="https://raw.githubusercontent.com/LabKey/RNASeqMatrixDemo/master/docs/img/create-matrix-form.png"/>

As the pipeline executes, files are written into the pipeline directory under a "RNASeqMatrixDemo" directory.  The pipeline will create two assay runs: one ExpressionMatrix assay run and one DifferentialExpression assay run.  As a part of the differential expression analysis, a series of images will also be generated.  Further analysis could be performed on the values imported into the DifferentialExpression assay.

The original input file and samples in the filtered matrix will be attached to an experiment run for the ExpressionMatrix assay import:

<img src="https://raw.githubusercontent.com/LabKey/RNASeqMatrixDemo/master/docs/img/matrix-exp-run.png" />

The filtered matrix file and the output files will be attached to an experiment run for the DifferentialExpression assay import:

<img src="https://raw.githubusercontent.com/LabKey/RNASeqMatrixDemo/master/docs/img/diff-expr-exp-run.png"/>

The output file values will be imported into the assay's results table:

<img src="https://raw.githubusercontent.com/LabKey/RNASeqMatrixDemo/master/docs/img/diff-expr-assay-run.png"/>
