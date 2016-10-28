##
#  Copyright (c) 2014-2016 LabKey Corporation
# 
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
##
# Add custom rpackage location to our search path
.libPaths(c("/Users/kevink/RNASeqDemo/rnaseq-rpackages", .libPaths()))

library(edgeR)

# read the filtered matrix
counts <- read.csv(file="${input.counts.tsv}", header = TRUE, sep="\t", row.names=1)

# Quality Control
# ---------------
# We can look at the sum of each column and observe its library size. Next we can look for those genes that have very low gene counts because this will effect our analysis downstream. Once we see how many low expressed genes we have we will filter out those genes that have too few reads. We will keep genes that have at least 1 read per million in at least 3 samples

colSums( counts ) # Library Sizes

colSums( counts ) / 1e06 # Library Sizes in millions of reads

table( rowSums( counts ) )[ 1:30 ] # Number of genes with low counts
head( counts )
group <- c(rep("Baseline", 3) , rep("Infected", 3))
cds <- DGEList( counts , group = group )
names( cds )
head(cds$counts) # original count matrix
cds$samples # contains a summary of your samples
sum( cds$all.zeros ) # How many genes have 0 counts across all samples

# Normalization
# -------------
# Here we will normalize the data by calculating normalization factors and effective libraries sizes. Normalization factors correct for the differences within the samples. Normalization factors are then multipled by the e???ective library size and the product becomes the library size which are displayed in the table below.


cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
dim( cds )
cds <- calcNormFactors( cds )
cds$samples

# effective library sizes
library_sizes <- cds$samples$lib.size * cds$samples$norm.factors


# Dimensionality Reduction and Biological Variance
# -------------------------------------------------
# 
# Lets make a MDS plot(Multi-Dimensional Scaling plot, like a PCA) this will take all my now normalized counts and plot them in geometric space so we can look at the relationship between the conditions, how replicates bind together, and the level of variability that exists between samples and conditions.

png(filename="${output-mds.png}", width=800, height=800)
plotMDS( cds , main = "MDS Plot for Count Data", labels = colnames( cds$counts ) )
dev.off()

# Next lets look at the Biological coefficient of variation. The reason why we are looking at this is 
# we will be applying a negative binomial distribution model to transcript counts to look at whats changing across our conditions. We know we will allow for some gene-specific variability when we feed this model and some genes have more biological variability than others. The BCV plot shows how much the variance from the counts exceeds the variance that would arise from a Poisson model.

#plotBCV(cds, cex=0.4, main="Biological coefficient of variation (BCV) vs abundance")
## Error: No dispersions to plot

# Estimating Dispersions
# ----------------------
# [Taken from the edger tutorial :http://cgrlucb.wikispaces.com/file/view/edgeR_Tutorial.pdf]By looking at the variance of the negative binomial model and possion we can estimate the common dispersion of our gene counts. Once the common dispersion is estimated we can estimate the tagwise dispersions. In this scenario, each gene will get its own unique dispersion estimate, but the common dispersion is still used in the calculation. The tagwise dispersions are squeezed toward the common value. The amount of squeezing is governed by the paramter prior.n. The higher prior.n, the closer the estimates will be to the common dispersion. The recommended value is the nearest integer to 50/(#samples ??? #groups). For this data set that???s 50/(7 ??? 2) = 10.

cds <- estimateCommonDisp( cds )
names( cds )
# The estimate
cds$common.dispersion

sqrt( 200 ) # poisson sd
sqrt( 200 + 200^2 * cds$common.dispersion ) # negative binomial sd
sqrt( 200 + 200^2 * cds$common.dispersion ) / sqrt( 200 ) #

cds <- estimateTagwiseDisp( cds , prior.df = 10 )
names( cds )
summary( cds$tagwise.dispersion )


# More shrinkage/sqeezing toward the common
cds <- estimateTagwiseDisp( cds , prior.df = 25 )
summary( cds$tagwise.dispersion ) # not much changed, but the ends got squeezed in quite a bit.

# The recommended setting for this data set is the default of 10. Let???s stick with that.
cds <- estimateTagwiseDisp( cds , prior.df = 10 )

# Mean Variance Plot
# ------------------
# This is basically test we perform on dispersion data. The dispersion parameters will show how well the data  ???ts in the mean-variance relationship which is important because that is what the negative binomial is based off of. Four things are shown in the plot: the raw variances of the counts (grey dots), the variances using the tagwise dispersions (light blue dots), the variances using the common dispersion (solid blue line), and the variance = mean a.k.a.poisson variance (solid black line). 

png(filename="${output-var.png}", width=800, height=800)
plotMeanVar( cds , show.raw.vars=TRUE ,
    show.tagwise.vars=TRUE ,
    show.binned.common.disp.vars=FALSE ,
    show.ave.raw.vars=FALSE ,
    NBline = TRUE ,
    nbins = 100 ,
    pch = 16 ,
    xlab ="Mean Expression (Log10 Scale)" ,
    ylab = "Variance (Log10 Scale)" ,
    main = "Mean-Variance Plot" )
dev.off()


# Differential Expression
# -----------------------
# Now we can calculate our differential expressed genes that we can feed into our functional analysis. This will allow us the ability to answer biological questions from the data by comparing the infected and baseline conditions to each other.To start this process we first run a series of pair-wise tests for di???erential expression between two groups using the extact tests function. Afterwards we will run Toptags function which will generate adjusted p values with our DE results.


de.cmn <- exactTest( cds , dispersion="common" , pair = c("Baseline" , "Infected" ) )
de.tgw <- exactTest( cds , dispersion="common" , pair = c("Baseline" , "Infected" ) )
de.poi <- exactTest( cds , dispersion = 1e-06 , pair = c("Baseline" , "Infected" ) )
names( de.tgw )
de.tgw$comparison # which groups have been compared
head( de.tgw$table ) # results table in order of your count matrix.
head( cds$counts )

options( digits = 3 ) # print only 3 digits
topTags( de.tgw , n = 20 , sort.by = "p.value" ) # top 20 DE genes
# Back to count matrix for tagwise analysis
cds$counts[ rownames( topTags( de.tgw , n = 15 )$table ) , ]

# Sort tagwise results by Fold-Change instead of p-value
resultsByFC.tgw <- topTags( de.tgw , n = nrow( de.tgw$table ) , sort.by = "logFC" )$table
head( resultsByFC.tgw )

# Store full topTags results table
resultsTbl.cmn <- topTags( de.cmn , n = nrow( de.cmn$table ) )$table
resultsTbl.tgw <- topTags( de.tgw , n = nrow( de.tgw$table ) )$table
resultsTbl.poi <- topTags( de.poi , n = nrow( de.poi$table ) )$table
head( resultsTbl.tgw )


##
# Output Results
#

# Change column names to be specific to the analysis, logConc and logFC are the same in both.
#colnames( resultsTbl.cmn ) <- c( "logConc" , "logFC" , "pVal.Cmn" , "adj.pVal.Cmn" )
#colnames( resultsTbl.tgw ) <- c( "logConc" , "logFC" , "pVal.Tgw" , "adj.pVal.Tgw" )
colnames( resultsTbl.tgw ) <- c( "logConc" , "logFC" , "pvalue" , "padj" )

# Below provides the info to re-order the count matrix to be in line with the order of the results.
wh.rows.tgw <- match( rownames( resultsTbl.tgw ) , rownames( cds$counts ) )
wh.rows.cmn <- match( rownames( resultsTbl.cmn ) , rownames( cds$counts ) )
head( wh.rows.tgw )

# change the column header of 
# write our final DE results to file below

# add the row names as a column in the dataset
feature_id <- rownames(resultsTbl.tgw)
res <- cbind(feature_id, resultsTbl.tgw)
colnames(res)[1] <- "FeatureId"

write.table(res, file = "${output.tsv}", sep = "\t", row.names=FALSE);

# read the task info
taskInfo <- read.table("${pipeline, taskInfo}",
                       col.names=c("name", "value", "type"),
                       header=FALSE, check.names=FALSE,
                       stringsAsFactors=FALSE, sep="\t", quote="",
                       fill=TRUE, na.strings="")

# write out run properties to attach the generated images to the run
outputParams <- data.frame(name=c("assay run property, multiDimensionalScalingPlot",
                                  "assay run property, meanVariancePlot"),
                           value=c(file.path(taskInfo$value[taskInfo$name == "analysisDirectory"], "${output-mds.png}"),
                                   file.path(taskInfo$value[taskInfo$name == "analysisDirectory"], "${output-var.png}")))

write.table(outputParams, file = "${pipeline, taskOutputParams}", sep = "\t", qmethod="double", col.names=TRUE, row.names=FALSE)

