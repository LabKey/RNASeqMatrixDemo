# Add custom rpackage location to our search path
.libPaths(c("/Users/kevink/RNASeqDemo/rnaseq-rpackages", .libPaths()))

library(edgeR)

# read the raw matrix
raw.data <- read.csv(file="${input.txt}", header = TRUE, sep="\t")


#we check to make sure our data loaded sucessfully
dim(raw.data)


#lets take a peek at the data set
head(raw.data)


#I am only interested in the Baseline samples and Infected samples at later time points
#I will select only those columns I am interested in.

counts <- raw.data[,c(3,4,5,15,17,19) ]

#I want to keep the gene names because I know I will need them later for analysis so I will store them
rownames( counts ) <- raw.data[ , 1 ] # gene names

#I want to simply the column names to either baseline and Infected
colnames( counts ) <- paste(c(rep("Baseline",3),rep("Infected",3)),c(1:3,1:3),sep="") # sample names

#Lets see how the data has changed 
dim( counts )

# write out the progress so far in a format that can be imported into LabKey
# add the row names as a column in the dataset
feature_id <- rownames(counts)
res <- cbind(feature_id, counts)
colnames(res)[1] <- "ID_REF"

write.table(res, file = "${output.counts.tsv}", sep = "\t", quote=FALSE, row.names=FALSE);

