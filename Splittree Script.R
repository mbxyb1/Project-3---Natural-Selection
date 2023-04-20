# IMPORT DATA (from SNP file)
library(adegenet)
library(vcfR)
library(StAMPP)

### vcf import

#first define a function for tetraploids
vcfR2genlight.tetra <- function (x) 
{
  bi <- is.biallelic(x)
  if (sum(!bi) > 0) {
    msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
    msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
    msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
    warning(msg)
    x <- x[bi, ]
  }
  CHROM <- x@fix[, "CHROM"]
  POS <- x@fix[, "POS"]
  ID <- x@fix[, "ID"]
  x <- extract.gt(x)
  x[x == "0|0"] <- 0
  x[x == "0|1"] <- 1
  x[x == "1|0"] <- 1
  x[x == "1|1"] <- 2
  x[x == "0/0"] <- 0
  x[x == "0/1"] <- 1
  x[x == "1/1"] <- 2
  x[x == "1/1/1/1"] <- 4
  x[x == "0/1/1/1"] <- 3
  x[x == "0/0/1/1"] <- 2
  x[x == "0/0/0/1"] <- 1
  x[x == "0/0/0/0"] <- 0
    x <- adegenet::as.genlight(t(x))
  adegenet::chromosome(x) <- CHROM
  adegenet::position(x) <- POS
  adegenet::locNames(x) <- ID
  x
}

vcf <- read.vcfR("reheadered_4dg_dips.clean_BI.ann.vcf")
aa.genlight <- vcfR2genlight.tetra(vcf)

## OPTIONAL add pop names: here pop names are first 5 chars of ind name
pop(aa.genlight)<-substr(indNames(aa.genlight),1,3)

# Calculate genetic distance between individuals/pops
aa.D.ind <- stamppNeisD(aa.genlight, pop = FALSE)  # Nei's 1972 distance between indivs
aa.D.pop <- stamppNeisD(aa.genlight, pop = TRUE)   # Nei's 1972 distance between pops

# Export the genetic distance matrix in Phylip format - for SplitsTree
stamppPhylip(aa.D.ind, file="aa.indiv_Neis_distance.phy.dst")
stamppPhylip(aa.D.pop, file="aa.pop_Neis_distance.phy.dst")

#Open in Splits tree !!!


###############################
#   OPTIONAL - SUBSETTING
##################################

### OPTION 1) extract specific indivs (but does not retain groupnames)
#get the names of indivs
aa.names<-row.names(as.matrix(aa.genlight))  # get just the vector of indiv names
write.table(aa.names, file="names_all.txt")
aa.matrix.BSW <- as.matrix(aa.genlight[c(10:16,25:27,29:31,33:49,58:85,93:119,125:132,135:204,210,212:279),]) # without weird samples
names.BSW <- names[c(10:16,25:27,29:31,33:49,58:85,93:119,125:132,135:204,210,212:279)]  # numbers = indices of the samples in vecotr of sample names
aa.genlight.BSW <- new("genlight",aa.matrix.BSW)
# OPTIONAL - pop names
names.BSW.pops<-substr(names.BSW,1,3)    # select first three characters from samplename = it will be a pop
aa.genlight.BSW@pop<-as.factor(names.BSW.pops)


### OPTION 2) separate by each population 
# generates a list of genlight objects, one per pop.
aa.genlight.sep <- seppop(aa.genlight)  # separate genlight per population





