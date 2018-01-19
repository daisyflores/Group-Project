working_directory <-'/Users/Daisy/Desktop/MultivariateStats/multivar'
#COMPUTATION OF FD USING MULTIDIM FD####
source("quality_funct_space.R")
source("plot_funct_space.R")
source("multidimFD.R")

# 1.1. ORDINATION OF SPECIES BASED ON MORPHOLOGICAL TRAITS####
# Import traits matrix
traits<- read.table(paste0(working_directory, '//urchingtraits.txt'), 
                    h=T, row.names = 1)
#Delete Species that is not in our abundance data set.
traits <- traits[-8,] #echinotrix.diademis
# Load the dataset
biomass <-(read.delim(file = paste0(working_directory, '//data2_species.txt'), 
                      h=T, sep ="\t"))

# imported the abundance instead of the bites and transform it into matrix
require(reshape)
biomass.wide <- cast(biomass, siteID~species, sum) 
#save table
write.table(biomass.wide, file="Exam.biomass.wide.txt",sep = "\t") 

biomass.wide.df <- data.frame(biomass.wide[,-1], row.names = biomass.wide[,1])

abundance <-as.matrix(biomass.wide.df)
#checking if colnames of the species is matching the other df
sum(row.names(traits) %in% colnames(abundance)) == ncol(abundance)
identical(rownames(traits), colnames(abundance))
# computing all functional spaces based on PCoA (up to 6 axes) BECAUSE WE ONLY HAVE 6 TRAITS SEE nbdim=6
# NOTE: you need to be connected to internet for a correct functioning of quality_funct_space function
qual_funct_space<-quality_funct_space(traits, 
                                      traits_weights=NULL, 
                                      nbdim=6, metric="Euclidean",
                                      dendro=F, plot="quality_funct_space")
qual_funct_space$meanSD


# extracting species coordinates in the 4D space
coord_4D<-qual_funct_space$details_funct_space$mat_coord[,1:4]
# saving them
save(coord_4D, file="coord_4D")

# computing functional diversity indices 
FD_baskets<-multidimFD(coord_4D,abundance, check_species_pool=F, verb=TRUE)

# printing results = rounded FD indices values
IndicesAll<-round(FD_baskets,3)
Indices<-IndicesAll[,c(1,19,21:22)]
Indices

# exporting to excel
setwd("/Users/Daisy/Desktop/MultivariateStats/multivar") 
write.table(Indices, file="ExamFIndicesMultidimFD.txt", sep="\t") 


