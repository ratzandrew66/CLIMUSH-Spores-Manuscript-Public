# CLIMUSH-Spores-Manuscript-Public
Processed data files and code used to produce all analyses and figures for the manuscript "Sporecasting biogeography across scales: continental patterns and local constraints on macrofungal dispersal"

CLIMUSH_ITS_OTU_QC was used to process the OTU table output from DADA2 and taxonomic curation detailed within the 'Molecular Quantification and Identification' section of the methods

CLIMUSH_metadata was used to generate a metadata file from the sample IDs encoded in the OTU table and join this with measured climatic variables

dat.csv is the output of this file, encoding all neccesary metadata by Sample and is used in conjunction with the long dataframe
macro_long is the long-format of the OTU table along with the metadata and is used for many of the analyses within the Results.Rmd


CLIMUSH_Results contains all analysss not including the iNaturalist dataframe

CLIMUSH_iNat contains all analyses that do include the community data pulled from iNaturalist, along with processing the raw pulled data into a workable dataframe


