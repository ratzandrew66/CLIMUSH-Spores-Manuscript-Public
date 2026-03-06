# Project: CLIMUSH: ITS OTU QC & Creating a Phyloseq Object
# Author: Peter Kennedy & Andrew Ratz
# Date: 10/21/24

# load libraries
library(here)
library(dplyr)
library(vegan)
library(phyloseq)
library(tidyverse)
library(stringr) ## for manipulating character strings
#install.packages("readxl")
library(readxl)

# set working directory
setwd("C:/Users/ratza/OneDrive/Desktop/Kennedy Lab Work/CLIMUSH")

# Load in ASV table (csv, txt, tab delimited file, etc.)
data <- read.delim("ITS_combined_sequences_clustered97_mumu_curated_with_taxonomy_PK.txt")

####Step 1 - tidy up the raw data table

# Looking at column names
colnames(data)

# Re-organize the columns to move "OTUId" column to last
data2 <- data %>%
  select(-c(1), everything(), c(1))

colnames(data2)

# To add an OTU_ID column, need to know number of rows
dim(data2)

# Add a column called "otu_id" and place it as the first column in the data frame
data3 <- data2 %>%
  mutate(otu_id = paste0("OTU", sprintf("%04d", 1:16114))) %>%  # the sprintf helps keeping the zeroes before the numbers
  select(otu_id, everything())

# looking at the column names to make adjustments
colnames(data3)

# Edit column names for mock communities
colnames(data3)[828] <- "synmock"
colnames(data3)[826] <- "ArnoldMock1"
colnames(data3)[827] <- "ArnoldMock2"

# Edit column names for negative controls
colnames(data3)[825] <- "neg1"

# Looking at the column names to check adjustments
colnames(data3)

# Eliminating extra columns

data3$Strand <- NULL

# Looking at the column names to check adjustments
colnames(data3)

# Extract the genus after g: and before the next comma and puts in new column "Genus_80"
data3 <- data3 %>%
  mutate(Genus_80 = str_extract(Tax_above_threshold, "(?<=g:)[^,]+"))
## these are all the taxa with >80% bootstrap support for the genus ID, which we feel comfortable giving a functional trait annotation

##### Step 2 - Assess mock community

##Look at the mock community sample(s) and see if it contains the number of OTUs expected.

# Counting non-zero cells for the synmock
sum(data3$synmock != 0)  #74 non-zero OTUs

# Sort synmock column in decreasing order and display top 50 values
sorted_synmock <- sort(data3$synmock, decreasing = TRUE)
top_50_synmock <- head(sorted_synmock, 50)
print(top_50_synmock)
temp_mocks <- data3 %>% select(otu_id, synmock, ArnoldMock1, ArnoldMock2)
##only 4 OTUs >1000 reads - that is 7 less than the # of OTUs expected. Next most abundant is 14 reads - good clustering, some index bleed.

# Eliminating the synmock column
data3$synmock <- NULL

# Counting non-zero cells for the synmock
sum(data3$ArnoldMock1 != 0)  #44 non-zero OTUs
sum(data3$ArnoldMock2 != 0)  #56 non-zero OTUs

# Sort ArnoldMock1 column in decreasing order and display top 50 values
sorted_ArnoldMock1 <- sort(data3$ArnoldMock1, decreasing = TRUE)
top_50_ArnoldMock1 <- head(sorted_ArnoldMock1, 50)
print(top_50_ArnoldMock1)
##only 2 OTUs >500 reads - that is 21 less than the # of OTUs expected. Next most abundant is 7 reads - good clustering, some index bleed.

# Sort ArnoldMock2 column in decreasing order and display top 50 values
sorted_ArnoldMock2 <- sort(data3$ArnoldMock2, decreasing = TRUE)
top_50_ArnoldMock2 <- head(sorted_ArnoldMock2, 50)
print(top_50_ArnoldMock2)
##only 2 OTUs >100 reads - that is 21 less than the # of OTUs expected. Next most abundant is 28 reads (next after that 8) - good clustering, some index bleed.

# Eliminating the other mock columns
data3$ArnoldMock1 <- NULL
data3$ArnoldMock2 <- NULL

# looking at the column names to check adjustments
colnames(data3)

#####Step 3 - Separate Sequence and Taxonomic Datasets in OTU table

# Before separating seq and tax data, sort to only keep Fungi

# List which kingdoms are present
table(data3$Kingdom)

# now subset out kingdom to only keep Fungi
data4.fungi <-data3 %>%
  dplyr::filter(Kingdom=="Fungi")

# confirm Fungi are the only kingdom left
table(data4.fungi$Kingdom)

# ID which columns to keep that have sequence data
colnames(data4.fungi)

# now subset out columns with the sequence data only, but make sure to include otu_id
data4.fungi.seq <- data4.fungi[,c(1:825)] %>%
  tibble::column_to_rownames("otu_id")

# confirm the number of OTUs and samples:
dim(data4.fungi.seq)

# ID which columns to keep that have taxonomy data
colnames(data4.fungi)
sum(data4.fungi[,c(2:824)])
# now subset out columns with the taxonomic data only, but make sure to include otu_id
data4.fungi.tax <- data4.fungi[,c(1,826:842)] %>%
  tibble::column_to_rownames("otu_id")

#check dimensions
dim(data4.fungi.tax) #should match sequence matrix for rows

####Step 3 - Assess the negative controls

# ID columns with negative controls
colnames(data4.fungi.seq)

# Examine read counts in neg controls
sum(data4.fungi.seq$neg1)
##contamination in neg controls: 1928 reads total

# Sort neg control column in decreasing order and display top 50 values
sorted_neg <- sort(data4.fungi.seq$neg1, decreasing = TRUE)
top_50_neg <- head(sorted_neg, 50)
print(top_50_neg)

#Adding controls together
data4.fungi.seq.neg <- data4.fungi.seq$neg1 + data4.fungi.seq$neg2

#checking sum higher than individual neg control columns
sum(data4.fungi.seq.neg)

#Subtracting negative control reads from each sample
data4.fungi.seq.nc <-
  sweep(data4.fungi.seq, 1, data4.fungi.seq$neg1, FUN="-")

#Converting any negative values to zero
data4.fungi.seq.nc[data4.fungi.seq.nc  < 0] <- 0

# Eliminating any samples that now have 0 total reads (this eliminates the negative control column
# plus any others with zero reads total)
data4.fungi.no.col.zero <- data4.fungi.seq.nc %>%
  select(-where(~ sum(.) == 0))

# Checking new dimensions to determine how many columns were eliminated
dim(data4.fungi.no.col.zero)
## Only the negative control column was eliminated

# Eliminating any OTUs that now have 0 reads
data4.fungi.no.zero <- data4.fungi.no.col.zero %>%
  filter(rowSums(across(everything())) != 0)

#Check new dimensions to determine how many rows were eliminated
dim(data4.fungi.no.zero)
## Went down from 13126 to 13124 OTUs

#Comparing the first 5 rows and columns of the raw and trimmed data frame to see if they are different
data4.fungi.seq[1:5,1:5]
data4.fungi.no.zero[1:5,1:5]
#---- Proportional Index bleed corrections vs Absolute ----

####Step 4 - Account for index bleed (absolute)

## Based on the synmock counts, it appears that any cell values less than 15 are suspect. The start to the low abundance OTUs is at 14, so I am going to set the minimum cell value at 14.

# Converting values in any cell with a value of <=14 to zero
data4.fungi.no.zero.thres <- data4.fungi.no.zero %>%
  mutate(across(everything(), ~ ifelse(. <= 14, 0, .)))

# Check if all values in each cell have counts above 14
counts_above_14 <- data4.fungi.no.zero.thres %>%
  summarize(across(everything(), ~ sum(. > 14, na.rm = TRUE) == n()))

# Print result (the outcome should say FALSE for all samples)
print(counts_above_14)

# Eliminating any OTUs that now have 0 reads
data4.fungi.no.zero.thres <- data4.fungi.no.zero.thres %>%
  filter(rowSums(across(everything())) != 0)

#Check new dimensions
dim(data4.fungi.no.zero.thres)
## Went down from 13124 to 7572 OTUs --> lost 5552 OTUs

#lets quickly check what kinds of OTUs were lost
lost.OTUs <- data4.fungi %>%
  filter(!otu_id %in% rownames(data4.fungi.no.zero.thres))

lost.OTUs.summary <- data.frame(Sum = rowSums(lost.OTUs[,c(2:825)]),
                   Max = apply(lost.OTUs[,c(2:825)], 1, max))
rownames(lost.OTUs.summary) <- lost.OTUs$otu_id

lost.OTUs.tax <- lost.OTUs[,c(1,826:842)]
view(lost.OTUs.tax)

# Renaming trimmed and minimum thresholded dataframe

data5.fungi.seq <- data4.fungi.no.zero.thres
sum(data5.fungi.seq)
## 12,318,854 total reads across 7572 OTUs


####Step 4(2) Proportional index bleed

#data4.fungi.no.zero.prop <- data4.fungi.no.zero
#Choosing a cutoff of 0.5% due to that being the max ratio of index bleed  between the three mock communities
#data4.fungi.no.zero.prop[data4.fungi.no.zero < (row_totals * 0.005)] <- 0


# Eliminating any OTUs that now have 0 reads
#data4.fungi.no.zero.prop <- data4.fungi.no.zero.prop %>%
  #filter(rowSums(across(everything())) != 0)

#data5.fungi.seq <- data4.fungi.no.zero.prop


#check total number of fungal reads


#####Step 5 - Assess sample sequence count totals

## Some samples can have very few total sequences and as such are not likely good to include in the final analyses. However, determining where to cut-off samples in terms of read totals is notably subjective.

# Create a list of column sums
column_sums <- data5.fungi.seq %>%
  summarize(across(everything(), sum))

# Visualize the distribution of counts across all samples
column_sums <-as.numeric(unlist(column_sums))

hist(log10(column_sums), breaks = 100)
abline(v = log2(1020), lty = 2, col = 'red')
abline(v = log2(2060), lty = 2, col = 'red')

#----rarefaction curve?----
flipped <- as.data.frame(t(data5.fungi.seq))
flipped_cleaned <- flipped[rowSums(flipped) > 1, ]

# Calculate the curves
# step = 100 (draws a point every 100 reads to save memory)
# sample = 1020 (draws a vertical line at your current cutoff)

rarecurve(flipped_cleaned, step = 200, sample = 2063, label = FALSE, 
          col = "blue", cex = 0.6,
          main = "Rarefaction Curves: Species Richness vs. Depth",
          xlab = "Sequencing Depth (Reads)", 
          ylab = "Number of OTUs")

# Add a vertical line at your current cutoff to see where it hits the curves
abline(v = 1020, col = "red", lty = 2, lwd = 2)
abline(v = 2063, col = "red", lty = 2, lwd = 2)
abline(v = 10000, col = "red", lty = 2, lwd = 2)


#another test of samples
sample_stats <- data5.fungi.seq %>%
  pivot_longer(cols = everything(), 
               names_to = "Sample", 
               values_to = "count") %>%
  group_by(Sample) %>%
  summarise(
    Reads = sum(count),
    nOTUs = sum(count > 0)
  ) %>%
  filter(Reads > 0)

ggplot(sample_stats_15, aes(x = log(Reads), y = nOTUs))+
  geom_point()+
  geom_smooth(method = 'lm')

summary(lm(nOTUs ~ log(Reads), data = sample_stats_15))

#----Continue with QC----
# Re-create a list of column sums
column_sums <- data5.fungi.seq %>%
  summarize(across(everything(), sum))

# Pivot the column sums object to a long list
column_sums_long <- column_sums %>%
  pivot_longer(cols = everything(), names_to = "Column_Name", values_to = "Sums")

# Print the ranked list of 50 lowest column sums
column_sums_ranked <- column_sums_long %>%
  mutate(rank = rank(Sums)) %>%
  arrange(Sums)
print(column_sums_ranked)

view(column_sums_ranked)

## There are 17 of the 823 samples that have 0 reads.
## Based on the read counts, if we cut-off all samples with < 2000 reads, we lose 70 samples.

#Will rarefy to 1020 --> 2063

# Set the threshold for column sums
threshold <- 2063

# Drop columns with sums less than the threshold
data6.fungi.seq <- data5.fungi.seq %>%
 select(where(~ sum(.) >= threshold))




dim(data6.fungi.seq)
##753 samples have >2000 reads/sample

####Step 5 - tidying up the taxonomy dataframe to match the final sequencing dataframe

# Convert taxonomy df row names to a column
df_tax <- rownames_to_column(data4.fungi.tax, var = "RowName")

# Convert sequence df row names to a column
df_seq <- rownames_to_column(data6.fungi.seq, var = "RowName")

# Perform the join based on df_seq
df_joined <- df_seq %>%
  inner_join(df_tax, by = "RowName")

# list column names
colnames(df_joined)

# change first column name to "otu_id"
colnames(df_joined)[1] <- "otu_id"

# list column names
colnames(df_joined)

# getting tax and seq parts of df split back into separate parts

data5.fungi.tax <- df_joined[,c(1,755:771)] %>%
  tibble::column_to_rownames("otu_id")

data6.fungi.seq <- df_joined[,c(1:754)] %>%
  tibble::column_to_rownames("otu_id")

#####Step 6 - add FungalTraits data to the taxonomic dataset (based off code from Talia)

dim(data5.fungi.tax)
colnames(data5.fungi.tax)

## Load in FungalTraits database:
traits <- read.csv("C:/Users/ratza/OneDrive/Desktop/Kennedy Lab Work/CLIMUSH/FungalTraitsDatabase.csv")

colnames(traits)

temp.traits <- traits %>% select(
  Phylum, Class, Order, Family, GENUS, primary_lifestyle, Fruitbody_type_template
)

##rename GENUS to Genus_80
traits <-traits %>%
  dplyr::rename("Genus_80"= "GENUS")
names(traits)

###
names(data5.fungi.tax)
table(data5.fungi.tax$Genus_80)

##Note There are a few cases where GENUS is NA:

traitsOfInterest <- select(traits, Genus_80, primary_lifestyle, Fruitbody_type_template)
traitsOfInterest
colnames(data5.fungi.tax)
class(data5.fungi.tax$Genus_80)
data5.fungi.tax$Genus_80<-as.factor(data5.fungi.tax$Genus_80)
class(data5.fungi.tax$Genus_80)
traitsOfInterest$Genus_80<-as.factor(traitsOfInterest$Genus_80)

# joining tax with traits database
data6.fungi.tax <- left_join(data5.fungi.tax, traitsOfInterest, by="Genus_80")

## unfortunately this loses the otu_ids as the row names, so we need to add those back.
data6.temp <- data6.fungi.tax %>% filter(!is.na(Genus_80), !str_detect(Genus_80, "_sedis"))
# Save the row names from data5.fungi.tax
row_names <- rownames(data5.fungi.tax)

# Apply the row names to data5.fungi.tax.trait
rownames(data6.fungi.tax) <- row_names

# viewing the join
View(data6.fungi.tax)

dim(data6.fungi.tax)
dim(data5.fungi.tax)

table(data6.fungi.tax$primary_lifestyle)
table(data6.fungi.tax$Kingdom)

####Step 7 - Loading the mapping file (i.e. metadata file)

## Will need to first modify the column names to retain only the sample information in the sequence file so that it can be matched with the SampleName or SampleID in the mapping file. The function below keeps everything in front of the period in the column name.


# take the column names and put in a vector
col.names.fungi <- colnames(data6.fungi.seq)

# write that vector to a table
write.table(col.names.fungi, "col.names.fungi.csv", row.names = F)

# import the colnames table as a data object
col.names.fungi <-read.csv("C:/Users/ratza/OneDrive/Desktop/Kennedy Lab Work/CLIMUSH/col.names.fungi.csv")
dim(col.names.fungi)

# add "Sample_ID" to column header
colnames(col.names.fungi)[1] <- "Sample_ID"

# read in mapping file
map<-read.csv("C:/Users/ratza/OneDrive/Desktop/Kennedy Lab Work/CLIMUSH/CLIMUSH_mapping.csv")
dim(map)

# confirm that mapping file and col.names have overlapping Sample_ID using left_join
data6.fungi.map <- left_join(col.names.fungi, map, by = "Sample_ID") %>%
  tibble::column_to_rownames("Sample_ID")

# look at new map to confirm everything works
View(data6.fungi.map)
View(data6.fungi.seq)

####Step 8 - Putting it all into a PhyloSeq Object

# Transform otu and tax tables intro matrices
OTU.seq <- as.matrix(data6.fungi.seq)
OTU.tax <- as.matrix(data6.fungi.tax)

# check dimensions
dim(OTU.seq)
dim(OTU.tax)

## Mapping file can be left as data frame

# check dimensions
dim(data6.fungi.map)

OTU = otu_table(OTU.seq, taxa_are_rows = TRUE) # set your OTU/ASV file
TAX = tax_table(OTU.tax) # set you Tax data file
MAP = sample_data(data6.fungi.map) # set your Map file

# Converting matrices into a phyloseq object!
CLIMUSH.phylo.otu <- phyloseq(OTU, TAX, MAP)

# Save this phyloseq object
saveRDS(CLIMUSH.phylo.otu, "CLIMUSH.phylo.otu.rds")

####Step 9 - Create a rarefied version of your data

#----Perform rarefaction to an even depth----
set.seed(123)  # Setting seed for reproducibility
CLIMUSH.phylo.rare.otu <- rarefy_even_depth(CLIMUSH.phylo.otu, rngseed = FALSE, sample.size = min(sample_sums(CLIMUSH.phylo.otu)))

# Get the sample sums
sample_sums_rarefied <- sample_sums(CLIMUSH.phylo.rare.otu)

# Print the sample sums
print(sample_sums_rarefied) #All samples rarefied to 1020 reads

##Save this final phyloseq object
saveRDS(CLIMUSH.phylo.rare.otu, "CLIMUSH.phylo.rare.otu.rds")

#----Finally, test Goods Coverage index for rarefied data----
otu_matrix <- as(otu_table(CLIMUSH.phylo.rare.otu), "matrix")

# 2. Convert to a dataframe and move OTU names to a column
otu_df <- as.data.frame(otu_matrix)
otu_df$OTU <- rownames(otu_df)

# 3. Pivot to long format
sample_stats_rare <- otu_df %>%
  pivot_longer(
    cols = -OTU, 
    names_to = "Sample", 
    values_to = "Abundance"
  )

sample_stats_rare <- sample_stats_rare %>% filter(Abundance > 0)

rare_Goods <- sample_stats_rare %>%
  group_by(Sample)%>%
  summarize(n1 = sum(Abundance == 1),
            N = sum(Abundance),
            Goods = 1-(n1/n))

rare_Goods <- rare_Goods %>%
  separate_wider_delim(Sample,".", names = c('sequencing_type','compartment','year',
                                                'month','day','domain','treatment'
                                                ,'subplot','tube_ID','OTU_placeholder'))%>%
  filter(domain %in% c("D19","D14","D13","D03","D06","D01","D05",'D16'))

ggplot(rare_Goods, aes(x = as.factor(domain), y = Goods))+
  geom_boxplot()+
  geom_hline(yintercept = 0.90, linetype = 'dashed', color = 'black')+
  geom_hline(yintercept = 0.95, linetype = 'dashed', color = 'red')+
  theme_classic()

