library(tidyverse)
library(ggplot2)

# ---- Initial Simple Mapping File ----

#reading in cleaned OTU table
DNA.dat <- read.delim("ITS_combined_sequences_clustered97_mumu_curated_with_taxonomy_PK.txt")

#creating a mapping file based on names found in the OTU data table
colnames(DNA.dat)
col.dat <- data.frame(colnames(DNA.dat[,-c(1,825:844)]))
colnames(col.dat) <- c("Sample_ID")

#breaking the name apart into each of the components 
test <- col.dat %>%
  separate_wider_delim(Sample_ID,".", names = c('sequencing_type','compartment','year',
                                                'month','day','domain','treatment'
                                                ,'subplot','tube_ID','OTU_placeholder'))
test$Sample_ID <- col.dat$Sample_ID
test <- test[,c(11,1,2,3,4,5,6,7,8,9,10)]

#fixing some odd Arizona datapoints - this will come up more....
test[760:776,7] <- "D14"

#creating burn history based on abbreiations
test <- test %>%
  mutate(burn_history = case_when ( treatment == 'BG' | treatment == 'BC' | treatment == 'BO' ~ "burned",
                                    treatment == 'UG' | treatment == 'UC' | treatment == 'UO' ~ "unburned"))
#creating habitat column
test <- test %>%
  mutate(habitat = case_when ( treatment == 'BG' | treatment == 'UG' ~ 'grassland',
                               treatment == 'BO' | treatment == 'UO' ~ 'oaks',
                               treatment == 'BC' | treatment == 'UC' ~ 'conifers'))
#making ecoregion and state more easy to determine based on the domain designation
test <- test %>%
  mutate(ecoregion = case_when ( domain == "D19" ~ 'TAIG',
                                 domain == "D14" ~ 'DSW',
                                 domain == "D13" ~ 'SROC',
                                 domain == "D03" ~ 'SE',
                                 domain == "D06" ~ 'PRAI',
                                 domain == "D01" ~ 'NE',
                                 domain == "D05" ~ 'GRTL',
                                 domain == "D16" ~ 'PNW'))

test <- test %>%
  mutate(state = case_when ( domain == "D19" ~ 'AK',
                             domain == "D14" ~ 'AZ',
                             domain == "D13" ~ 'CO',
                             domain == "D03" ~ 'FL',
                             domain == "D06" ~ 'KS',
                             domain == "D01" ~ 'MA',
                             domain == "D05" ~ 'MN',
                             domain == "D16" ~ 'OR'))


test$subplot <- gsub("0","",test$subplot)

#adding site metadata info from official CLIMUSH documentation
site_metadata <- read.csv("CLIMUSH_site_metadata.csv")

site_metadata <- site_metadata[c(1:279),c(2:4,6:9,14:16,33:35,40:42)] #selectng columns I want

colnames(site_metadata) <- c("site","state","ecoregion","subplot","habitat","fire","yrsincefire",
                             "Lat_dec","Long_dec","Elevation.m","soil_N",'soil_C','soil_pH','litt_N','litt_C','litt_pH')

site_metadata$subplot <- gsub("S","",site_metadata$subplot)

#creating the columns to join by
site_metadata <- site_metadata %>%
  mutate(burn_history = case_when(fire == "UB" ~ "U",fire == "B" ~ "B"))%>%
  mutate(treatment = paste(burn_history,habitat, sep = ""))

#removing unused columns before joining
site_metadata$habitat <- NULL
site_metadata$fire <- NULL
site_metadata$burn_history <- NULL

mapping <- left_join(test,site_metadata, by = c('state','subplot','treatment','ecoregion'))

write.csv(mapping, file = "CLIMUSH_mapping.csv")

#---- Alpha Diversity ----

library(vegan)

#Load in the rarefied phyloseq object:
df.rare<-readRDS("CLIMUSH.phylo.rare.otu.rds")
df.rare

###Step 1. 

#Filtering only Macrofungi
macro_included <- read.csv("Macrofungi_Defined.csv")
macro_excluded <- read.csv("Macrofungal_exclusions.csv")

df.rare.macrofungi <- df.rare%>%
  subset_taxa(Class %in% macro_included$Only_included_class & 
              (Order %in% macro_included$Only_included | Family %in% macro_included$Only_included) &
             !Family %in% macro_excluded$Excluded_Fam_Gen & 
             !Genus %in% macro_excluded$Excluded_Fam_Gen)


get_taxa_unique(df.rare.macrofungi, "Order")

# 1. Remove OTUs that now have 0 reads across all samples
df.rare.macrofungi <- prune_taxa(taxa_sums(df.rare.macrofungi) > 0, df.rare.macrofungi)

# 2. Remove Samples that now have 0 reads (because all their fungi were excluded)
df.rare.macrofungi <- prune_samples(sample_sums(df.rare.macrofungi) > 0, df.rare.macrofungi)



#'  ## First calculate alpha diversity metrics for ALL samples
#'  
#' ############################################################ 
#'  #' ## ALPHA DIVERSITY CALCULATION 
#'  
#' separate out the OTU table and transpose it:
OTU.df.rare<-df.rare.macrofungi@otu_table
OTU.df.rare.t<-t(OTU.df.rare)
OTU.df.rare.t<-as.matrix(OTU.df.rare.t)

rownames(OTU.df.rare.t)
colnames(OTU.df.rare.t)
dim(OTU.df.rare.t)
#3216 macrofungal OTUs across 744 samples

#' #' #Calculating observed OTU richness per sample:
Sobs <-specnumber(OTU.df.rare.t)
##Shannon Diversity
Shannon<-diversity(OTU.df.rare.t)
#calculate Pielous evenness:
evenness<- Shannon/log(Sobs)
#Simpson diversity
Simpson <- diversity(OTU.df.rare.t, index = "simpson")
##cbind these three alpha diversity metrics:
alpha.diversity<-cbind(Sobs,Shannon,
                       Simpson,evenness)
#" Set as data frame:
alpha.diversity<-data.frame(alpha.diversity)
rownames(alpha.diversity)

# Convert alpha div df row names to a column
alpha.div2 <- rownames_to_column(alpha.diversity, var = "Sample_ID")

#Bring in mapping file created earlier
mapping <- read.csv("CLIMUSH_mapping.csv")

# Left join the mapping file and alpha div df
alpha.div.table <- left_join(alpha.div2, mapping, by="Sample_ID") 

alpha.div.table$X <- NULL

colnames(alpha.div.table)

#Set working directory
#setwd("~/Library/CloudStorage/GoogleDrive-kennedyp@umn.edu/My Drive/Research_Active/CLIMUSH/Analyses")
setwd("C:/Users/ratza/OneDrive/Desktop/Kennedy Lab Work/CLIMUSH")
#write div table to .csv
write.csv(alpha.div.table, file = "alpha.div.table.csv")


#---- Adding small metadata ----

##---- qPCR and Diversity ----
#dealing with the qPCR and Div data

#reading in qPCR and Diversity tables made previously
QPCR.df <- read.csv("QC_report_master.csv")
alpha.div <- read.csv("alpha.div.table.csv")

#checking column names
colnames(alpha.div)
colnames(QPCR.df)

#changing the Sample_ID to have hyphens (correct naming convention) 
alpha.div$X <- NULL
QPCR.df$Sample_ID <- gsub("[.-]","_",QPCR.df$Sample_ID)
alpha.div$Sample_ID <- gsub("[.-]","_", alpha.div$Sample_ID)

#Remove the OTU placeholder after the last underscore
alpha.div$Sample_ID <- sub("_[^_]*$", "", alpha.div$Sample_ID)

#joining data frames
dat <- left_join(alpha.div, QPCR.df, by = 'Sample_ID')  

# Remove commas and convert to numeric for easier analysis
dat$copy_number <- as.numeric(gsub(",", "", dat$copy_number))

#testing
ggplot(data = dat,aes(x = log(copy_number)))+geom_histogram() #looks decent

#removing an odd dupe column
dat$copy_number_.molecules.ul. <- NULL

#creating an n-day column to look at seasonal effects
# Combine month and day into a date
dates <- as.Date(paste(dat$year, dat$month, dat$day, sep = "-"), format = "%Y-%m-%d")

# Calculate the number of days since the beginning of the year
start_of_year <- as.Date(paste(dat$year, "01", "01", sep = "-"), format = "%Y-%m-%d")
days_since_start <- as.numeric(dates - start_of_year)
dat$nday <- days_since_start

##---- Climatic Data ----

#climatic data from CLIMUSH_climate_data.R

clim.data <- read.csv("climate_dat.csv")
clim.data$X <- NULL

subplot_T <- read.csv("subplot_tempdata.csv")
subplot_T$X <- NULL

dat <- dat%>%
  left_join(clim.data, by = c("site","year","month","day"), 
            relationship = "many-to-many", multiple = "any")

##----Correcting Diversity by Number of days the traps were left out ----
collect_dat <- read.csv("CollectionDays.csv")
collect_dat <- collect_dat[,1:4]

dat <- left_join(dat, collect_dat, by = "Sample_ID")

dat$Nr.days.out <- as.numeric(dat$Nr.days.out)
dat$log.days.out <- log(dat$Nr.days.out)

## ---- Phenology ----

#because each state has a different fruiting season, I will have to do the season for each state differently


#after looking at each site's sampling dates, plus the fruiting graphs, I can split all of the sites into two groups that have different split days

dat <- dat%>%
  mutate(season = case_when(
    site %in% c("HARV","NWT") & nday < 195 ~ "early",
    site %in% c("HARV","NWT") & nday > 195 ~ "late",
    site %in% c("BNZ","SRE","KON","CDR","HJA","PIS") & nday < 155 ~ "early",
    site %in% c("BNZ","SRE","KON","CDR","HJA","PIS") & nday > 155 ~ "late",
    (site == "ORD" & nday < 155) | (site == "ORD" & nday > 310) ~ "early",
    site == "ORD" & nday > 155 & nday < 310 ~ "late",
  ))


##---- Long dataframe ----

#adding the metadata for the samples and QPCR counts
df.rare.long <- df.rare.macrofungi %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  psmelt() %>%
  filter(Abundance > 0)

df.rare.long <- df.rare.long %>%
  rename("Sample_ID" = "Sample")%>%
  mutate(Sample_ID = gsub("[.]","_",Sample_ID))%>%
  mutate(Sample_ID = sub("_[^_]*$","", Sample_ID))
         
macro.data.long <- df.rare.long %>%
  left_join(dat[,c(1:5,34:50)], by = "Sample_ID", relationship = "many-to-many", multiple = "any")%>%
  select(-X)

#Fixing an incorrect Florida sample
dat <- dat %>% filter(!Sample_ID == "illumina_spore_23_06_07_D03_UG_08_79")
macro.data.long <- macro.data.long %>% filter(!Sample_ID == "illumina_spore_23_06_07_D03_UG_08_79")


#----Writing the final metadata frames----
write.csv(dat, "dat.csv")
write.csv(macro.data.long, "macro_long.csv")

