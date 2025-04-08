##########################################
## MICROBIOME TRANSFER IN A POISON FROG ##
##           ANGIOLANI ET AL            ##
##                SCRIPT                ##
##########################################

#GLOSARY----
# In cases when Posit_ctrl/cpos/ is used, we refer to the "Normal transport" treatment
# In cases when Neg_ctrl/cneg/ is used, we refer to the "No transport" treatment
# In cases when Test is used, we refer to the "Foster transport" treatment

#LOAD PACKAGES ----
library(vegan)
library(ggplot2)
library(veganEx)
library(phyloseq)
library(microbiome) 
library(tidyverse)
library(microViz)
library(decontam)
library(MiscMetabar)
library(ComplexUpset)
library(dplyr)
library(DHARMa)
library(glmmTMB)
library(emmeans)
library(car)
library(metagMisc)

# 1.- Test differences between adult samples (ventral vs. dorsal) ----

#load specific data
otufile     <- read.table("initial_otu_table.txt", h=T, row.names = 1)
mapfile     <- read.table("initial_metadata.txt", h=T, row.names = 1)
taxfile     <- read.table("initial_tax_table.txt", h=T, row.names = 1) 

otufile <- as.matrix(otufile)
taxfile <- as.matrix(taxfile)

OTU = otu_table(otufile, taxa_are_rows = TRUE)
TAX = tax_table(taxfile)
samples = sample_data(mapfile)

#create pso
database.pso <- phyloseq(OTU, TAX, samples)

# Cleaning ps object
database.pso <- tax_fix(database.pso, unknowns = 'NA' )

# PERMANOVA FOR DIFFERENCES BY SIDE ON ADULT SWABS

Males <- subset_samples(database.pso, Adult == "Male")
Females <- subset_samples(database.pso, Adult =="Female")

#MALES
adonis_pq(Males, "Side", na_remove = TRUE)

#FEMALES
adonis_pq(Females, "Side", na_remove = TRUE)


# 2.-Remove contamination from database ----

#Identify Contaminants with witness with a treshold of 0.1

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant)) #which OTUs are contaminants ?

#with a treshold of 0.5
contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
head(which(contamdf.prev05$contaminant)) 

write.table(contamdf.prev05, "contaminant-0.5.txt", sep="\t")
write.table(contamdf.prev, "contaminant-0.1.txt", sep="\t")


# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

## ELIMINATE CONTAMINANT zOTUs FROM DATABASES WITH 0.1 THRESHOLD AND CREATE NEW PSO

# - Load general data ----

otufile     <- read.table("ZOTU_c97_Count_GOOD.txt")
mapfile     <- read.table("METADATA_GOOD.txt", sep = "\t", h=T, row.names = 1)
taxfile     <- read.table("OTU_tax_GOOD.txt", h=T, sep = "\t", row.names = 1) 

otufile <- as.matrix(otufile)
taxfile <- as.matrix(taxfile)

OTU = otu_table(otufile, taxa_are_rows = TRUE)
TAX = tax_table(taxfile)
samples = sample_data(mapfile)

database.pso <- phyloseq(OTU, TAX, samples)
# cleaning object
database.pso <- tax_fix(database.pso, unknowns = 'NA' )


# 3.-Sample type similarity analyses ----
OTU = as.data.frame(t(otufile))

#compute bray-curtis distances
bray <- vegdist(OTU,method="bray",binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE) 

#create bray-curtis dataframe
braydf = (as.data.frame(as.matrix(bray))) 
write.table(braydf, "bray-curtis.txt", sep="\t")
mapfile <- read.table("METADATA_GOOD.txt")

#mMDS
mMDS<- metaMDS(bray,k=2,stress=0.1)

# anosim to test the significance
anosim(bray, mapfile$GruopID, permutations = 999, distance = "bray")

#pairwise post hoc
anosim.pairwise(OTU, mapfile$GruopID, sim.function = "vegdist",
  sim.method = "bray",p.adjust.m = "bonferroni", perm = 999)

# 4.- Alpha-diversity analyses ----
richness <- estimate_richness(database.pso)

#### CHAO1

# Check normality
hist(richness$Chao1, main = "Chao1 index", xlab = "")
shapiro.test(richness$Chao1)

# test
kruskal.test(richness$Chao1 ~ sample_data(database.pso)$GruopID)

# Perform pairwise Wilcoxon tests and extract W statistics
unique_groups <- unique(sample_data(database.pso)$GruopID)
W_matrixCHAO <- matrix(NA, nrow = length(unique_groups), ncol = length(unique_groups))
rownames(W_matrixCHAO) <- unique_groups
colnames(W_matrixCHAO) <- unique_groups

for (i in 1:(length(unique_groups) - 1)) {
  for (j in (i + 1):length(unique_groups)) {
    group1 <- unique_groups[i]
    group2 <- unique_groups[j]
    
    # Subset the Chao1 diversity values for the two groups
    data1 <- richness$Chao1[sample_data(database.pso)$GruopID == group1]
    data2 <- richness$Chao1[sample_data(database.pso)$GruopID == group2]
    
    # Perform the Wilcoxon test
    test_result <- wilcox.test(data1, data2)
    
    # Store the W statistic in the matrix
    W_matrixCHAO[group1, group2] <- test_result$statistic
    W_matrixCHAO[group2, group1] <- test_result$statistic
  }
}

# Print the matrix of W statistics
print(W_matrixCHAO)

# Perform pairwise Wilcoxon test with p-value adjustment
pairwise_result <- pairwise.wilcox.test(richness$Chao1, 
                                        sample_data(database.pso)$GruopID, 
                                        p.adj = "bonf")
print(pairwise_result)

#### SHANNON

hist(richness$Shannon, main = "Shannon index", xlab = "")
shapiro.test(richness$Shannon)

kruskal.test(richness$Shannon ~ sample_data(database.pso)$GruopID)

unique_groups <- unique(sample_data(database.pso)$GruopID)
W_matrixSH <- matrix(NA, nrow = length(unique_groups), ncol = length(unique_groups))
rownames(W_matrixSH) <- unique_groups
colnames(W_matrixSH) <- unique_groups

for (i in 1:(length(unique_groups) - 1)) {
  for (j in (i + 1):length(unique_groups)) {
    group1 <- unique_groups[i]
    group2 <- unique_groups[j]

    data1 <- richness$Shannon[sample_data(database.pso)$GruopID == group1]
    data2 <- richness$Shannon[sample_data(database.pso)$GruopID == group2]
    
    test_result <- wilcox.test(data1, data2)

    W_matrixSH[group1, group2] <- test_result$statistic
    W_matrixSH[group2, group1] <- test_result$statistic
  }
}

print(W_matrixSH)

pairwise_result <- pairwise.wilcox.test(richness$Shannon, 
                                        sample_data(database.pso)$GruopID, 
                                        p.adj = "bonf")
print(pairwise_result)


##### INVERSE SIMPSON

hist(richness$InvSimpson, main = "InvSimpson index", xlab = "")
shapiro.test(richness$InvSimpson)

kruskal.test(richness$InvSimpson ~ sample_data(database.pso)$GruopID)

unique_groups <- unique(sample_data(database.pso)$GruopID)
W_matrixIS <- matrix(NA, nrow = length(unique_groups), ncol = length(unique_groups))
rownames(W_matrixIS) <- unique_groups
colnames(W_matrixIS) <- unique_groups

for (i in 1:(length(unique_groups) - 1)) {
  for (j in (i + 1):length(unique_groups)) {
    group1 <- unique_groups[i]
    group2 <- unique_groups[j]

    data1 <- richness$InvSimpson[sample_data(database.pso)$GruopID == group1]
    data2 <- richness$InvSimpson[sample_data(database.pso)$GruopID == group2]

    test_result <- wilcox.test(data1, data2)

    W_matrixIS[group1, group2] <- test_result$statistic
    W_matrixIS[group2, group1] <- test_result$statistic
  }
}

print(W_matrixIS)

pairwise_result <- pairwise.wilcox.test(richness$InvSimpson, 
                                        sample_data(database.pso)$GruopID, 
                                        p.adj = "bonf")
print(pairwise_result)

#### ESTIMATING ALPHA DIVERSITY PER OFFSPRING SAMPLE TYPE

offspring <- subset_samples(database.pso, Var == "offspring")

richness_offspring <- estimate_richness(offspring)
head(richness_offspring)

#Clutch week 1
CW1 <- subset_samples(offspring, GruopID == "Clutch_W1")
richness_CW1 <- estimate_richness(CW1)

hist(richness_CW1$Shannon, main="Shannon index", xlab="")
shapiro.test(richness_CW1$Shannon)
summary(aov(richness_CW1$Shannon ~ sample_data(CW1)$mockT))

hist(richness_CW1$Chao1, main="Chao1 index", xlab="")
shapiro.test(richness_CW1$Chao1)
summary(aov(richness_CW1$Chao1 ~ sample_data(CW1)$mockT))

hist(richness_CW1$InvSimpson, main="Inverse Simpson index", xlab="")
shapiro.test(richness_CW1$InvSimpson)
kruskal.test(richness_CW1$InvSimpson ~ sample_data(CW1)$mockT)

#Clutch week 2
CW2 <- subset_samples(offspring, GruopID == "Clutch_W2")
richness_CW2 <- estimate_richness(CW2)

hist(richness_CW2$Shannon, main="Shannon index", xlab="")
shapiro.test(richness_CW2$Shannon)
summary(aov(richness_CW2$Shannon ~ sample_data(CW2)$mockT))

hist(richness_CW2$Chao1, main="Chao1 index", xlab="")
shapiro.test(richness_CW2$Chao1)
summary(aov(richness_CW2$Chao1 ~ sample_data(CW2)$mockT))

hist(richness_CW2$InvSimpson, main="Inverse Simpson index", xlab="")
shapiro.test(richness_CW2$InvSimpson)
summary(aov(richness_CW2$InvSimpson ~ sample_data(CW2)$mockT))


#Tadpoles
TP <- subset_samples(offspring, GruopID == "Tadpole")
richness_TP <- estimate_richness(TP)

hist(richness_TP$Shannon, main="Shannon index", xlab="")
shapiro.test(richness_TP$Shannon)
kruskal.test(richness_TP$Shannon ~ sample_data(TP)$mockT)
pairwise.wilcox.test(richness_TP$Shannon, sample_data(TP)$mockT, p.adj = "bonf")

unique_groups <- unique(sample_data(TP)$mockT)
WTPmatrix <- matrix(NA, nrow = length(unique_groups), ncol = length(unique_groups))
rownames(WTPmatrix) <- unique_groups
colnames(WTPmatrix) <- unique_groups

for (i in 1:(length(unique_groups) - 1)) {
  for (j in (i + 1):length(unique_groups)) {
    group1 <- unique_groups[i]
    group2 <- unique_groups[j]

    data1 <- richness$Shannon[sample_data(TP)$mockT == group1]
    data2 <- richness$Shannon[sample_data(TP)$mockT == group2]

    test_result <- wilcox.test(data1, data2)

    WTPmatrix[group1, group2] <- test_result$statistic
    WTPmatrix[group2, group1] <- test_result$statistic
  }
}

print(WTPmatrix)


hist(richness_TP$Chao1, main="Chao1 index", xlab="")
shapiro.test(richness_TP$Chao1)
kruskal.test(richness_TP$Chao1 ~ sample_data(TP)$mockT)

hist(richness_TP$InvSimpson, main="Inverse Simpson index", xlab="")
shapiro.test(richness_TP$InvSimpson)
kruskal.test(richness_TP$InvSimpson ~ sample_data(TP)$mockT)


#Metamorphs
Meta <- subset_samples(offspring, GruopID == "Metamorph")
richness_Meta <- estimate_richness(Meta)

hist(richness_Meta$Shannon, main="Shannon index", xlab="")
shapiro.test(richness_Meta$Shannon)
kruskal.test(richness_Meta$Shannon ~ sample_data(Meta)$mockT)

hist(richness_Meta$Chao1, main="Chao1 index", xlab="")
shapiro.test(richness_Meta$Chao1)
kruskal.test(richness_Meta$Chao1 ~ sample_data(Meta)$mockT)

hist(richness_Meta$InvSimpson, main="Inverse Simpson index", xlab="")
shapiro.test(richness_Meta$InvSimpson)
kruskal.test(richness_Meta$InvSimpson ~ sample_data(Meta)$mockT)


# 5.- Intersection analyses ----
#### subset by treatment
ps_cpos <- subset_samples(database.pso, Posit_ctrl == "1")
ps_cneg <- subset_samples(database.pso, Neg_ctrl == "1")
ps_test <- subset_samples(database.pso, Test == "1")

#### UpsetR graph

#normal transport treatment

intersectionP1 <-
  upset_pq(
    ps_cpos,
    fact = "GruopID", width_ratio = 0.2,
    #  taxa_fill = "Phyllum", 
    min_nb_seq = 10,
    na_remove = TRUE,
    numeric_fonction = sum,
    rarefy_after_merging = FALSE,
    min_size = 7,
    
    guides='over',
    queries=list(
      upset_query(
        intersect = c('Father', 'Tadpole'),
        color='green4', fill='green4',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Substrate', 'Tadpole'),
        color='yellow3', fill='yellow3',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Female', 'Tadpole'),
        color='maroon', fill='maroon',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Water', 'Tadpole'),
        color='dodgerblue3', fill='dodgerblue3',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Father', 'Substrate'),
        color='red3', fill='red3',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Father', 'Metamorph'),
        color='green3', fill='green3',
        only_components = 'intersections_matrix',
      )),
    set_sizes=(
      upset_set_size(
        geom=geom_bar(),
        position='right'
      )))
intersectionP1

#no transport treatment

intersectionN1 <-
  upset_pq(
    ps_cneg,
    fact = "GruopID", width_ratio = 0.2,
    queries=list(
      upset_query(
        intersect = c('Father', 'Tadpole'),
        color='green4', fill='green4',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Substrate', 'Tadpole'),
        color='yellow3', fill='yellow3',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Female', 'Tadpole'),
        color='maroon', fill='maroon',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Water', 'Tadpole'),
        color='dodgerblue3', fill='dodgerblue3',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Father', 'Substrate'),
        color='red3', fill='red3',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Father', 'Metamorph'),
        color='green3', fill='green3',
        only_components = 'intersections_matrix',
      )),
    taxa_fill = "Phyllum",
    min_nb_seq = 5,
    na_remove = TRUE,
    numeric_fonction = sum,
    rarefy_after_merging = FALSE,
    min_size = 7,
    guides='over',
    set_sizes=(
      upset_set_size(
        geom=geom_bar(),
        position='right'
      )))
intersectionN1



#foster transport treatment

intersectionF1 <-
  upset_pq(
    ps_test,
    fact = "GruopID", width_ratio = 0.2,
    taxa_fill = "Phyllum",
    min_nb_seq = 5,
    na_remove = TRUE,
    numeric_fonction = sum,
    rarefy_after_merging = FALSE,
    min_size = 7,
    guides='over',
    queries=list(
      upset_query(
        intersect = c('Father', 'Tadpole'),
        color='green4', fill='green4',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Substrate', 'Tadpole'),
        color='yellow3', fill='yellow3',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Female', 'Tadpole'),
        color='maroon', fill='maroon',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Water', 'Tadpole'),
        color='dodgerblue3', fill='dodgerblue3',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Father', 'Substrate'),
        color='red3', fill='red3',
        only_components = 'intersections_matrix',
      ),
      upset_query(
        intersect = c('Father', 'Metamorph'),
        color='green3', fill='green3',
        only_components = 'intersections_matrix',
      )),
    set_sizes=(
      upset_set_size(
        geom=geom_bar(),
        position='right'
      )))
intersectionF1

# 6.- Extraction of acquired taxa by transportation ----

#For TADPOLES - Subset by trial of FOSTER TRANSPORT treatment
#s3
set3full <- subset_samples(database.pso, S3 == "1")
set3 <- subset_samples(set3full, Var == "var")

#extract non-shared taxa of foster parents (unique foster father taxa)
ns3 <- phyloseq_extract_non_shared_otus(set3, samp_names = c("A018"))

#extract shared taxa of foster father with tadpoles
set3tp <- subset_samples(set3full, GruopID == "Tadpole")
mergetp <- merge_phyloseq(ns3,set3tp)
s3_transferred <- phyloseq_extract_shared_otus(mergetp)

#extract names of taxa
names_s3_transferred <- tax_table(s3_transferred)

#s23
set23full <- subset_samples(database.pso, S23 == "1")
set23 <- subset_samples(set23full, Var == "var")
ns23 <- phyloseq_extract_non_shared_otus(set23, samp_names = c("A054"))
set23tp <- subset_samples(set23full, GruopID == "Tadpole")
mergetp23 <- merge_phyloseq(ns23,set23tp)
s23_transferred <- phyloseq_extract_shared_otus(mergetp23)
names_s23_transferred <- tax_table(s23_transferred)

#s24
set24full <- subset_samples(database.pso, S24 == "1")
set24 <- subset_samples(set24full, Var == "var")
ns24 <- phyloseq_extract_non_shared_otus(set24, samp_names = c("A002"))
set24tp <- subset_samples(set24full, GruopID == "Tadpole")
mergetp24 <- merge_phyloseq(ns24,set24tp)
s24_transferred <- phyloseq_extract_shared_otus(mergetp24)
names_s24_transferred <- tax_table(s24_transferred)

#s26
set26full <- subset_samples(database.pso, S26 == "1")
set26 <- subset_samples(set26full, Var == "var")
ns26 <- phyloseq_extract_non_shared_otus(set26, samp_names = c("A028"))
set26tp <- subset_samples(set26full, GruopID == "Tadpole")
mergetp26 <- merge_phyloseq(ns26,set26tp)
s26_transferred <- phyloseq_extract_shared_otus(mergetp26)
names_s26_transferred <- tax_table(s26_transferred)

#s52
set52full <- subset_samples(database.pso, S52 == "1")
set52 <- subset_samples(set52full, Var == "var")
ns52 <- phyloseq_extract_non_shared_otus(set52, samp_names = c("A014"))
set52tp <- subset_samples(set52full, GruopID == "Tadpole")
mergetp52 <- merge_phyloseq(ns52,set52tp)
s52_transferred <- phyloseq_extract_shared_otus(mergetp52)
names_s52_transferred <- tax_table(s52_transferred)

#For METAMORPHS - Subset by trial of FOSTER TRANSPORT treatment
#s3
set3m <- subset_samples(set3full, GruopID == "Metamorph")
mergem <- merge_phyloseq(ns3,set3m)
s3m_transferred <- phyloseq_extract_shared_otus(mergem)
names_s3m_transferred <- tax_table(s3m_transferred)

#s23
set23m <- subset_samples(set23full, GruopID == "Metamorph")
mergem23 <- merge_phyloseq(ns23,set23m)
s23m_transferred <- phyloseq_extract_shared_otus(mergem23)
names_s23m_transferred <- tax_table(s23m_transferred)

#s26
set26m <- subset_samples(set26full, GruopID == "Metamorph")
mergem26 <- merge_phyloseq(ns26,set26m)
s26m_transferred <- phyloseq_extract_shared_otus(mergem26)
names_s26m_transferred <- tax_table(s26m_transferred)

#s52
set52m <- subset_samples(set52full, GruopID == "Metamorph")
mergem52 <- merge_phyloseq(ns52,set52m)
s52m_transferred <- phyloseq_extract_shared_otus(mergem52)
names_s52m_transferred <- tax_table(s52m_transferred)


##For TADPOLES - Subset by trial of NORMAL TRANSPORT treatment
#s1
set1full <- subset_samples(database.pso, S1 == "1")
set1 <- subset_samples(set1full, Var == "var")
ns1 <- phyloseq_extract_non_shared_otus(set1, samp_names = c("A019"))
set1tp <- subset_samples(set1full, GruopID == "Tadpole")
mergetp <- merge_phyloseq(ns1,set1tp)
s1_transferred <- phyloseq_extract_shared_otus(mergetp)
names_s1_transferred <- tax_table(s1_transferred)

#s5
set5full <- subset_samples(database.pso, S5 == "1")
set5 <- subset_samples(set5full, Var == "var")
ns5 <- phyloseq_extract_non_shared_otus(set5, samp_names = c("A005"))
set5tp <- subset_samples(set5full, GruopID == "Tadpole")
mergetp5 <- merge_phyloseq(ns5,set5tp)
s5_transferred <- phyloseq_extract_shared_otus(mergetp5)
names_s5_transferred <- tax_table(s5_transferred)

#s9
set9full <- subset_samples(database.pso, S9 == "1")
set9 <- subset_samples(set9full, Var == "var")
ns9 <- phyloseq_extract_non_shared_otus(set9, samp_names = c("A035"))
set9tp <- subset_samples(set9full, GruopID == "Tadpole")
mergetp9 <- merge_phyloseq(ns9,set9tp)
s9_transferred <- phyloseq_extract_shared_otus(mergetp9)
names_s9_transferred <- tax_table(s9_transferred)

#s15
set15full <- subset_samples(database.pso, S15 == "1")
set15 <- subset_samples(set15full, Var == "var")
ns15 <- phyloseq_extract_non_shared_otus(set15, samp_names = c("A038"))
set15tp <- subset_samples(set15full, GruopID == "Tadpole")
mergetp15 <- merge_phyloseq(ns15,set15tp)
s15_transferred <- phyloseq_extract_shared_otus(mergetp15)
names_s15_transferred <- tax_table(s15_transferred)

#s19
set19full <- subset_samples(database.pso, S19 == "1")
set19 <- subset_samples(set19full, Var == "var")
ns19 <- phyloseq_extract_non_shared_otus(set19, samp_names = c("A056"))
set19tp <- subset_samples(set19full, GruopID == "Tadpole")
mergetp19 <- merge_phyloseq(ns19,set19tp)
s19_transferred <- phyloseq_extract_shared_otus(mergetp19)
names_s19_transferred <- tax_table(s19_transferred)

#s34
set34full <- subset_samples(database.pso, S34 == "1")
set34 <- subset_samples(set34full, Var == "var")
ns34 <- phyloseq_extract_non_shared_otus(set34, samp_names = c("A065b"))
set34tp <- subset_samples(set34full, GruopID == "Tadpole")
mergetp34 <- merge_phyloseq(ns34,set34tp)
s34_transferred <- phyloseq_extract_shared_otus(mergetp34)
names_s34_transferred <- tax_table(s34_transferred)

#s37
set37full <- subset_samples(database.pso, S37 == "1")
set37 <- subset_samples(set37full, Var == "var")
ns37 <- phyloseq_extract_non_shared_otus(set37, samp_names = c("A047a"))
set37tp <- subset_samples(set37full, GruopID == "Tadpole")
mergetp37 <- merge_phyloseq(ns37,set37tp)
s37_transferred <- phyloseq_extract_shared_otus(mergetp37)
names_s37_transferred <- tax_table(s37_transferred)

#s43
set43full <- subset_samples(database.pso, S43 == "1")
set43 <- subset_samples(set43full, Var == "var")
ns43 <- phyloseq_extract_non_shared_otus(set43, samp_names = c("A078"))
set43tp <- subset_samples(set43full, GruopID == "Tadpole")
mergetp43 <- merge_phyloseq(ns43,set43tp)
s43_transferred <- phyloseq_extract_shared_otus(mergetp43)
names_s43_transferred <- tax_table(s43_transferred)

#For METAMORPHS - Subset by trial of NORMAL TRANSPORT treatment
#s1m
set1mfull <- subset_samples(database.pso, S1 == "1")
set1m <- subset_samples(set1mfull, Var == "var")
ns1m <- phyloseq_extract_non_shared_otus(set1m, samp_names = c("A019"))
set1m <- subset_samples(set1full, GruopID == "Metamorph")
mergem <- merge_phyloseq(ns1m,set1m)
s1m_transferred <- phyloseq_extract_shared_otus(mergem)
names_s1m_transferred <- tax_table(s1m_transferred)

#s5
set5mfull <- subset_samples(database.pso, S5 == "1")
set5m <- subset_samples(set5mfull, Var == "var")
ns5m <- phyloseq_extract_non_shared_otus(set5m, samp_names = c("A005"))
set5m <- subset_samples(set5mfull, GruopID == "Metamorph")
mergem5 <- merge_phyloseq(ns5m,set5m)
s5m_transferred <- phyloseq_extract_shared_otus(mergem5)
names_s5m_transferred <- tax_table(s5m_transferred)

#s9
set9mfull <- subset_samples(database.pso, S9 == "1")
set9m <- subset_samples(set9mfull, Var == "var")
ns9m <- phyloseq_extract_non_shared_otus(set9m, samp_names = c("A035"))
set9m <- subset_samples(set9mfull, GruopID == "Metamorph")
mergem9 <- merge_phyloseq(ns9m,set9m)
s9m_transferred <- phyloseq_extract_shared_otus(mergem9)
names_s9m_transferred <- tax_table(s9m_transferred)

#s12
set12mfull <- subset_samples(database.pso, S12 == "1")
set12m <- subset_samples(set12mfull, Var == "var")
ns12m <- phyloseq_extract_non_shared_otus(set12m, samp_names = c("A030"))
set12m <- subset_samples(set12mfull, GruopID == "Metamorph")
mergem12 <- merge_phyloseq(ns12m,set12m)
s12m_transferred <- phyloseq_extract_shared_otus(mergem12)
names_s12m_transferred <- tax_table(s12m_transferred)

#s15
set15mfull <- subset_samples(database.pso, S15 == "1")
set15m <- subset_samples(set15mfull, Var == "var")
ns15m <- phyloseq_extract_non_shared_otus(set15m, samp_names = c("A038"))
set15m <- subset_samples(set15mfull, GruopID == "Metamorph")
mergem15 <- merge_phyloseq(ns15m,set15m)
s15m_transferred <- phyloseq_extract_shared_otus(mergem15)
names_s15m_transferred <- tax_table(s15m_transferred)

#s19
set19mfull <- subset_samples(database.pso, S19 == "1")
set19m <- subset_samples(set19mfull, Var == "var")
ns19m <- phyloseq_extract_non_shared_otus(set19m, samp_names = c("A056"))
set19m <- subset_samples(set19mfull, GruopID == "Metamorph")
mergem19 <- merge_phyloseq(ns19m,set19m)
s19m_transferred <- phyloseq_extract_shared_otus(mergem19)
names_s19m_transferred <- tax_table(s19m_transferred)

#s34
set34mfull <- subset_samples(database.pso, S34 == "1")
set34m <- subset_samples(set34mfull, Var == "var")
ns34m <- phyloseq_extract_non_shared_otus(set34m, samp_names = c("A065b"))
set34m <- subset_samples(set34mfull, GruopID == "Metamorph")
mergem34 <- merge_phyloseq(ns34m,set34m)
s34m_transferred <- phyloseq_extract_shared_otus(mergem34)
names_s34m_transferred <- tax_table(s34m_transferred)

#s37
set37mfull <- subset_samples(database.pso, S37 == "1")
set37m <- subset_samples(set37mfull, Var == "var")
ns37m <- phyloseq_extract_non_shared_otus(set37m, samp_names = c("A047a"))
set37m <- subset_samples(set37mfull, GruopID == "Metamorph")
merge37m <- merge_phyloseq(ns37m,set37m)
s37m_transferred <- phyloseq_extract_shared_otus(merge37m)
names_s37m_transferred <- tax_table(s37m_transferred)

#s43
set43mfull <- subset_samples(database.pso, S43 == "1")
set43m <- subset_samples(set43mfull, Var == "var")
ns43m <- phyloseq_extract_non_shared_otus(set43m, samp_names = c("A078"))
set43m <- subset_samples(set43mfull, GruopID == "Metamorph")
mergem43 <- merge_phyloseq(ns43m,set43m)
s43m_transferred <- phyloseq_extract_shared_otus(mergem43)
names_s43m_transferred <- tax_table(s43m_transferred)

#For TADPOLES - Subset by trial of NO TRANSPORT treatment
#s2
set2full <- subset_samples(database.pso, S2 == "1")
set2 <- subset_samples(set2full, Var == "var")
ns2 <- phyloseq_extract_non_shared_otus(set2, samp_names = c("A019"))
set2tp <- subset_samples(set2full, GruopID == "Tadpole")
mergetp <- merge_phyloseq(ns2,set2tp)
s2_transferred <- phyloseq_extract_shared_otus(mergetp)
names_s2_transferred <- tax_table(s2_transferred)

#s6
set6full <- subset_samples(database.pso, S6 == "1")
set6 <- subset_samples(set6full, Var == "var")
ns6 <- phyloseq_extract_non_shared_otus(set6, samp_names = c("A027"))
set6tp <- subset_samples(set6full, GruopID == "Tadpole")
mergetp6 <- merge_phyloseq(ns6,set6tp)
s6_transferred <- phyloseq_extract_shared_otus(mergetp6)
names_s6_transferred <- tax_table(s6_transferred)

#s10
set10full <- subset_samples(database.pso, S10 == "1")
set10 <- subset_samples(set10full, Var == "var")
ns10 <- phyloseq_extract_non_shared_otus(set10, samp_names = c("A035"))
set10tp <- subset_samples(set10full, GruopID == "Tadpole")
mergetp10 <- merge_phyloseq(ns10,set10tp)
s10_transferred <- phyloseq_extract_shared_otus(mergetp10)
names_s10_transferred <- tax_table(s10_transferred)

#s17
set17full <- subset_samples(database.pso, S17 == "1")
set17 <- subset_samples(set17full, Var == "var")
ns17 <- phyloseq_extract_non_shared_otus(set17, samp_names = c("A047a"))
set17tp <- subset_samples(set17full, GruopID == "Tadpole")
mergetp17 <- merge_phyloseq(ns17,set17tp)
s17_transferred <- phyloseq_extract_shared_otus(mergetp17)
names_s17_transferred <- tax_table(s17_transferred)

#s18
set18full <- subset_samples(database.pso, S18 == "1")
set18 <- subset_samples(set18full, Var == "var")
ns18 <- phyloseq_extract_non_shared_otus(set18, samp_names = c("A048"))
set18tp <- subset_samples(set18full, GruopID == "Tadpole")
mergetp18 <- merge_phyloseq(ns18,set18tp)
s18_transferred <- phyloseq_extract_shared_otus(mergetp18)
names_s18_transferred <- tax_table(s18_transferred)

#s35
set35full <- subset_samples(database.pso, S35 == "1")
set35 <- subset_samples(set35full, Var == "var")
ns35 <- phyloseq_extract_non_shared_otus(set35, samp_names = c("A069"))
set35tp <- subset_samples(set35full, GruopID == "Tadpole")
mergetp35 <- merge_phyloseq(ns35,set35tp)
s35_transferred <- phyloseq_extract_shared_otus(mergetp35)
names_s35_transferred <- tax_table(s35_transferred)

#s40
set40full <- subset_samples(database.pso, S40 == "1")
set40 <- subset_samples(set40full, Var == "var")
ns40 <- phyloseq_extract_non_shared_otus(set40, samp_names = c("A018"))
set40tp <- subset_samples(set40full, GruopID == "Tadpole")
mergetp40 <- merge_phyloseq(ns40,set40tp)
s40_transferred <- phyloseq_extract_shared_otus(mergetp40)
names_s40_transferred <- tax_table(s40_transferred)

#s42
set42full <- subset_samples(database.pso, S42 == "1")
set42 <- subset_samples(set42full, Var == "var")
ns42 <- phyloseq_extract_non_shared_otus(set42, samp_names = c("S303"))
set42tp <- subset_samples(set42full, GruopID == "Tadpole")
mergetp42 <- merge_phyloseq(ns42,set42tp)
s42_transferred <- phyloseq_extract_shared_otus(mergetp42)
names_s42_transferred <- tax_table(s42_transferred)

#For METAMORPHS - Subset by trial of NO TRANSPORT treatment
#s2
set2full <- subset_samples(database.pso, S2 == "1")
set2 <- subset_samples(set2full, Var == "var")
ns2 <- phyloseq_extract_non_shared_otus(set2, samp_names = c("A019"))
set2m <- subset_samples(set2full, GruopID == "Metamorph")
mergem <- merge_phyloseq(ns2,set2m)
s2m_transferred <- phyloseq_extract_shared_otus(mergem)
names_s2m_transferred <- tax_table(s2m_transferred)

#s6
set6full <- subset_samples(database.pso, S6 == "1")
set6 <- subset_samples(set6full, Var == "var")
ns6 <- phyloseq_extract_non_shared_otus(set6, samp_names = c("A027"))
set6m <- subset_samples(set6full, GruopID == "Metamorph")
mergem6 <- merge_phyloseq(ns6,set6m)
s6m_transferred <- phyloseq_extract_shared_otus(mergem6)
names_s6m_transferred <- tax_table(s6m_transferred)

#s10
set10full <- subset_samples(database.pso, S10 == "1")
set10 <- subset_samples(set10full, Var == "var")
ns10 <- phyloseq_extract_non_shared_otus(set10, samp_names = c("A035"))
set10m <- subset_samples(set10full, GruopID == "Metamorph")
mergem10 <- merge_phyloseq(ns10,set10m)
s10m_transferred <- phyloseq_extract_shared_otus(mergem10)
names_s10m_transferred <- tax_table(s10m_transferred)

#s16
set16full <- subset_samples(database.pso, S16 == "1")
set16 <- subset_samples(set16full, Var == "var")
ns16 <- phyloseq_extract_non_shared_otus(set16, samp_names = c("A047a"))
set16m <- subset_samples(set16full, GruopID == "Metamorph")
mergem16 <- merge_phyloseq(ns16,set16m)
s16m_transferred <- phyloseq_extract_shared_otus(mergem16)
names_s16m_transferred <- tax_table(s16m_transferred)

#s17
set17full <- subset_samples(database.pso, S17 == "1")
set17 <- subset_samples(set17full, Var == "var")
ns17 <- phyloseq_extract_non_shared_otus(set17, samp_names = c("A047a"))
set17m <- subset_samples(set17full, GruopID == "Metamorph")
mergem17 <- merge_phyloseq(ns17,set17m)
s17m_transferred <- phyloseq_extract_shared_otus(mergem17)
names_s17m_transferred <- tax_table(s17m_transferred)

#s35
set35full <- subset_samples(database.pso, S35 == "1")
set35 <- subset_samples(set35full, Var == "var")
ns35 <- phyloseq_extract_non_shared_otus(set35, samp_names = c("A069"))
set35m <- subset_samples(set35full, GruopID == "Metamorph")
mergem35 <- merge_phyloseq(ns35,set35m)
s35m_transferred <- phyloseq_extract_shared_otus(mergem35)
names_s35m_transferred <- tax_table(s35m_transferred)

#s40
set40full <- subset_samples(database.pso, S40 == "1")
set40 <- subset_samples(set40full, Var == "var")
ns40 <- phyloseq_extract_non_shared_otus(set40, samp_names = c("A018"))
set40m <- subset_samples(set40full, GruopID == "Metamorph")
mergem40 <- merge_phyloseq(ns40,set40m)
s40m_transferred <- phyloseq_extract_shared_otus(mergem40)
names_s40m_transferred <- tax_table(s40m_transferred)

#s42
set42full <- subset_samples(database.pso, S42 == "1")
set42 <- subset_samples(set42full, Var == "var")
ns42 <- phyloseq_extract_non_shared_otus(set42, samp_names = c("S303"))
set42m <- subset_samples(set42full, GruopID == "Metamorph")
mergem42 <- merge_phyloseq(ns42,set42m)
s42m_transferred <- phyloseq_extract_shared_otus(mergem42)
names_s42m_transferred <- tax_table(s42m_transferred)


# 7.- Analyses of proportion of contribution ----
# Obtained data from SourceTracker to use those proportions for following analyses

#7.1.- model 1
#Load data with proportions of contribution by all source types
Q1 <- read.table("Contributions_sources.txt", h=T)
Q1$Source <- as.factor(Q1$Source) 

#model for testing differences in contribution by type of source
Q1lmm2 <- glmmTMB(contribution ~ Treatment * Source * Sink  + (1|Trial),
                  data = Q1, family=beta_family())

glmmTMB:::Anova.glmmTMB(Q1lmm2)
emmeans(Q1lmm2, pairwise ~ Treatment | Source | Sink)

#diagnostics of model
testDispersion(simulateResiduals(Q1lmm2, plot = TRUE))

#7.2.- model 2
#Load data with proportions of contribution from developmental steps
devo   <- read.table("contributions_development.txt",h=T)

#model for testing differences in contribution from developmental stages
Q2m <- glmmTMB(Source_pr ~ Treatment + Sink  + (1|Trial),
               data = devo, family=beta_family() )

glmmTMB:::Anova.glmmTMB(Q2m)
emmeans(Q2m, pairwise ~  Sink)

simulateResiduals(Q2m, plot = TRUE)

#7.3.- model 3
#Load data with number of putative Bd inhibitory zOTUs acquired from fathers
transfer <- read.table("transferred_zOTU_model.txt", h=T)
transfer$condition <- as.factor(transfer$condition) 
transfer$sample_type <- as.factor(transfer$sample_type) 
transfer$Trial <- as.factor(transfer$Trial) 

#subset by offspring type
transferTP <- transfer[ which(transfer$sample_type=='tadpole'), ]
transferM <- transfer[ which(transfer$sample_type=='metamorph'), ]

#models for testing differences in number of inhibitor Bd taxa in tadpoles and metamorphs
#tadpoles
TP <- glm(zOTU_putative_Bd ~ condition, data = transferTP, family=poisson)

Anova(TP)
emmeans(TP, pairwise ~ condition)

simulateResiduals(TP, plot = TRUE)

#metamorphs
M <- glmmTMB(zOTU_putative_Bd ~ condition, data = transferM, 
             family= nbinom2())

glmmTMB:::Anova.glmmTMB(M)
emmeans(M, pairwise ~ condition)

testDispersion(simulateResiduals(M, plot = TRUE))
