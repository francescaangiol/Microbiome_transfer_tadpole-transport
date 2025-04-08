# Microbiome_transfer_tadpole-transport
R code for analysis linked to the project (research article) "Paternal care as a source of key skin microbiota in a poison frog".

Description of the data:
This data was collected for a microbiome tracking experiment where we describe the relevance of verticla transmission during tadpole transportation behaviours (parental care performed by the adult male). We set three treatments (normal transport, no transport and foster transport) to determine the effect of transportation on the seeding of skin microbiome of the offspring. 

Files:

-METADATA_GOOD.txt = clean sample metadata table (used for downstream analyses)
-OTU_tax_GOOD.txt = clean tax table assigned to zOTUs found in our experiment (used for dwonstream analyses)
-ZOTU_c97_Count_GOOD.txt = clean otu table (used for dwonstream analyses)
-initial_metadata.txt = raw sample metadata table
-initial_otu_table.txt = raw otu table
-initial_tax_table.txt = raw tax table
-Contributions_sources.txt = dataframe with mean proportion of contribution by all sources to all sinks (obtained from SourceTracker)
-contributions_development.txt = dataframe with mean proportion of contribution by developmental stages (obtained from SourceTracker)
-transferred_zOTU_model.txt = dataframe of the putative Bd inhibitory taxa from our experimental setup
-Microbiome_transfer_script.R = r code used for analyses


Description: Code used in RStudio for analyzing the data
Code/software
Software used: R

All packages are mentioned within the code
