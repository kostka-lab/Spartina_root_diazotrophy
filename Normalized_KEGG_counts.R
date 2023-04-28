#Set working directory
setwd("~/")


#Function to make normalized count reads mapping from KO numbers for metaG and metaT short reads.
#Folder containing annotation files has to have same name as eggnog-mapper annotation files

run_short_read_eggnog <- function(folder_name){
  library(dplyr)
  folder_direction <- paste0("./",folder_name)
  setwd(folder_direction)
  eggnog_df_1 <- read.delim2(paste0(folder_name, "_L001.emapper.annotations"), skip = 4, header = T, sep = "\t")
  eggnog_df_1$X.query <- gsub("_0", "_1", eggnog_df_1$X.query)
  eggnog_df_2 <- read.delim2(paste0(folder_name, "_L002.emapper.annotations"), skip = 4, header = T, sep = "\t")
  names(eggnog_df_2) <- names(eggnog_df_1)
  eggnog_df_2$X.query <- gsub("_0", "_2", eggnog_df_2$X.query)
  eggnog_df <- rbind(eggnog_df_1, eggnog_df_2)
  rm(list = c("eggnog_df_1", "eggnog_df_2"))
  eggnog_df$OG_1 <- sub("\\@.*", "", eggnog_df$eggNOG_OGs)
  eggnog_df_euk <- eggnog_df[grep("Eukaryota", eggnog_df$eggNOG_OGs),]
  eggnog_df_bac <- eggnog_df[grep("Bacteria", eggnog_df$eggNOG_OGs),]
  eggnog_df_arch <- eggnog_df[grep("Archaea", eggnog_df$eggNOG_OGs),]
  eggnog_df_prok <- rbind(eggnog_df_bac, eggnog_df_arch)
  eggnog_df_virus <- eggnog_df[grep("Viruses", eggnog_df$eggNOG_OGs),]
  
  Domain <- c("Total", "Eukaryota", "Bacteria", "Archaea", "Prokaryota", "Virus") 
  Total_reads <- c(length(eggnog_df$X.query), 
                   length(eggnog_df_euk$X.query),
                   length(eggnog_df_bac$X.query), 
                   length(eggnog_df_arch$X.query),
                   length(eggnog_df_prok$X.query),
                   length(eggnog_df_virus$X.query))
  
  Domain_reads_OG_df <- tibble(Domain = Domain, Total_reads =Total_reads)
  Domain_reads_OG_df$Percentage <- 100*Domain_reads_OG_df$Total_reads/length(eggnog_df$X.query)
  
  #Write table calculating the percentage of Eukaryotic, Prokaryotic and Viral short reads
  write.table(x = Domain_reads_OG_df, paste0(folder_name, "_Kingdom_reads.tsv"), sep = "\t", 
              row.names = F, quote = F)
  
  #Make count profiles based on eggnog OG and Kegg KO numbers
  OG_profile <- as.data.frame(sort(table(eggnog_df_prok$OG_1)))
  KEGG_profile <- as.data.frame(sort(table(eggnog_df_prok$KEGG_ko)))
  
  
  #Calculate median count of 10 universal single-copy marker genes for normalization 
  MG_OGs <- OG_profile[grep("^COG0012|^COG0016|^COG0018|^COG0172|^COG0215|^COG0495|^COG0525|^COG0533|^COG0541|^COG0552", OG_profile$Var1),]
  median_MG_OGs <- median(MG_OGs$Freq)
  #Normalize eggnog OG count profile
  OG_profile$Freq_normalized <- OG_profile$Freq/median_MG_OGs
  write.table(MG_OGs, paste0(folder_name, "OG_MGs_sum.tsv"), sep = "\t") #Export the count of marker genes by eggnog OG
  
  #Calculate median count of 10 universal single-copy marker genes for normalization 
  MG_KEGG <- KEGG_profile[grep("^ko:K06942$|^ko:K01889$|^ko:K01887$|^ko:K01875$|^ko:K01883$|^ko:K01869$|^ko:K01873$|^ko:K01409$|^ko:K03106$|^ko:K03110$", KEGG_profile$Var1),]
  median_MG_KEGG <- median(MG_KEGG$Freq)
  #Normalize eggnog KEGG count profile
  KEGG_profile$Freq_normalized <- KEGG_profile$Freq/median_MG_KEGG
  write.table(MG_KEGG, paste0(folder_name, "_MGs_sum.tsv"), sep = "\t")  #Export the count of marker genes by KEGG KO
  
  names(KEGG_profile) <- c("KEGG_ko", "Count", "Count_normalized")
  
  KEGG_profile$KEGG_ko <- gsub("^-$", "Undetermined", KEGG_profile$KEGG_ko)
  
  #Merge OGs counts with the same KO number
  KEGG_profile_unique <-  eggnog_df_prok[match(unique(eggnog_df_prok$KEGG_ko), eggnog_df_prok$KEGG_ko),]
  KEGG_profile_unique <- KEGG_profile_unique[,c(9,11,12,13,14,15,16,8)]
  KEGG_profile <- merge(KEGG_profile, KEGG_profile_unique, by = "KEGG_ko")
  
  #Write normalized count database with annotation
  write.table(KEGG_profile, paste0(folder_name, "_KEGG_profile.tsv"), sep = "\t")
  setwd("../")
  rm(list = ls())
  gc()
}