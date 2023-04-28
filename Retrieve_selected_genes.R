library(dplyr)

#Set working directory
setwd("~/p-jk129-0/rich_project_bio-konstantinidis/Sapelo_2020_MetaT/eggnog_short_reads/")

selected_genes <- c('ko:K02586', 'ko:K02591', 'ko:K02588', 'ko:K01601', 
           'ko:K01602 ', 'ko:K11180', 'ko:K11181', 'ko:K00394', 
           'ko:K00395', 'ko:K00958', 'ko:K17222', 'ko:K17224', 
           'ko:K17225', 'ko:K17223', 'ko:K22622', 'ko:K17226', 
           'ko:K17227', 'ko:K10944', 'ko:K10945', 'ko:K10946', 
           'ko:K10535', 'ko:K00370', 'ko:K00371', 'ko:K00374', 
           'ko:K02567', 'ko:K02568', 'ko:K00368', 'ko:K15864', 
           'ko:K04561', 'ko:K02305', 'ko:K00376', 'ko:K00362', 
           'ko:K00363', 'ko:K03385', 'ko:K15876', 'ko:K02297', 
           'ko:K02298', 'ko:K02299', 'ko:K02300', 'ko:K00425', 
           'ko:K00426', 'ko:K00424', 'ko:K22501', 'ko:K02275', 
           'ko:K02274', 'ko:K02276', 'ko:K15408', 'ko:K02277', 
           'ko:K00404', 'ko:K00405', 'ko:K15862', 'ko:K00407', 
           'ko:K00406', 'ko:K17067', 'ko:K00399', 'ko:K00195', 
           'ko:K05299', 'ko:K15022', 'ko:K22015', 'ko:K25123', 
           'ko:K25124', 'ko:K01938', 'ko:K01491', 'ko:K01500', 
           'ko:K00297', 'ko:K25007', 'ko:K25008', 'ko:K15023', 
           'ko:K14138', 'ko:K00197', 'ko:K00194')

#Function to retrieve reads mapping  to KO numbers. Folder containing annotation files 
#has to have same name as eggnog-mapper annotation files

run_eggnog_selected_genes <- function(folder_name){
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
  eggnog_df_bac <- eggnog_df[grep("Bacteria", eggnog_df$eggNOG_OGs),]
  eggnog_df_arch <- eggnog_df[grep("Archaea", eggnog_df$eggNOG_OGs),]
  eggnog_df_prok <- rbind(eggnog_df_bac, eggnog_df_arch)

  eggnog_df_genes <- eggnog_df[eggnog_df$KEGG_ko %in% selected_genes,]

  write.table(eggnog_df_genes, paste0(folder_name, "_selected_genes.tsv"), sep = "\t")
  setwd("../")
  rm(list = ls())
  gc()
}


