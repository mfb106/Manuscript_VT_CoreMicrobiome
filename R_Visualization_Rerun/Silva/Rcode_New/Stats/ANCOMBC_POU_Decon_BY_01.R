# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ANCOMBC")
#BiocManager::install("mia")

rm(list = ls())
library(mia)
library(ANCOMBC)
library(tidyverse)
library(DT)
library(devtools)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(lubridate)
library(broom)
library(lmtest)
library(car)
library(xlsx)
library(tidyverse)
library(qiime2R)
library(reshape2)
library(readxl)
library(dplyr)
library(openxlsx)
library(DataCombine)
library(ggplot2)
library(tidyr)
library(vegan)
library(tidyverse)
library(ape)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(phangorn)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)
library(phyloseq)
library(plyr)
library(hrbrthemes)
library(viridis)
library(gridExtra)
library(shadowtext)
library(ggpubr)
library(gghighlight)
library(ANCOMBC)
library(tidyverse)
library(DT)
setwd("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures")
getwd()

#Theme 
#Themes
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1.2)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            #axis.text.x = element_text(angle = 45),
            axis.text = element_text(size=rel(1.2)), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.box="vertical", 
            legend.margin=margin(),
            legend.direction = "horizontal",
            legend.key.size= unit(2, "cm"),
            legend.text = element_text(size=rel(1.2)),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic",size=rel(1.2)),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#7fc97f","#fdb462","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#7fc97f","#fdb462","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

####
####None rarified inputs
#from final table that was created... note this is not the rarefied data, if rarification is need can use vegan package - https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/rarefy
table_all<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Removed_classsilva_clust99_Comparison__table-no-mitochondria-chloroplast.qza")

#taxonomy, classification
taxonomy_all<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Classify_silva_clust99_Comparison_classifications.qza")

#alpha diversity evenness
evenness_RL<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/betadiversity_3500-results/evenness_vector.qza")

#alpha diversity faithpd
faith_RL<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/phylo-betadiversity_3500-results/faith_pd_vector.qza")

#alpha diversity observed otus
observedotus_RL<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/phylo-betadiversity_3500-results/observed_features_vector.qza")

#alpha diversity shannon
shannon_RL<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/betadiversity_3500-results/shannon_vector.qza")

#beta diversity weighted and unweighted unifrac
NMDS_unweighted_unifrac_RL<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/phylo-betadiversity_3500-results/unweighted_unifrac_NMDS_results.qza")

NMDS_weighted_unifrac_RL<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/phylo-betadiversity_3500-results/weighted_unifrac_NMDS_results.qza")

#rooted tree
comp_rooted_tree<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/Comparison_Rerun_silva-rooted-tree.qza")

#Alpha diversity
alpha_diversity<-read_qza("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/Qiime_Outputs_Rerun/Silva/Analysis/alpha-diversity/alpha_diversity.qza")


#for all comparisons 
rownames(metadata_raw) <- metadata_raw$SampleID
metadata_raw$SampleID <- NULL

#extract data version, one for RA and one for PA 
table_all_data<-as.data.frame(table_all$data)
taxonomy_all_data<-taxonomy_all$data

#Filter OTUs from Deontam
#Decontam
Decontam_OTU<-read.xlsx("C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/QAQC/Decontam_all.xlsx", sheet= 'Sheet1')
Decontam_OTU_filter<-Decontam_OTU[Decontam_OTU$contaminant == 'TRUE',]#can add Blank

table_all_data <- tibble::rownames_to_column(table_all_data, "OTU")
table_all_data <- table_all_data[!(table_all_data$OTU %in% Decontam_OTU_filter$OTU),]
rownames(table_all_data) <- table_all_data[,1]
table_all_data[,1] <- NULL

table_all_data_PA<-table_all_data
table_all_data_PA[]<-as.integer((table_all_data_PA[]!=0))

#parse out taxonomy data so that all columns are individual 
taxonomy_all_data<-parse_taxonomy(taxonomy_all_data)

#merge taxonomy into data table 
Combined_Tax_Data<-transform(merge(taxonomy_all_data,table_all_data,by=0), row.names=Row.names, Row.names=NULL)
Combined_Tax_Data_PA<-transform(merge(taxonomy_all_data,table_all_data_PA,by=0), row.names=Row.names, Row.names=NULL)

#Removed Outliers, note could maybe provide a normalization especially for PR1WEEF
Test.ConnorPrimer58<-as.data.frame(table_all_data$ConnerPrimer58)  ###Removed Connorprimer58 due to low number of reads (3) 

Test.PR1WEEF<-as.data.frame(table_all_data$PR1WEEF)  ###Removed PR1WEEF due to High number of reads (256687) 
sum(Test.PR1WEEF$`table_all_data$PR1WEEF`)


#Read Meta Data 
Metadata_POCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POCPOU[Metadata_POCPOU$Sample_Local == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POCPOU[Metadata_POCPOU$Sample_Local == "POU",]

#Order Meta Data 
Metadata_POCPOU$Classification_1 = factor(Metadata_POCPOU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POU
Meta.All.POU<-Metadata_POU
Meta.All.POU<-Meta.All.POU[Meta.All.POU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POU<-Meta.All.POU[Meta.All.POU$Description != 'PR1WEEF',] #removed because of high read count 

Meta.All.POU.Filter<-Meta.All.POU[Meta.All.POU$Classification_1 == 'Potable'| Meta.All.POU$Classification_1 =='Non-potable Reuse'|Meta.All.POU$Classification_1 =='Potable Reuse'|Meta.All.POU$Classification_1 =='',]#can add Blank
#Meta.All.POU.Filter<-Meta.All.POU.Filter[Meta.All.POU.Filter$Matrix == 'water'| Meta.All.POU.Filter$Matrix == 'EXTRACTIONBLANK'| Meta.All.POU.Filter$Matrix == 'FIELDBLANK'| Meta.All.POU.Filter$Matrix == 'PCRBLANK'| Meta.All.POU.Filter$Matrix == 'NA',]
Meta.All.POU.RowNames<-as.data.frame(rownames(Meta.All.POU.Filter))
Meta.All.POU.List<-dplyr::pull(Meta.All.POU.RowNames,1)

OTU.Table.Master<-table_all_data
OTU.Table.Master.Trans<- as.data.frame(t(as.matrix(OTU.Table.Master)))
OTU.All.POU.Filter<- OTU.Table.Master.Trans %>% filter(row.names(OTU.Table.Master.Trans) %in% Meta.All.POU.List)
Meta.All.POU.Filter<-Meta.All.POU.Filter %>% filter(row.names(Meta.All.POU.Filter) %in% row.names(OTU.All.POU.Filter))

Tax.Master<-taxonomy_all$data
Tax.Master<-parse_taxonomy(Tax.Master)

OTU.phy.All.POU=otu_table(as.matrix(OTU.All.POU.Filter), taxa_are_rows=FALSE)
Tax.phy.All.POU = tax_table(as.matrix(Tax.Master))
Meta.phy.All.POU = sample_data(Meta.All.POU.Filter)

physeq.All.POU = phyloseq(OTU.phy.All.POU, Tax.phy.All.POU, Meta.phy.All.POU)

Tree.Master<-comp_rooted_tree$data
edges=phy_tree(Tree.Master)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has

Tree.Master.V2 <- ape::multi2di(Tree.Master)
edges=phy_tree(Tree.Master.V2)$edge
mycounts = table(edges[,1]) # Source nodes; 1st column of edge matrix
length(mycounts[mycounts ==2]) # Number of nodes with exactly 2 children
length(mycounts[mycounts !=2]) # Number of nodes with more or fewer children
mycounts[mycounts !=2] # How many nodes each of the above has


physeq.combined.All.POU = merge_phyloseq(physeq.All.POU,Tree.Master.V2)
physeq.combined.All.POU

#remove if needed
phy_tree(physeq.combined.All.POU)<-phangorn::midpoint(phy_tree(physeq.combined.All.POU))

physeq.combined.All.POU.rarified3500<-rarefy_even_depth(physeq.combined.All.POU, sample.size = 3500,
                                                        rngseed = 711, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)



##Tutorials 
# https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
# https://microbiome.github.io/course_2021_radboud/differential-abundance-analysis.html
# https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html

#ancombc also supports importing data in phyloseq format, Class_1 on POU

physeq.All.POU_level2_Class1_Pot_NonPot_alt <- tax_glom(physeq.All.POU, taxrank=rank_names(physeq.All.POU)[2], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

physeq.All.POU_level2_Class1_Pot_NonPot_alt@sam_data$Class1_Pot_NonPot = factor(physeq.All.POU_level2_Class1_Pot_NonPot_alt@sam_data$Class1_Pot_NonPot, levels = c("Potable", "Non-potable Reuse"))

OUT.physeq.All.POU_level2_Class1_Pot_NonPot_alt = ancombc(data = NULL, assay_name = NULL,
              tax_level = "Phylum", phyloseq = physeq.All.POU_level2_Class1_Pot_NonPot_alt,
              formula = "Class1_Pot_NonPot",
              p_adj_method = "BY", prv_cut = 0.01, lib_cut = 0,
              group = "Class1_Pot_NonPot", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
              max_iter = 100, conserve = TRUE, alpha = 0.01, global = TRUE,
              n_cl = 1, verbose = TRUE)


res.physeq.All.POU_level2_Class1_Pot_NonPot_alt = OUT.physeq.All.POU_level2_Class1_Pot_NonPot_alt$res
res_global.physeq.All.POU_level2_Class1_Pot_NonPot_alt = OUT.physeq.All.POU_level2_Class1_Pot_NonPot_alt$res_global



##Genus, level6
physeq.All.POU_level6_Class1_Pot_NonPot_alt <- tax_glom(physeq.All.POU, taxrank=rank_names(physeq.All.POU)[6], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

physeq.All.POU_level6_Class1_Pot_NonPot_alt@sam_data$Class1_Pot_NonPot = factor(physeq.All.POU_level6_Class1_Pot_NonPot_alt@sam_data$Class1_Pot_NonPot, levels = c("Potable", "Non-potable Reuse"))

OUT.physeq.All.POU_level6_Class1_Pot_NonPot_alt = ancombc(data = NULL, assay_name = NULL,
                                                             tax_level = "Genus", phyloseq = physeq.All.POU_level6_Class1_Pot_NonPot_alt,
                                                             formula = "Class1_Pot_NonPot",
                                                             p_adj_method = "BY", prv_cut = 0.01, lib_cut = 0,
                                                             group = "Class1_Pot_NonPot", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                                                             max_iter = 100, conserve = TRUE, alpha = 0.01, global = TRUE,
                                                             n_cl = 1, verbose = TRUE)


res.physeq.All.POU_level6_Class1_Pot_NonPot_alt = OUT.physeq.All.POU_level6_Class1_Pot_NonPot_alt$res
res_global.physeq.All.POU_level6_Class1_Pot_NonPot_alt = OUT.physeq.All.POU_level6_Class1_Pot_NonPot_alt$res_global


##OTU
physeq.All.POU_OTU_Class1_Pot_NonPot_alt <- physeq.All.POU

physeq.All.POU_OTU_Class1_Pot_NonPot_alt@sam_data$Class1_Pot_NonPot = factor(physeq.All.POU_OTU_Class1_Pot_NonPot_alt@sam_data$Class1_Pot_NonPot, levels = c("Potable", "Non-potable Reuse"))

OUT.physeq.All.POU_OTU_Class1_Pot_NonPot_alt = ancombc(data = NULL, assay_name = NULL,
                                                             tax_level = NULL, phyloseq = physeq.All.POU_OTU_Class1_Pot_NonPot_alt,
                                                             formula = "Class1_Pot_NonPot",
                                                             p_adj_method = "BY", prv_cut = 0.01, lib_cut = 0,
                                                             group = "Class1_Pot_NonPot", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                                                             max_iter = 100, conserve = TRUE, alpha = 0.01, global = TRUE,
                                                             n_cl = 1, verbose = TRUE)


res.physeq.All.POU_OTU_Class1_Pot_NonPot_alt = OUT.physeq.All.POU_OTU_Class1_Pot_NonPot_alt$res
res_global.physeq.All.POU_OTU_Class1_Pot_NonPot_alt = OUT.physeq.All.POU_OTU_Class1_Pot_NonPot_alt$res_global



###Classification_1
#ancombc also supports importing data in phyloseq format, Class_1 on POU
physeq.All.POU_level2_Classification_1_alt <- tax_glom(physeq.All.POU, taxrank=rank_names(physeq.All.POU)[2], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

physeq.All.POU_level2_Classification_1_alt@sam_data$Classification_1 = factor(physeq.All.POU_level2_Classification_1_alt@sam_data$Classification_1, levels = c("Potable","Potable Reuse", "Non-potable Reuse"))

OUT.physeq.All.POU_level2_Classification_1_alt = ancombc(data = NULL, assay_name = NULL,
                                                             tax_level = "Phylum", phyloseq = physeq.All.POU_level2_Classification_1_alt,
                                                             formula = "Classification_1",
                                                             p_adj_method = "BY", prv_cut = 0.01, lib_cut = 0,
                                                             group = "Classification_1", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                                                             max_iter = 100, conserve = TRUE, alpha = 0.01, global = TRUE,
                                                             n_cl = 1, verbose = TRUE)


res.physeq.All.POU_level2_Classification_1_alt = OUT.physeq.All.POU_level2_Classification_1_alt$res
res_global.physeq.All.POU_level2_Classification_1_alt = OUT.physeq.All.POU_level2_Classification_1_alt$res_global



##Genus, level6
physeq.All.POU_level6_Classification_1_alt <- tax_glom(physeq.All.POU, taxrank=rank_names(physeq.All.POU)[6], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

physeq.All.POU_level6_Classification_1_alt@sam_data$Classification_1 = factor(physeq.All.POU_level6_Classification_1_alt@sam_data$Classification_1, levels = c("Potable","Potable Reuse","Non-potable Reuse"))

OUT.physeq.All.POU_level6_Classification_1_alt = ancombc(data = NULL, assay_name = NULL,
                                                             tax_level = "Genus", phyloseq = physeq.All.POU_level6_Classification_1_alt,
                                                             formula = "Classification_1",
                                                             p_adj_method = "BY", prv_cut = 0.01, lib_cut = 0,
                                                             group = "Classification_1", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                                                             max_iter = 100, conserve = TRUE, alpha = 0.01, global = TRUE,
                                                             n_cl = 1, verbose = TRUE)


res.physeq.All.POU_level6_Classification_1_alt = OUT.physeq.All.POU_level6_Classification_1_alt$res
res_global.physeq.All.POU_level6_Classification_1_alt = OUT.physeq.All.POU_level6_Classification_1_alt$res_global


##OTU
physeq.All.POU_OTU_Classification_1_alt <- physeq.All.POU

physeq.All.POU_OTU_Classification_1_alt@sam_data$Classification_1 = factor(physeq.All.POU_OTU_Classification_1_alt@sam_data$Classification_1, levels = c("Potable", "Non-potable Reuse"))

OUT.physeq.All.POU_OTU_Classification_1_alt = ancombc(data = NULL, assay_name = NULL,
                                                          tax_level = NULL, phyloseq = physeq.All.POU_OTU_Classification_1_alt,
                                                          formula = "Classification_1",
                                                          p_adj_method = "BY", prv_cut = 0.01, lib_cut = 0,
                                                          group = "Classification_1", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                                                          max_iter = 100, conserve = TRUE, alpha = 0.01, global = TRUE,
                                                          n_cl = 1, verbose = TRUE)
res.physeq.All.POU_OTU_Classification_1_alt = OUT.physeq.All.POU_OTU_Classification_1_alt$res
res_global.physeq.All.POU_OTU_Classification_1_alt = OUT.physeq.All.POU_OTU_Classification_1_alt$res_global


###Creating Dataframe for output. 

Final_All.POU_level2_Class1_Pot_NonPot <- res.physeq.All.POU_level2_Class1_Pot_NonPot_alt$diff_abn %>%
  dplyr::full_join(res.physeq.All.POU_level2_Class1_Pot_NonPot_alt$q_val, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level2_Class1_Pot_NonPot_alt$p_val, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level2_Class1_Pot_NonPot_alt$lfc, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level2_Class1_Pot_NonPot_alt$se, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level2_Class1_Pot_NonPot_alt$W, by='taxon')

Final_All.POU_level6_Class1_Pot_NonPot <- res.physeq.All.POU_level6_Class1_Pot_NonPot_alt$diff_abn %>%
  dplyr::full_join(res.physeq.All.POU_level6_Class1_Pot_NonPot_alt$q_val, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level6_Class1_Pot_NonPot_alt$p_val, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level6_Class1_Pot_NonPot_alt$lfc, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level6_Class1_Pot_NonPot_alt$se, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level6_Class1_Pot_NonPot_alt$W, by='taxon')

Final_All.POU_OTU_Class1_Pot_NonPot <- res.physeq.All.POU_OTU_Class1_Pot_NonPot_alt$diff_abn %>%
  dplyr::full_join(res.physeq.All.POU_OTU_Class1_Pot_NonPot_alt$q_val, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_OTU_Class1_Pot_NonPot_alt$p_val, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_OTU_Class1_Pot_NonPot_alt$lfc, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_OTU_Class1_Pot_NonPot_alt$se, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_OTU_Class1_Pot_NonPot_alt$W, by='taxon')

### Classification_1

Final_All.POU_level2_Classification_1 <- res.physeq.All.POU_level2_Classification_1_alt$diff_abn %>%
  dplyr::full_join(res.physeq.All.POU_level2_Classification_1_alt$q_val, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level2_Classification_1_alt$p_val, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level2_Classification_1_alt$lfc, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level2_Classification_1_alt$se, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level2_Classification_1_alt$W, by='taxon')

Final_All.POU_level6_Classification_1 <- res.physeq.All.POU_level6_Classification_1_alt$diff_abn %>%
  dplyr::full_join(res.physeq.All.POU_level6_Classification_1_alt$q_val, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level6_Classification_1_alt$p_val, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level6_Classification_1_alt$lfc, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level6_Classification_1_alt$se, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_level6_Classification_1_alt$W, by='taxon')

Final_All.POU_OTU_Classification_1 <- res.physeq.All.POU_OTU_Classification_1_alt$diff_abn %>%
  dplyr::full_join(res.physeq.All.POU_OTU_Classification_1_alt$q_val, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_OTU_Classification_1_alt$p_val, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_OTU_Classification_1_alt$lfc, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_OTU_Classification_1_alt$se, by='taxon') %>%
  dplyr::full_join(res.physeq.All.POU_OTU_Classification_1_alt$W, by='taxon')


#Filter for Adjusted-Pvalue < 0.05, Play with outputs first. 


#Export
write_xlsx(Final_All.POU_level2_Class1_Pot_NonPot,"C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Stats/ANCOM/POU/Final_All.POU_level2_Class1_Pot_NonPot.xlsx")
write_xlsx(Final_All.POU_level6_Class1_Pot_NonPot,"C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Stats/ANCOM/POU/Final_All.POU_level6_Class1_Pot_NonPot.xlsx")
write_xlsx(Final_All.POU_OTU_Class1_Pot_NonPot,"C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Stats/ANCOM/POU/Final_All.POU_OTU_Class1_Pot_NonPot.xlsx")

write_xlsx(Final_All.POU_level2_Classification_1,"C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Stats/ANCOM/POU/Final_All.POU_level2_Classification_1.xlsx")
write_xlsx(Final_All.POU_level6_Classification_1,"C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Stats/ANCOM/POU/Final_All.POU_level6_Classification_1.xlsx")
write_xlsx(Final_All.POU_OTU_Classification_1,"C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Stats/ANCOM/POU/Final_All.POU_OTU_Classification_1.xlsx")

