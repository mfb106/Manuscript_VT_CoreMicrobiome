rm(list = ls())
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
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            #axis.text.x = element_text(angle = 45),
            axis.text = element_text(size=rel(1)), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.box="vertical", 
            legend.direction = "horizontal",
            legend.key.size= unit(0.1, "cm"),
            legend.text = element_text(size=rel(1)),
            legend.margin=margin(),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic",size=rel(1)),
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

#Read Meta Data 
Metadata_POCPOU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_DS_Special_V2', rowNames = TRUE)
#Metadata_POC<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'WTP_Final', rowNames = TRUE)
Metadata_POC<-Metadata_POCPOU[Metadata_POCPOU$Sample_Local == "POC",]
#Metadata_POU<-read.xlsx("MapingFile_Updated_2024.xlsx", sheet= 'DS_Final', rowNames = TRUE)
Metadata_POU<-Metadata_POCPOU[Metadata_POCPOU$Sample_Local == "POU",]

#Order Meta Data 
Metadata_POCPOU$Classification_1 = factor(Metadata_POCPOU$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POC$Classification_1 = factor(Metadata_POC$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))
Metadata_POU$Classification_1 = factor(Metadata_POU$Classification_1, levels=c("Potable Conventional","Potable Reuse","Non-potable Reuse","Blank"))

######All Samples, POCPOU
Meta.All.POCPOU<-Metadata_POCPOU
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'ConnerPrimer58',] #removed because of low read count
Meta.All.POCPOU<-Meta.All.POCPOU[Meta.All.POCPOU$Description != 'PR1WEEF',] #removed because of high read count 
Meta.All.POCPOU.Filter<-Meta.All.POCPOU[Meta.All.POCPOU$Classification_1 == 'Potable Conventional'| Meta.All.POCPOU$Classification_1 =='Non-potable Reuse'|Meta.All.POCPOU$Classification_1 =='Potable Reuse'|Meta.All.POCPOU$Classification_1 =='',]#can add Blank


#Alpha Diversity
df <- tibble::rownames_to_column(table_all_data, "OTU")
mydf1 <- data.frame(df %>%
                      pivot_longer(-`OTU`, names_to = "Sample", values_to = "Count"))

colnames(mydf1)<-c("OTU","Group","value")

mydf2<-mydf1[mydf1$value != 0, ]

rand <- mydf2 %>%
  uncount(value) %>%
  dplyr::mutate(name = sample(OTU)) %>%
  dplyr::count(Group, OTU, name="value")

Alpha_Diversity_All<-rand %>%
  dplyr::group_by(Group) %>%
  dplyr::summarize(richness = specnumber(value),
                   shannon = vegan::diversity(value, index="shannon"),
                   simpson = vegan::diversity(value, index="simpson"),
                   #invsimpson = 1/simpson,
                   n = sum(value))
Alpha_Diversity_All

Alpha_Diversity_All_longer<-pivot_longer(Alpha_Diversity_All, cols = c(richness, shannon, simpson),  #invsimpson,
                                        names_to = "Index")

ggplot(Alpha_Diversity_All_longer,aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Index, nrow=4, scales="free_y")



###Bringing in Metadata, Testing Potability
Combined_Meta_Data<-merge(Alpha_Diversity_All,Meta.All.POCPOU.Filter,by.y = 0, by.x = "Group")
Combined_Meta_Data_Longer<-merge(Alpha_Diversity_All_longer,Meta.All.POCPOU.Filter,by.y = 0, by.x = "Group")


#Plot, Class1_Pot_NonPot
All_potability<-ggplot(Combined_Meta_Data_Longer, aes(x=Class1_Pot_NonPot,y = value, fill = Class1_Pot_NonPot)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free_y")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
All_potability

#Plot, Classification_1 
All_Class1<-ggplot(Combined_Meta_Data_Longer, aes(x=Classification_1,y = value, fill = Classification_1)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free_y")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
All_Class1

#Plot, climate - Ignore
All_climate<-ggplot(Combined_Meta_Data, aes(x=Climate,y = richness, fill = Climate)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Classification_1,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
All_climate

#Plot, sample_local
All_POCPOU_1<-ggplot(Combined_Meta_Data, aes(x=Sample_Local,y = richness, fill = Sample_Local)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Classification_1,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
All_POCPOU_1

All_POCPOU_2<-ggplot(Combined_Meta_Data_Longer, aes(x=Sample_Local,y = value, fill = Sample_Local)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free_y")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
All_POCPOU_2

#Conventional Potable
Combined_Meta_Data_Conventional_Potable<-Combined_Meta_Data[Combined_Meta_Data$Classification_1 == "Potable Conventional", ]
Combined_Meta_Data_Conventional_Potable_Longer<-Combined_Meta_Data_Longer[Combined_Meta_Data_Longer$Classification_1 == "Potable Conventional", ]

#Sample Local
Conventional_Potable_POCPOU<-ggplot(Combined_Meta_Data_Conventional_Potable_Longer, aes(x=Sample_Local,y = value, fill = Sample_Local)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
Conventional_Potable_POCPOU

#Climate
Conventional_Potable_Climate<-ggplot(Combined_Meta_Data_Conventional_Potable_Longer, aes(x=Climate,y = value, fill = Climate)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
Conventional_Potable_Climate

#Region
Conventional_Potable_Region<-ggplot(Combined_Meta_Data_Conventional_Potable_Longer, aes(x=Region,y = value, fill = Region)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
Conventional_Potable_Region

#Final_Disinfection_Residual
Conventional_Potable_Disinfection<-ggplot(Combined_Meta_Data_Conventional_Potable_Longer, aes(x=Final_Disinfection_Residual,y = value, fill = Final_Disinfection_Residual)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
Conventional_Potable_Disinfection

#Pot._Source_Water
Conventional_Potable_Source<-ggplot(Combined_Meta_Data_Conventional_Potable_Longer, aes(x=Potable_Source_Water,y = value, fill = Potable_Source_Water)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
Conventional_Potable_Source

#Treatment_Classification_Number
Conventional_Potable_Treatment<-ggplot(Combined_Meta_Data_Conventional_Potable_Longer, aes(x=Treatment_Classification_Number,y = value, fill = Treatment_Classification_Number)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
Conventional_Potable_Treatment

#Treatment_Classification_Cat
Conventional_Potable_Treatment_Cat<-ggplot(Combined_Meta_Data_Conventional_Potable_Longer, aes(x=Treatment_Classification_Cat,y = value, fill = Treatment_Classification_Cat)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
Conventional_Potable_Treatment_Cat

##Potable Reuse 
Combined_Meta_Data_Potable_Reuse<-Combined_Meta_Data[Combined_Meta_Data$Classification_1 == "Potable Reuse", ]
Combined_Meta_Data_Potable_Reuse_Longer<-Combined_Meta_Data_Longer[Combined_Meta_Data_Longer$Classification_1 == "Potable Reuse", ]

#Sample Local
Potable_Reuse_POCPOU<-ggplot(Combined_Meta_Data_Potable_Reuse_Longer, aes(x=Sample_Local,y = value, fill = Sample_Local)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
Potable_Reuse_POCPOU

#Climate
Potable_Reuse_Climate<-ggplot(Combined_Meta_Data_Potable_Reuse_Longer, aes(x=Climate,y = value, fill = Climate)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
Potable_Reuse_Climate

#Region
Potable_Reuse_Region<-ggplot(Combined_Meta_Data_Potable_Reuse_Longer, aes(x=Region,y = value, fill = Region)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
Potable_Reuse_Region

#Final_Disinfection_Residual
Potable_Reuse_Disinfection<-ggplot(Combined_Meta_Data_Potable_Reuse_Longer, aes(x=Final_Disinfection_Residual,y = value, fill = Final_Disinfection_Residual)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
Potable_Reuse_Disinfection

#Pot._Source_Water
Potable_Reuse_Source<-ggplot(Combined_Meta_Data_Potable_Reuse_Longer, aes(x=Potable_Source_Water,y = value, fill = Potable_Source_Water)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
Potable_Reuse_Source

#Treatment_Classification_Number
Potable_Reuse_Treatment<-ggplot(Combined_Meta_Data_Potable_Reuse_Longer, aes(x=Treatment_Classification_Number,y = value, fill = Treatment_Classification_Number)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
Potable_Reuse_Treatment

#Treatment_Classification_Cat
Potable_Reuse_Treatment_Cat<-ggplot(Combined_Meta_Data_Potable_Reuse_Longer, aes(x=Treatment_Classification_Cat,y = value, fill = Treatment_Classification_Cat)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
Potable_Reuse_Treatment_Cat




##Non-Potable Reuse 
Combined_Meta_Data_NonPotable_Reuse<-Combined_Meta_Data[Combined_Meta_Data$Classification_1 == "Non-potable Reuse", ]
Combined_Meta_Data_NonPotable_Reuse_Longer<-Combined_Meta_Data_Longer[Combined_Meta_Data_Longer$Classification_1 == "Non-potable Reuse", ]

#Sample Local
NonPotable_Reuse_POCPOU<-ggplot(Combined_Meta_Data_NonPotable_Reuse_Longer, aes(x=Sample_Local,y = value, fill = Sample_Local)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
NonPotable_Reuse_POCPOU

#Climate
NonPotable_Reuse_Climate<-ggplot(Combined_Meta_Data_NonPotable_Reuse_Longer, aes(x=Climate,y = value, fill = Climate)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
NonPotable_Reuse_Climate

#Region
NonPotable_Reuse_Region<-ggplot(Combined_Meta_Data_NonPotable_Reuse_Longer, aes(x=Region,y = value, fill = Region)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
NonPotable_Reuse_Region

#Final_Disinfection_Residual
NonPotable_Reuse_Disinfection<-ggplot(Combined_Meta_Data_NonPotable_Reuse_Longer, aes(x=Final_Disinfection_Residual,y = value, fill = Final_Disinfection_Residual)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
NonPotable_Reuse_Disinfection

#Pot._Source_Water
NonPotable_Reuse_Source<-ggplot(Combined_Meta_Data_NonPotable_Reuse_Longer, aes(x=Potable_Source_Water,y = value, fill = Potable_Source_Water)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
NonPotable_Reuse_Source

#Treatment_Classification_Number
NonPotable_Reuse_Treatment<-ggplot(Combined_Meta_Data_NonPotable_Reuse_Longer, aes(x=Treatment_Classification_Number,y = value, fill = Treatment_Classification_Number)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
NonPotable_Reuse_Treatment

#Treatment_Classification_Cat
NonPotable_Reuse_Treatment_Cat<-ggplot(Combined_Meta_Data_NonPotable_Reuse_Longer, aes(x=Treatment_Classification_Cat,y = value, fill = Treatment_Classification_Cat)) +
  geom_violin()+
  geom_pwc(method = "wilcox.test", label = "p.signif", p.adjust.method="BY", hide.ns = FALSE, remove.bracket = FALSE )+
  geom_boxplot(width = 0.1, outlier.shape = 8, outlier.size = 2, outlier.colour = "red")+ #, show.legend = FALSE
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  facet_wrap(~Index,scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  scale_colour_Publication()+
  scale_x_discrete(labels = label_wrap(10))
NonPotable_Reuse_Treatment_Cat

##Export
ggsave("All_potability.jpeg",plot = All_potability,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("All_Class1.jpeg",plot = All_Class1,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("All_climate.jpeg",plot = All_climate,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("All_POCPOU_1.jpeg",plot = All_POCPOU_1,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("All_POCPOU_2.jpeg",plot = All_POCPOU_2,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)

ggsave("Conventional_Potable_POCPOU.jpeg",plot = Conventional_Potable_POCPOU,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("Conventional_Potable_Region.jpeg",plot = Conventional_Potable_Region,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("Conventional_Potable_Disinfection.jpeg",plot = Conventional_Potable_Disinfection,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("Conventional_Potable_Source.jpeg",plot = Conventional_Potable_Source,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("Conventional_Potable_Treatment.jpeg",plot = Conventional_Potable_Treatment,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("Conventional_Potable_Treatment_Cat.jpeg",plot = Conventional_Potable_Treatment_Cat,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)

ggsave("NonPotable_Reuse_POCPOU.jpeg",plot = NonPotable_Reuse_POCPOU,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("NonPotable_Reuse_Region.jpeg",plot = NonPotable_Reuse_Region,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("NonPotable_Reuse_Disinfection.jpeg",plot = NonPotable_Reuse_Disinfection,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("NonPotable_Reuse_Source.jpeg",plot = NonPotable_Reuse_Source,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("NonPotable_Reuse_Treatment.jpeg",plot = NonPotable_Reuse_Treatment,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("NonPotable_Reuse_Treatment_Cat.jpeg",plot = NonPotable_Reuse_Treatment_Cat,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)

ggsave("Potable_Reuse_POCPOU.jpeg",plot = Potable_Reuse_POCPOU,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("Potable_Reuse_Region.jpeg",plot = Potable_Reuse_Region,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("Potable_Reuse_Disinfection.jpeg",plot = Potable_Reuse_Disinfection,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("Potable_Reuse_Source.jpeg",plot = Potable_Reuse_Source,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("Potable_Reuse_Treatment.jpeg",plot = Potable_Reuse_Treatment,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("Potable_Reuse_Treatment_Cat.jpeg",plot = Potable_Reuse_Treatment_Cat,path="C:/Users/matth/OneDrive - Virginia Tech/Manuscripts/Amplicon_Comparison/Figures/R_Visualization_Rerun/Silva/Alpha",scale = 1, width = 12, height = 8, units = "in", dpi = 300)

