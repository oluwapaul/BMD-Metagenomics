#Script to get TIGRFAMs annotations and larger roles in one file

#Load library
library(tidyverse)

setwd("../OutputFiles/")

#load all files
#see comments in files for directions on formatting
tigrfam <- read.table("../TIGRfam_Annotations_Roles/TIGRfam_annotations.txt", sep = "\t", header = T)
tigr_link <- read.table("../TIGRfam_Annotations_Roles/TIGRfam_ROLE_LINK.txt", sep = "\t", header = T)
mainrole <- read.table("../TIGRfam_Annotations_Roles/TIGRfam_mainroles.txt", sep = "\t", header = T)
subrole <- read.table("../TIGRfam_Annotations_Roles/TIGRfam_subroles.txt", sep = "\t", header = T)

tigr_main <- merge(tigr_link, mainrole, by = "link", all.x=T, all.y =T)
tigr_sub <- merge(tigr_main, subrole, by = "link", all.x=T, all.y =T)
tigr_func <- merge(tigr_sub, tigrfam, by = "TIGRFAM", all.x=T, all.y =T)

write.table(tigr_func, file = "TIGR_func_full.tsv", sep = "\t", row.names = F, quote = F)
