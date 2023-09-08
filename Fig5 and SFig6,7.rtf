#### Conjugation and Phage differentially abundant genes
#Stating home directory
setwd(dir = "~/Desktop/BMD/DESeq/Conjugation_phage")

my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)

#Loading data
d7_Subtherapeutic <- read.delim("d7_Subtherapeutic.tsv", sep = "\t", header = T,check.names=FALSE, row.names = 1)
d7_Therapeutic <- read.delim("d7_Therapeutic.tsv", sep = "\t", header = T,check.names=FALSE, row.names = 1)
d35_Subtherapeutic <- read.delim("d35_Subtherapeutic.tsv", sep = "\t", header = T,check.names=FALSE, row.names = 1)
d35_Therapeutic <- read.delim("d35_Therapeutic.tsv", sep = "\t", header = T,check.names=FALSE, row.names = 1)
d78_Subtherapeutic <- read.delim("d78_Subtherapeutic.tsv", sep = "\t", header = T,check.names=FALSE, row.names = 1)
d78_Therapeutic <- read.delim("d78_Therapeutic.tsv", sep = "\t", header = T,check.names=FALSE, row.names = 1)

#Sub-setting conjugation and phage
d7_Subtherapeutic_conj <- subset(d7_Subtherapeutic, MGE=="conjugation")
d7_Subtherapeutic_phage <- subset(d7_Subtherapeutic, MGE=="phage")
d7_Subtherapeutic_transposon <- subset(d7_Subtherapeutic, MGE=="transposon")

d7_Therapeutic_conj <- subset(d7_Therapeutic, MGE=="conjugation")
d7_Therapeutic_phage <- subset(d7_Therapeutic, MGE=="phage")
d7_Therapeutic_transposon <- subset(d7_Therapeutic, MGE=="transposon")

d35_Subtherapeutic_conj <- subset(d35_Subtherapeutic, MGE=="conjugation")
d35_Subtherapeutic_phage <- subset(d35_Subtherapeutic, MGE=="phage")
d35_Subtherapeutic_transposon <- subset(d35_Subtherapeutic, MGE=="transposon")


d35_Therapeutic_conj <- subset(d35_Therapeutic, MGE=="conjugation")
d35_Therapeutic_phage <- subset(d35_Therapeutic, MGE=="phage")
d35_Therapeutic_transposon <- subset(d35_Therapeutic, MGE=="transposon")

d78_Subtherapeutic_conj <- subset(d78_Subtherapeutic, MGE=="conjugation")
d78_Subtherapeutic_phage <- subset(d78_Subtherapeutic, MGE=="phage")
d78_Subtherapeutic_transposon <- subset(d78_Subtherapeutic, MGE=="transposon")

d78_Therapeutic_conj <- subset(d78_Therapeutic, MGE=="conjugation")
d78_Therapeutic_phage <- subset(d78_Therapeutic, MGE=="phage")
d78_Therapeutic_transposon <- subset(d78_Therapeutic, MGE=="transposon")

#sigtab_pos <- subset(sigtab, log2FoldChange>0)
#sigtab_pos$subrole <- factor(sigtab_pos$subrole, levels = unique(sigtab[order(sigtab$mainrole), "subrole"]))

#Plotting the bar graph
#Day 7 Subtherapeutic

#mutate(A = fct_reorder(as.character(B), A))

d7_Subtherapeutic_conj$subsystem <- factor(d7_Subtherapeutic_conj$subsystem, levels = c("Others", "Mating Pair Stabilization", "Pilus Assembly/Extension", "Pilus Retraction", "ProPilin Maturation", "Relaxase", "Relaxosome Accessory", "T4SS Core Proteins", "Transcription Regulation"))
d7_Subtherapeutic_conj <- d7_Subtherapeutic_conj[order(d7_Subtherapeutic_conj$subsystem, d7_Subtherapeutic_conj$Name),] 
d7_Subtherapeutic_conj %>%
  mutate(Name = fct_reorder(as.character(subsystem), Name)) %>%
 ggplot( aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Type IV Conjugation System", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()

d7_subther_conj <- d7_subther_conj + coord_flip()
ggsave("d7_subthera_conj.png", height 
       = 6, width = 7)


d7_Subtherapeutic_conj$subsystem <- factor(d7_Subtherapeutic_conj$subsystem, levels = c("Others", "Relaxosome Accessory", "Relaxase", "Mating Pair Stabilization", "Pilus Retraction", "Pilus Assembly/Extension", "T4SS Core Proteins", "ProPilin Maturation", "Transcription Regulation"))
d7_Subtherapeutic_conj$Name <- factor(d7_Subtherapeutic_conj$Name, levels = c("VirB11", "TraD_Ftype", "SXT_TraD", "DotA_TraY", "relax_trwC", "TraN_Ftype", "TrbI_Ftype", "TrbC_Ftype", "TrbB", "TraW", "TraL_TIGR", "TraF", "TraE_TIGR", "TraC-F-type", "VirB9", "TraV", "TraK_Ftype", "TraX_Ftype", "TraQ", "TraA_TIGR", "plasmid_TraJ"))

#mutate(d7_Subtherapeutic_conj,Name = fct_reorder(as.character(subsystem), Name))

d7_subther_conj <-ggplot(d7_Subtherapeutic_conj, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Type IV Conjugation System", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
  
d7_subther_conj <- d7_subther_conj + coord_flip()
ggsave("d7_subthera_conj.png", height 
       = 6, width = 7)


d7_Subtherapeutic_phage$subsystem <- factor(d7_Subtherapeutic_phage$subsystem, levels = c("Others", "Bacterial Defence", "Bacterial Lysis", "Phage Tail", "Phage Head", "Phage DNA Packaging","Prophage Replication & Expression"))
d7_Subtherapeutic_phage$Name <- factor(d7_Subtherapeutic_phage$Name, levels = c("phg_TIGR02220", "phge_rel_HI1409", "C4_traR_proteo", "RNA_lig_RNL2",  "RNA_lig_T4_1", "tail_comp_S", "LGT_TIGR03299", "phageshock_pspE", "phageshock_pspA", "phageshock_pspB", "phageshock_pspC", "phageshock_pspG", "phageshock_pspD","holin_lambda",
                                                                                "phage_LysB", "maj_tail_phi13", "tape_meas_lam_C", "phi3626_gp14_N", "phage_tail_L", "phage_lam_T", "tail_TIGR02242", "tail_P2_I", "major_cap_HK97","major_capsid_P2", "proheadase_HK97", "portal_HK97", "portal_lambda", "A118_put_portal",
                                                                                "put_anti_recept", "phage_pRha", "phage_rep_org_N", "phage_term_2", "phage_arpU", "phage_rinA", "rep_II_X", "phage_P_loop","bet_lambda"))

d7_subther_phage <- ggplot(d7_Subtherapeutic_phage, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Phage Subsystem", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE))+
  theme_classic()
d7_subther_phage <- d7_subther_phage + coord_flip()
ggsave("d7_subthera_phage.png", height 
       = 6, width = 7)

d7_subther_trans <- ggplot(d7_Subtherapeutic_transposon, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Conjugation transposon", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d7_subther_trans <- d7_subther_trans + coord_flip()
ggsave("d7_subthera_trans2.png", height 
       = 1.5, width = 5)


#Day 7 Therapeutic
d7_Therapeutic_conj$subsystem <- factor(d7_Therapeutic_conj$subsystem, levels = c("Others", "Relaxosome Accessory", "Relaxase", "Mating Pair Stabilization", "Pilus Retraction", "Pilus Assembly/Extension", "T4SS Core Proteins", "ProPilin Maturation", "Transcription Regulation"))
d7_Therapeutic_conj$Name <- factor(d7_Therapeutic_conj$Name, levels = c("VirB11", "TraD_Ftype", "SXT_TraD", "DotA_TraY", "relax_trwC", "TraN_Ftype", "TrbI_Ftype", "TrbC_Ftype", "TrbB", "TraW", "TraL_TIGR", "TraF", "TraE_TIGR", "TraC-F-type", "VirB9", "TraV", "TraK_Ftype", "VirB4_CagE", "TraX_Ftype", "TraQ", "TraA_TIGR", "plasmid_TraJ"))
d7_ther_conj <- ggplot(d7_Therapeutic_conj, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Type IV Conjugation System", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d7_ther_conj <- d7_ther_conj + coord_flip()
ggsave("d7_thera_conj.png", height 
       = 6, width = 7)


d7_Therapeutic_phage$subsystem <- factor(d7_Therapeutic_phage$subsystem, levels = c("Others", "Bacterial Defence", "Bacterial Lysis", "Phage Tail", "Phage Head", "Phage DNA Packaging","Prophage Replication & Expression"))
d7_Therapeutic_phage$Name <- factor(d7_Therapeutic_phage$Name, levels = c("phg_TIGR02220", "phge_rel_HI1409", "C4_traR_proteo", "RNA_lig_RNL2",  "RNA_lig_T4_1", "tail_comp_S", "LGT_TIGR03299", "phageshock_pspE", "phageshock_pspA", "phageshock_pspB", "phageshock_pspC", "phageshock_pspG", "phageshock_pspD","holin_lambda",
                                                                                "phage_LysB", "maj_tail_phi13", "tape_meas_lam_C", "phi3626_gp14_N", "phage_tail_L", "phage_lam_T", "tail_TIGR02242", "tail_P2_I", "major_cap_HK97","major_capsid_P2", "proheadase_HK97", "portal_HK97", "portal_lambda", "A118_put_portal",
                                                                                "put_anti_recept", "phage_pRha", "phage_rep_org_N", "phage_term_2", "phage_arpU", "phage_rinA", "rep_II_X", "phage_P_loop", "sm_term_P27", "bet_lambda"))

d7_ther_phage <- ggplot(d7_Therapeutic_phage, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Phage Subsystem", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d7_ther_phage <- d7_ther_phage + coord_flip()
ggsave("d7_thera_phage.png", height 
       = 6, width = 7)

d7_ther_trans <- ggplot(d7_Therapeutic_transposon, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Conjugation transposon", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d7_ther_trans <- d7_ther_trans + coord_flip()
ggsave("d7_thera_trans2.png", height 
       = 1.5, width = 5)


#Day 35 Subtherapeutic
d35_Subtherapeutic_conj$subsystem <- factor(d35_Subtherapeutic_conj$subsystem, levels = c("Others", "Relaxosome Accessory", "Relaxase", "Mating Pair Stabilization", "Pilus Retraction", "Pilus Assembly/Extension", "T4SS Core Proteins", "ProPilin Maturation", "Transcription Regulation"))
d35_Subtherapeutic_conj$Name <- factor(d35_Subtherapeutic_conj$Name, levels = c("VirB11", "TraD_Ftype", "SXT_TraD", "DotA_TraY", "relax_trwC", "TraN_Ftype", "TrbI_Ftype", "TrbC_Ftype", "TrbB", "TraW", "TraL_TIGR", "TraF", "TraE_TIGR", "TraF_Ti", "TraC-F-type", "VirB9", "TraV", "TraK_Ftype", "TraX_Ftype", "TraQ", "TraA_TIGR", "plasmid_TraJ"))
d35_subther_conj <- ggplot(d35_Subtherapeutic_conj, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge", width = 1) +
  scale_fill_manual("Type IV Conjugation System", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d35_subther_conj <- d35_subther_conj + coord_flip()
ggsave("d35_subthera_conj.png", height 
       = 6, width = 7)


d35_Subtherapeutic_phage$subsystem <- factor(d35_Subtherapeutic_phage$subsystem, levels = c("Others", "Bacterial Defence", "Bacterial Lysis", "Phage Tail", "Phage Head", "Phage DNA Packaging","Prophage Replication & Expression"))
d35_Subtherapeutic_phage$Name <- factor(d35_Subtherapeutic_phage$Name, levels = c("phg_TIGR02220", "phge_rel_HI1409", "C4_traR_proteo", "RNA_lig_RNL2",  "phage_chp_gp8", "psiM2_ORF9", "tail_comp_S", "LGT_TIGR03299", "phageshock_pspE", "phageshock_pspA", "phageshock_pspB", "phageshock_pspC", "phageshock_pspG", "phageshock_pspD","holin_SPP1",
                                                                                "holin_phiLC3", "phgtail_TP901_1", "tape_meas_lam_C", "phi3626_gp14_N", "phage_tail_L", "phage_lam_T", "tail_TIGR02242", "tail_P2_I", "major_cap_HK97","phageSPP1_gp7", "proheadase_HK97", "portal_HK97", "portal_lambda", "portal_SPP1", "put_DNA_pack",
                                                                                "put_anti_recept", "phage_pRha", "phage_rep_org_N", "phage_term_2", "phage_arpU", "phage_rinA", "rep_II_X", "phage_P_loop","bet_lambda"))

d35_subther_phage <- ggplot(d35_Subtherapeutic_phage, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Phage Subsystem", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d35_subther_phage <- d35_subther_phage + coord_flip()
ggsave("d35_subthera_phage.png", height 
       = 6, width = 7)

d35_subther_trans <- ggplot(d35_Subtherapeutic_transposon, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Conjugation transposon", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d35_subther_trans <- d35_subther_trans + coord_flip()
ggsave("d35_subthera_trans2.png", height 
       = 1.5, width = 5)

#Day 35 Therapeutic
d35_Therapeutic_conj$subsystem <- factor(d35_Therapeutic_conj$subsystem, levels = c("Others", "Relaxosome Accessory", "Relaxase", "Mating Pair Stabilization", "Pilus Retraction", "Pilus Assembly/Extension", "T4SS Core Proteins", "ProPilin Maturation", "Transcription Regulation"))
d35_Therapeutic_conj$Name <- factor(d35_Therapeutic_conj$Name, levels = c("VirB11", "TraD_Ftype", "SXT_TraD", "DotA_TraY", "relax_trwC", "TraN_Ftype", "TrbI_Ftype", "TrbC_Ftype", "TrbB", "TraW","TraE_TIGR", "TraL_TIGR", "TraF-like", "TraB", "VirB9", "TraV", "TraK_Ftype", "TraX_Ftype", "TraQ", "TraA_TIGR", "plasmid_TraJ"))
d35_ther_conj <- ggplot(d35_Therapeutic_conj, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Type IV Conjugation System", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d35_ther_conj <- d35_ther_conj + coord_flip()
ggsave("d35_thera_conj.png", height 
       = 6, width = 7)

d35_Therapeutic_phage$subsystem <- factor(d35_Therapeutic_phage$subsystem, levels = c("Others", "Bacterial Defence", "Bacterial Lysis", "Phage Tail", "Phage Head", "Phage DNA Packaging","Prophage Replication & Expression"))
d35_Therapeutic_phage$Name <- factor(d35_Therapeutic_phage$Name, levels = c("phg_TIGR02220", "phge_rel_HI1409", "C4_traR_proteo", "RNA_lig_RNL2",  "phage_chp_gp8", "psiM2_ORF9", "tail_comp_S", "LGT_TIGR03299", "phageshock_pspE", "phageshock_pspA", "phageshock_pspB", "phageshock_pspC", "phageshock_pspG", "phageshock_pspD","holin_SPP1",
                                                                                  "holin_phiLC3", "maj_tail_phi13", "phgtail_TP901_1", "phage_lambda_G", "tape_meas_lam_C", "phi3626_gp14_N", "phage_tail_L", "phage_lam_T", "tail_TIGR02242", "tail_P2_I", "holin_LLH", "major_cap_HK97","phageSPP1_gp7", "gp16_SPP1", "major_capsid_P2", "proheadase_HK97", "portal_lambda", "A118_put_portal", "portal_SPP1", "portal_HK97", "put_DNA_pack",
                                                                                  "put_anti_recept", "phage_pRha", "phage_rep_org_N", "phage_term_2", "phage_arpU", "phage_rinA", "phage_O_Nterm", "rep_II_X", "phage_P_loop","bet_lambda"))
d35_ther_phage <- ggplot(d35_Therapeutic_phage, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Phage Subsystem", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d35_ther_phage <- d35_ther_phage + coord_flip()
ggsave("d35_thera_phage.png", height 
       = 6, width = 7)

d35_ther_trans <- ggplot(d35_Therapeutic_transposon, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Conjugation transposon", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d35_ther_trans <- d35_ther_trans + coord_flip()
ggsave("d35_thera_trans2.png", height 
       = 1.5, width = 5)


#Day78 Subtherapeutic
d78_Subtherapeutic_conj$subsystem <- factor(d78_Subtherapeutic_conj$subsystem, levels = c("Others", "Relaxosome Accessory", "Relaxase", "Mating Pair Stabilization", "Pilus Retraction", "Pilus Assembly/Extension", "T4SS Core Proteins", "ProPilin Maturation", "Transcription Regulation"))
d78_Subtherapeutic_conj$Name <- factor(d78_Subtherapeutic_conj$Name, levels = c("VirB11", "TraD_Ftype", "SXT_TraD", "DotA_TraY", "relax_trwC", "TraN_Ftype", "TrbI_Ftype", "TrbC_Ftype", "TrbB", "TraW", "TraL_TIGR", "TraF", "TraE_TIGR", "TraC-F-type", "VirB9", "TraV", "TraK_Ftype", "TraX_Ftype", "TraQ", "TraA_TIGR", "plasmid_TraJ"))
d78_subther_conj <- ggplot(d78_Subtherapeutic_conj, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Type IV Conjugation System", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d78_subther_conj <- d78_subther_conj + coord_flip()
ggsave("d78_subthera_conj.png", height 
       = 6, width = 7)

d78_Subtherapeutic_phage$subsystem <- factor(d78_Subtherapeutic_phage$subsystem, levels = c("Others", "Bacterial Defence", "Bacterial Lysis", "Phage Tail", "Phage Head", "Phage DNA Packaging","Prophage Replication & Expression"))
d78_Subtherapeutic_phage$Name <- factor(d78_Subtherapeutic_phage$Name, levels = c("phg_TIGR02220", "primase_Cterm", "phge_rel_HI1409", "C4_traR_proteo", "RNA_lig_RNL2",  "phage_chp_gp8", "psiM2_ORF9", "tail_comp_S", "LGT_TIGR03299", "phageshock_pspE", "phageshock_pspA", "phageshock_pspB", "phageshock_pspC", "phageshock_pspG", "phageshock_pspD","holin_lambda",
                                                                            "phage_LysB", "maj_tail_phi13", "phgtail_TP901_1", "phage_lambda_G", "tape_meas_lam_C", "phi3626_gp14_N", "phage_tail_L", "phage_lam_T", "tail_TIGR02242", "tail_P2_I", "holin_LLH", "major_cap_HK97","phageSPP1_gp7", "gp16_SPP1", "major_capsid_P2", "proheadase_HK97", "portal_lambda", "A118_put_portal", "portal_SPP1", "portal_HK97", "put_DNA_pack",
                                                                            "put_anti_recept", "phage_pRha", "phage_rep_org_N", "phage_term_2", "phage_arpU", "phage_rinA", "phage_O_Nterm", "rep_II_X", "phage_P2_V", "phage_P_loop","bet_lambda"))

d78_subther_phage <- ggplot(d78_Subtherapeutic_phage, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Phage Subsystem", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d78_subther_phage <- d78_subther_phage + coord_flip()
ggsave("d78_subthera_phage.png", height 
       = 6, width = 7)

d78_subther_trans <- ggplot(d78_Subtherapeutic_transposon, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Conjugation transposon", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d78_subther_trans <- d78_subther_trans + coord_flip()
ggsave("d78_subthera_trans2.png", height 
       = 1.5, width = 5)


#Day 78 Therapeutic
d78_Therapeutic_conj$subsystem <- factor(d78_Subtherapeutic_conj$subsystem, levels = c("Others", "Relaxosome Accessory", "Relaxase", "Mating Pair Stabilization", "Pilus Retraction", "Pilus Assembly/Extension", "T4SS Core Proteins", "ProPilin Maturation", "Transcription Regulation"))
d78_Therapeutic_conj$Name <- factor(d78_Therapeutic_conj$Name, levels = c("VirB11", "TraD_Ftype", "SXT_TraD", "DotA_TraY", "relax_trwC", "TraN_Ftype", "TrbI_Ftype", "TrbC_Ftype", "TraF-like", "TraW", "TraL_TIGR", "TraF", "TraE_TIGR", "TraC-F-type", "VirB9", "TraV", "TraK_Ftype", "TraX_Ftype", "TraQ", "TraA_TIGR", "plasmid_TraJ"))
d78_ther_conj <- ggplot(d78_Therapeutic_conj, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Type IV Conjugation System", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d78_ther_conj <- d78_ther_conj + coord_flip()
ggsave("d78_thera_conj.png", height = 6, width = 7)

d78_Therapeutic_phage$subsystem <- factor(d78_Therapeutic_phage$subsystem, levels = c("Others", "Bacterial Defence", "Bacterial Lysis", "Phage Tail", "Phage Head", "Phage DNA Packaging","Prophage Replication & Expression"))
d78_Therapeutic_phage$Name <- factor(d78_Therapeutic_phage$Name, levels = c("phg_TIGR02220", "primase_Cterm", "phge_rel_HI1409", "C4_traR_proteo", "RNA_lig_RNL2", "RNA_lig_T4_1",  "phage_chp_gp8", "psiM2_ORF9", "tail_comp_S", "LGT_TIGR03299", "phageshock_pspE", "phageshock_pspA", "phageshock_pspB", "phageshock_pspC", "phageshock_pspG", "phageshock_pspD","holin_SPP1","holin_lambda", "phage_LysB",
                                                                                  "holin_phiLC3", "maj_tail_phi13", "phgtail_TP901_1", "phage_lambda_G", "tape_meas_lam_C", "phi3626_gp14_N", "phage_tail_L", "phage_lam_T", "tail_TIGR02242", "tail_P2_I", "holin_LLH", "major_cap_HK97","phageSPP1_gp7", "gp16_SPP1", "major_capsid_P2", "proheadase_HK97", "portal_lambda", "A118_put_portal", "portal_SPP1", "portal_HK97", "put_DNA_pack",
                                                                                  "put_anti_recept", "phage_pRha", "phage_rep_org_N", "phage_term_2", "phage_arpU", "phage_rinA", "phage_O_Nterm", "rep_II_X", "phage_P2_V", "phage_P_loop","bet_lambda"))
d78_ther_phage <- ggplot(d78_Therapeutic_phage, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge") +
  scale_fill_manual("Phage Subsystem", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d78_ther_phage <- d78_ther_phage + coord_flip()
ggsave("d78_thera_phage.png", height 
       = 6, width = 7)

d78_ther_trans <- ggplot(d78_Therapeutic_transposon, aes(x=Name, y=log2FoldChange)) +
  geom_col(aes(fill = subsystem), position = "dodge", width = 1) +
  scale_fill_manual("Conjugation transposon", values = my_colors) +
  xlab("") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_classic()
d78_ther_trans <- d78_ther_trans + coord_flip()
ggsave("d78_thera_trans2.png", height 
       = 1.5, width = 5)
