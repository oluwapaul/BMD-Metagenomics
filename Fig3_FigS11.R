setwd(dir = "~/R_script")

Function <- read.table("indole_function_rpm_unannotated_normalized.tsv", header=TRUE, sep="\t",check.names=FALSE, row.names = 1)
meta <- read.csv("mg_metadata.csv")
meta <- meta[-which(meta$sample %in% c("d07sub045")), ] 

Function2 <- merge(meta, Function, by.x = "sample", by.y = 0)
my_colors <- c("red", "green", "blue")


#Function2 <- subset(Function2, day =="7")

Function2$antibiotic <- factor(Function2$antibiotic, levels = c("Control", "Subtherapeutic", "Therapeutic"))
Function2$antibiotic <- as.factor(Function2$antibiotic)
Function2$day <- as.factor(Function2$day)

Function2_d7 <- subset(Function2, day =="7")
Function2_d35 <- subset(Function2, day =="35")
Function2_d78 <- subset(Function2, day =="78")

#Difference between treatment group and donor at d0
#Compare each treatment mean to the donor mean
#options(contrasts=c("contr.sum", "contr.poly"))
options(contrasts=c("contr.treatment", "contr.poly"))
#Tryptophan synthase alpha
#d7
m1_tsynthase <- lm(tryptophan_synthase_alpha_subunit_ ~ antibiotic, data = Function2_d7)
summary(m1_tsynthase)
m1_tsynthaseb <- aov(tryptophan_synthase_alpha_subunit_ ~ antibiotic, data = Function2_d7)
summary(m1_tsynthaseb)
TukeyHSD(m1_tsynthaseb)
#d35
m2_tsynthase <- lm(tryptophan_synthase_alpha_subunit_ ~ antibiotic, data = Function2_d35)
summary(m2_tsynthase)
m2_tsynthaseb <- aov(tryptophan_synthase_alpha_subunit_ ~ antibiotic, data = Function2_d35)
summary(m2_tsynthaseb)
TukeyHSD(m2_tsynthaseb)
#d35
m3_tsynthase <- lm(tryptophan_synthase_alpha_subunit_ ~ antibiotic, data = Function2_d78)
summary(m3_tsynthase)
m3_tsynthaseb <- aov(tryptophan_synthase_alpha_subunit_ ~ antibiotic, data = Function2_d78)
summary(m3_tsynthaseb)
TukeyHSD(m3_tsynthaseb)

#Tryptophan synthase beta
#d7
m1_tsynthasebet <- lm(tryptophan_synthase_beta_subunit_ ~ antibiotic, data = Function2_d7)
summary(m1_tsynthasebet)
m1_tsynthasebetb <- aov(tryptophan_synthase_beta_subunit_ ~ antibiotic, data = Function2_d7)
summary(m1_tsynthasebetb)
TukeyHSD(m1_tsynthasebetb)
#d35
m2_tsynthasebet <- lm(tryptophan_synthase_beta_subunit_ ~ antibiotic, data = Function2_d35)
summary(m2_tsynthasebet)
m2_tsynthasebetb <- aov(tryptophan_synthase_beta_subunit_ ~ antibiotic, data = Function2_d35)
summary(m2_tsynthasebetb)
TukeyHSD(m2_tsynthasebetb)
#d35
m3_tsynthasebet <- lm(tryptophan_synthase_beta_subunit_ ~ antibiotic, data = Function2_d78)
summary(m3_tsynthasebet)
m3_tsynthasebetb <- aov(tryptophan_synthase_beta_subunit_ ~ antibiotic, data = Function2_d78)
summary(m3_tsynthasebetb)
TukeyHSD(m3_tsynthasebetb)

#aminodeoxychorismate_synthase
#d7
m1_chorsynthase <- lm(glutamine_amidotransferase_anthranilate_synthase_aminodeoxychorismate_synthase ~ antibiotic, data = Function2_d7)
summary(m1_chorsynthase)
m1_chorsynthaseb <- aov(glutamine_amidotransferase_anthranilate_synthase_aminodeoxychorismate_synthase ~ antibiotic, data = Function2_d7)
summary(m1_chorsynthaseb)
TukeyHSD(m1_chorsynthaseb)
#d35
m2_chorsynthase <- lm(glutamine_amidotransferase_anthranilate_synthase_aminodeoxychorismate_synthase ~ antibiotic, data = Function2_d35)
summary(m2_chorsynthase)
m2_chorsynthaseb <- aov(glutamine_amidotransferase_anthranilate_synthase_aminodeoxychorismate_synthase ~ antibiotic, data = Function2_d35)
summary(m2_chorsynthaseb)
TukeyHSD(m2_chorsynthaseb)
#d78
m3_chorsynthase <- lm(glutamine_amidotransferase_anthranilate_synthase_aminodeoxychorismate_synthase ~ antibiotic, data = Function2_d78)
summary(m3_chorsynthase)
m3_chorsynthaseb <- aov(glutamine_amidotransferase_anthranilate_synthase_aminodeoxychorismate_synthase ~ antibiotic, data = Function2_d78)
summary(m3_chorsynthaseb)
TukeyHSD(m3_chorsynthaseb)

#anthranilase synthase component 1
#d7
m1_ansynthase <- lm(anthranilate_synthase_component_I~ antibiotic, data = Function2_d7)
summary(m1_ansynthase)
m1_ansynthaseb <- aov(anthranilate_synthase_component_I~ antibiotic, data = Function2_d7)
summary(m1_ansynthaseb)
TukeyHSD(m1_ansynthaseb)
#d35
m2_ansynthase <- lm(anthranilate_synthase_component_I~ antibiotic, data = Function2_d35)
summary(m2_ansynthase)
m2_ansynthaseb <- aov(anthranilate_synthase_component_I~ antibiotic, data = Function2_d35)
summary(m2_ansynthaseb)
TukeyHSD(m2_ansynthaseb)
#d78
m3_ansynthase <- lm(anthranilate_synthase_component_I~ antibiotic, data = Function2_d78)
summary(m3_ansynthase)
m3_ansynthaseb <- aov(anthranilate_synthase_component_I~ antibiotic, data = Function2_d78)
summary(m3_ansynthaseb)
TukeyHSD(m3_ansynthaseb)

#anthranilate phosphoribosyltransferase
#d7
m1_anribtransferase <- lm(anthranilate_phosphoribosyltransferase~ antibiotic, data = Function2_d7)
summary(m1_anribtransferase)
m1_anribtransferaseb <- aov( anthranilate_phosphoribosyltransferase~ antibiotic, data = Function2_d7)
summary(m1_anribtransferaseb)
TukeyHSD(m1_anribtransferaseb)
#d35
m2_anribtransferase <- lm(anthranilate_phosphoribosyltransferase~ antibiotic, data = Function2_d35)
summary(m2_anribtransferase)
m2_anribtransferaseb <- aov( anthranilate_phosphoribosyltransferase~ antibiotic, data = Function2_d35)
summary(m2_anribtransferaseb)
TukeyHSD(m2_anribtransferaseb)
#d78
m3_anribtransferase <- lm(anthranilate_phosphoribosyltransferase~ antibiotic, data = Function2_d78)
summary(m3_anribtransferase)
m3_anribtransferaseb <- aov(anthranilate_phosphoribosyltransferase~ antibiotic, data = Function2_d78)
summary(m3_anribtransferaseb)
TukeyHSD(m3_anribtransferaseb)

#chorismate_mutase
#d7
m1_chormutase <- lm(chorismate_mutase~ antibiotic, data = Function2_d7)
summary(m1_chormutase)
m1_chormutaseb <- aov(chorismate_mutase~ antibiotic, data = Function2_d7)
summary(m1_chormutaseb)
TukeyHSD(m1_chormutaseb)
#d35
m2_chormutase <- lm(chorismate_mutase~ antibiotic, data = Function2_d35)
summary(m2_chormutase)
m2_chormutaseb <- aov(chorismate_mutase~ antibiotic, data = Function2_d35)
summary(m2_chormutaseb)
TukeyHSD(m2_chormutaseb)
#d78
m3_chormutase <- lm(chorismate_mutase~ antibiotic, data = Function2_d78)
summary(m3_chormutase)
m3_chormutaseb <- aov(chorismate_mutase~ antibiotic, data = Function2_d78)
summary(m3_chormutaseb)
TukeyHSD(m3_chormutaseb)

#chorismate_synthase
#d7
m1_chorissynthase <- lm(chorismate_synthase~ antibiotic, data = Function2_d7)
summary(m1_chorissynthase)
m1_chorissynthaseb <- aov(chorismate_synthase~ antibiotic, data = Function2_d7)
summary(m1_chorissynthaseb)
TukeyHSD(m1_chorissynthaseb)
#d35
m2_chorissynthase <- lm(chorismate_synthase~ antibiotic, data = Function2_d35)
summary(m2_chorissynthase)
m2_chorissynthaseb <- aov(chorismate_synthase~ antibiotic, data = Function2_d35)
summary(m2_chorissynthaseb)
TukeyHSD(m2_chorissynthaseb)
#78
m3_chorissynthase <- lm(chorismate_synthase~ antibiotic, data = Function2_d78)
summary(m3_chorissynthase)
m3_chorissynthaseb <- aov(chorismate_synthase~ antibiotic, data = Function2_d78)
summary(m3_chorissynthaseb)
TukeyHSD(m3_chorissynthaseb)

#isochorismate_synthase
#d7
m1_isochorsynthase <- lm( isochorismate_synthase~ antibiotic, data = Function2_d7)
summary(m1_isochorsynthase)
m1_isochorsynthaseb <- aov( isochorismate_synthase~ antibiotic, data = Function2_d7)
summary(m1_isochorsynthaseb)
TukeyHSD(m1_isochorsynthaseb)
#d35
m2_isochorsynthase <- lm( isochorismate_synthase~ antibiotic, data = Function2_d35)
summary(m2_isochorsynthase)
m2_isochorsynthaseb <- aov( isochorismate_synthase~ antibiotic, data = Function2_d35)
summary(m2_isochorsynthaseb)
TukeyHSD(m2_isochorsynthaseb)
#d78
m3_isochorsynthase <- lm( isochorismate_synthase~ antibiotic, data = Function2_d78)
summary(m3_isochorsynthase)
m3_isochorsynthaseb <- aov( isochorismate_synthase~ antibiotic, data = Function2_d78)
summary(m3_isochorsynthaseb)
TukeyHSD(m3_isochorsynthaseb)

#aminodeoxychorismate_lyase
#d7
m1_aminodeoxychorismate_lyase <- lm(aminodeoxychorismate_lyase ~ antibiotic, data = Function2_d7)
summary(m1_aminodeoxychorismate_lyase)
m1_aminodeoxychorismate_lyaseb <- aov(aminodeoxychorismate_lyase ~ antibiotic, data = Function2_d7)
summary(m1_aminodeoxychorismate_lyaseb)
TukeyHSD(m1_aminodeoxychorismate_lyaseb)
#d35
m2_aminodeoxychorismate_lyase <- lm(aminodeoxychorismate_lyase ~ antibiotic, data = Function2_d35)
summary(m2_aminodeoxychorismate_lyase)
m2_aminodeoxychorismate_lyaseb <- aov(aminodeoxychorismate_lyase ~ antibiotic, data = Function2_d35)
summary(m2_aminodeoxychorismate_lyaseb)
TukeyHSD(m2_aminodeoxychorismate_lyaseb)
#d78
m3_aminodeoxychorismate_lyase <- lm(aminodeoxychorismate_lyase ~ antibiotic, data = Function2_d78)
summary(m3_aminodeoxychorismate_lyase)
m3_aminodeoxychorismate_lyaseb <- aov(aminodeoxychorismate_lyase ~ antibiotic, data = Function2_d78)
summary(m3_aminodeoxychorismate_lyaseb)
TukeyHSD(m3_aminodeoxychorismate_lyaseb)

#shikimate_dehydrogenase
#d7
m1_shikimate_dehydrogenase <- lm(shikimate_dehydrogenase ~ antibiotic, data = Function2_d7)
summary(m1_shikimate_dehydrogenase)
m1_shikimate_dehydrogenaseb <- aov(shikimate_dehydrogenase ~ antibiotic, data = Function2_d7)
summary(m1_shikimate_dehydrogenaseb)
TukeyHSD(m1_shikimate_dehydrogenaseb)
#d35
m2_shikimate_dehydrogenase <- lm(shikimate_dehydrogenase ~ antibiotic, data = Function2_d35)
summary(m2_shikimate_dehydrogenase)
m2_shikimate_dehydrogenaseb <- aov(shikimate_dehydrogenase ~ antibiotic, data = Function2_d35)
summary(m2_shikimate_dehydrogenaseb)
TukeyHSD(m2_shikimate_dehydrogenaseb)
#d78
m3_shikimate_dehydrogenase <- lm(shikimate_dehydrogenase ~ antibiotic, data = Function2_d78)
summary(m3_shikimate_dehydrogenase)
m3_shikimate_dehydrogenaseb <- aov(shikimate_dehydrogenase ~ antibiotic, data = Function2_d78)
summary(m3_shikimate_dehydrogenaseb)
TukeyHSD(m3_shikimate_dehydrogenaseb)

#shikimate_5_dehydrogenase
#d7
m1_shikimate_5_dehydrogenase <- lm(shikimate_5_dehydrogenase ~ antibiotic, data = Function2_d7)
summary(m1_shikimate_5_dehydrogenase)
m1_shikimate_5_dehydrogenaseb <- aov( shikimate_5_dehydrogenase~ antibiotic, data = Function2_d7)
summary(m1_shikimate_5_dehydrogenaseb)
TukeyHSD(m1_shikimate_5_dehydrogenaseb)
#d35
m2_shikimate_5_dehydrogenase <- lm(shikimate_5_dehydrogenase ~ antibiotic, data = Function2_d35)
summary(m2_shikimate_5_dehydrogenase)
m2_shikimate_5_dehydrogenaseb <- aov( shikimate_5_dehydrogenase~ antibiotic, data = Function2_d35)
summary(m2_shikimate_5_dehydrogenaseb)
TukeyHSD(m2_shikimate_5_dehydrogenaseb)
#d78
m3_shikimate_5_dehydrogenase <- lm(shikimate_5_dehydrogenase ~ antibiotic, data = Function2_d78)
summary(m3_shikimate_5_dehydrogenase)
m3_shikimate_5_dehydrogenaseb <- aov( shikimate_5_dehydrogenase~ antibiotic, data = Function2_d78)
summary(m3_shikimate_5_dehydrogenaseb)
TukeyHSD(m3_shikimate_5_dehydrogenaseb)

#3_phosphoshikimate_1_carboxyvinyltransferase
#d7
m1_3_phosphoshikimate_1_carboxyvinyltransferase <- lm(phosphoshikimate_1_carboxyvinyltransferase ~ antibiotic, data = Function2_d7)
summary(m1_3_phosphoshikimate_1_carboxyvinyltransferase)
m1_3_phosphoshikimate_1_carboxyvinyltransferaseb <- aov(phosphoshikimate_1_carboxyvinyltransferase~ antibiotic, data = Function2_d7)
summary(m1_3_phosphoshikimate_1_carboxyvinyltransferaseb)
TukeyHSD(m1_3_phosphoshikimate_1_carboxyvinyltransferaseb)
#d35
m2_3_phosphoshikimate_1_carboxyvinyltransferase <- lm(phosphoshikimate_1_carboxyvinyltransferase ~ antibiotic, data = Function2_d35)
summary(m2_3_phosphoshikimate_1_carboxyvinyltransferase)
m2_3_phosphoshikimate_1_carboxyvinyltransferaseb <- aov(phosphoshikimate_1_carboxyvinyltransferase~ antibiotic, data = Function2_d35)
summary(m2_3_phosphoshikimate_1_carboxyvinyltransferaseb)
TukeyHSD(m2_3_phosphoshikimate_1_carboxyvinyltransferaseb)
#d78
m3_3_phosphoshikimate_1_carboxyvinyltransferase <- lm(phosphoshikimate_1_carboxyvinyltransferase ~ antibiotic, data = Function2_d78)
summary(m3_3_phosphoshikimate_1_carboxyvinyltransferase)
m3_3_phosphoshikimate_1_carboxyvinyltransferaseb <- aov(phosphoshikimate_1_carboxyvinyltransferase~ antibiotic, data = Function2_d78)
summary(m3_3_phosphoshikimate_1_carboxyvinyltransferaseb)
TukeyHSD(m3_3_phosphoshikimate_1_carboxyvinyltransferaseb)

#PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family
#d7
m1_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family <- lm(membrane_bound_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family ~ antibiotic, data = Function2_d7)
summary(m1_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family)
m1_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_familyb <- aov(membrane_bound_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family ~ antibiotic, data = Function2_d7)
summary(m1_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_familyb)
TukeyHSD(m1_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_familyb)
#d35
m2_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family <- lm(membrane_bound_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family ~ antibiotic, data = Function2_d35)
summary(m2_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family)
m2_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_familyb <- aov(membrane_bound_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family ~ antibiotic, data = Function2_d35)
summary(m2_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_familyb)
TukeyHSD(m2_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_familyb)
#d78
m3_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family <- lm(membrane_bound_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family ~ antibiotic, data = Function2_d78)
summary(m3_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family)
m3_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_familyb <- aov(membrane_bound_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family ~ antibiotic, data = Function2_d78)
summary(m3_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_familyb)
TukeyHSD(m3_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_familyb)

#deoxy_7_phosphoheptulonate_synthase
#d7
m1_deoxy_7_phosphoheptulonate_synthase <- lm(deoxy_7_phosphoheptulonate_synthase ~ antibiotic, data = Function2_d7)
summary(m1_deoxy_7_phosphoheptulonate_synthase)
m1_deoxy_7_phosphoheptulonate_synthaseb <- aov(deoxy_7_phosphoheptulonate_synthase ~ antibiotic, data = Function2_d7)
summary(m1_deoxy_7_phosphoheptulonate_synthaseb)
TukeyHSD(m1_deoxy_7_phosphoheptulonate_synthaseb)
#d35
m2_deoxy_7_phosphoheptulonate_synthase <- lm(deoxy_7_phosphoheptulonate_synthase ~ antibiotic, data = Function2_d35)
summary(m2_deoxy_7_phosphoheptulonate_synthase)
m2_deoxy_7_phosphoheptulonate_synthaseb <- aov(deoxy_7_phosphoheptulonate_synthase ~ antibiotic, data = Function2_d35)
summary(m2_deoxy_7_phosphoheptulonate_synthaseb)
TukeyHSD(m2_deoxy_7_phosphoheptulonate_synthaseb)
#d78
m3_deoxy_7_phosphoheptulonate_synthase <- lm(deoxy_7_phosphoheptulonate_synthase ~ antibiotic, data = Function2_d78)
summary(m3_deoxy_7_phosphoheptulonate_synthase)
m3_deoxy_7_phosphoheptulonate_synthaseb <- aov(deoxy_7_phosphoheptulonate_synthase ~ antibiotic, data = Function2_d78)
summary(m3_deoxy_7_phosphoheptulonate_synthaseb)
TukeyHSD(m3_deoxy_7_phosphoheptulonate_synthaseb)

#indolepyruvate_ferredoxin_oxidoreductase_beta_subunit
#d7
m1_indolepyruvate_ferredoxin_oxidoreductase_beta_subunit <- lm(indolepyruvate_ferredoxin_oxidoreductase_beta_subunit ~ antibiotic, data = Function2_d7)
summary(m1_indolepyruvate_ferredoxin_oxidoreductase_beta_subunit)
m1_indolepyruvate_ferredoxin_oxidoreductase_beta_subunitb <- aov(indolepyruvate_ferredoxin_oxidoreductase_beta_subunit ~ antibiotic, data = Function2_d7)
summary(m1_indolepyruvate_ferredoxin_oxidoreductase_beta_subunitb)
TukeyHSD(m1_indolepyruvate_ferredoxin_oxidoreductase_beta_subunitb)
#d35
m2_indolepyruvate_ferredoxin_oxidoreductase_beta_subunit <- lm(indolepyruvate_ferredoxin_oxidoreductase_beta_subunit ~ antibiotic, data = Function2_d35)
summary(m2_indolepyruvate_ferredoxin_oxidoreductase_beta_subunit)
m2_indolepyruvate_ferredoxin_oxidoreductase_beta_subunitb <- aov(indolepyruvate_ferredoxin_oxidoreductase_beta_subunit ~ antibiotic, data = Function2_d35)
summary(m2_indolepyruvate_ferredoxin_oxidoreductase_beta_subunitb)
TukeyHSD(m2_indolepyruvate_ferredoxin_oxidoreductase_beta_subunitb)
#d78
m3_indolepyruvate_ferredoxin_oxidoreductase_beta_subunit <- lm(indolepyruvate_ferredoxin_oxidoreductase_beta_subunit ~ antibiotic, data = Function2_d78)
summary(m3_indolepyruvate_ferredoxin_oxidoreductase_beta_subunit)
m3_indolepyruvate_ferredoxin_oxidoreductase_beta_subunitb <- aov(indolepyruvate_ferredoxin_oxidoreductase_beta_subunit ~ antibiotic, data = Function2_d78)
summary(m3_indolepyruvate_ferredoxin_oxidoreductase_beta_subunitb)
TukeyHSD(m3_indolepyruvate_ferredoxin_oxidoreductase_beta_subunitb)

#indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit
#d7
m1_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit <- lm(indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit ~ antibiotic, data = Function2_d7)
summary(m1_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit)
m1_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunitb <- aov(indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit ~ antibiotic, data = Function2_d7)
summary(m1_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunitb)
TukeyHSD(m1_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunitb)
#d35
m2_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit <- lm(indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit ~ antibiotic, data = Function2_d35)
summary(m2_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit)
m2_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunitb <- aov(indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit ~ antibiotic, data = Function2_d35)
summary(m2_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunitb)
TukeyHSD(m2_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunitb)
#d78
m3_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit <- lm(indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit ~ antibiotic, data = Function2_d78)
summary(m3_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit)
m3_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunitb <- aov(indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit ~ antibiotic, data = Function2_d78)
summary(m3_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunitb)
TukeyHSD(m3_indolepyruvate_ferredoxin_oxidoreductase_alpha_subunitb)

#phenylacetate_CoA_oxygenase_PaaG_subunit
#d7
m1_phenylacetate_CoA_oxygenase_PaaG_subunit <- lm(phenylacetate_CoA_oxygenase_PaaG_subunit ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetate_CoA_oxygenase_PaaG_subunit)
m1_phenylacetate_CoA_oxygenase_PaaG_subunitb <- aov(phenylacetate_CoA_oxygenase_PaaG_subunit ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetate_CoA_oxygenase_PaaG_subunitb)
TukeyHSD(m1_phenylacetate_CoA_oxygenase_PaaG_subunitb)
#d35
m2_phenylacetate_CoA_oxygenase_PaaG_subunit <- lm(phenylacetate_CoA_oxygenase_PaaG_subunit ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetate_CoA_oxygenase_PaaG_subunit)
m2_phenylacetate_CoA_oxygenase_PaaG_subunitb <- aov(phenylacetate_CoA_oxygenase_PaaG_subunit ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetate_CoA_oxygenase_PaaG_subunitb)
TukeyHSD(m2_phenylacetate_CoA_oxygenase_PaaG_subunitb)
#d78
m3_phenylacetate_CoA_oxygenase_PaaG_subunit <- lm(phenylacetate_CoA_oxygenase_PaaG_subunit ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetate_CoA_oxygenase_PaaG_subunit)
m3_phenylacetate_CoA_oxygenase_PaaG_subunitb <- aov(phenylacetate_CoA_oxygenase_PaaG_subunit ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetate_CoA_oxygenase_PaaG_subunitb)
TukeyHSD(m3_phenylacetate_CoA_oxygenase_PaaG_subunitb)

#phenylacetate_CoA_oxygenase_PaaH_subunit
#d7
m1_phenylacetate_CoA_oxygenase_PaaH_subunit <- lm(phenylacetate_CoA_oxygenase_PaaH_subunit ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetate_CoA_oxygenase_PaaH_subunit)
m1_phenylacetate_CoA_oxygenase_PaaH_subunitb <- aov(phenylacetate_CoA_oxygenase_PaaH_subunit ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetate_CoA_oxygenase_PaaH_subunitb)
TukeyHSD(m1_phenylacetate_CoA_oxygenase_PaaH_subunitb)
#d35
m2_phenylacetate_CoA_oxygenase_PaaH_subunit <- lm(phenylacetate_CoA_oxygenase_PaaH_subunit ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetate_CoA_oxygenase_PaaH_subunit)
m2_phenylacetate_CoA_oxygenase_PaaH_subunitb <- aov(phenylacetate_CoA_oxygenase_PaaH_subunit ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetate_CoA_oxygenase_PaaH_subunitb)
TukeyHSD(m2_phenylacetate_CoA_oxygenase_PaaH_subunitb)
#d78
m3_phenylacetate_CoA_oxygenase_PaaH_subunit <- lm(phenylacetate_CoA_oxygenase_PaaH_subunit ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetate_CoA_oxygenase_PaaH_subunit)
m3_phenylacetate_CoA_oxygenase_PaaH_subunitb <- aov(phenylacetate_CoA_oxygenase_PaaH_subunit ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetate_CoA_oxygenase_PaaH_subunitb)
TukeyHSD(m3_phenylacetate_CoA_oxygenase_PaaH_subunitb)

#phenylacetate_CoA_oxygenase_PaaI_subunit
#d7
m1_phenylacetate_CoA_oxygenase_PaaI_subunit <- lm(phenylacetate_CoA_oxygenase_PaaI_subunit ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetate_CoA_oxygenase_PaaI_subunit)
m1_phenylacetate_CoA_oxygenase_PaaI_subunitb <- aov(phenylacetate_CoA_oxygenase_PaaI_subunit ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetate_CoA_oxygenase_PaaI_subunitb)
TukeyHSD(m1_phenylacetate_CoA_oxygenase_PaaI_subunitb)
#d35
m2_phenylacetate_CoA_oxygenase_PaaI_subunit <- lm(phenylacetate_CoA_oxygenase_PaaI_subunit ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetate_CoA_oxygenase_PaaI_subunit)
m2_phenylacetate_CoA_oxygenase_PaaI_subunitb <- aov(phenylacetate_CoA_oxygenase_PaaI_subunit ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetate_CoA_oxygenase_PaaI_subunitb)
TukeyHSD(m2_phenylacetate_CoA_oxygenase_PaaI_subunitb)
#d78
m3_phenylacetate_CoA_oxygenase_PaaI_subunit <- lm(phenylacetate_CoA_oxygenase_PaaI_subunit ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetate_CoA_oxygenase_PaaI_subunit)
m3_phenylacetate_CoA_oxygenase_PaaI_subunitb <- aov(phenylacetate_CoA_oxygenase_PaaI_subunit ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetate_CoA_oxygenase_PaaI_subunitb)
TukeyHSD(m3_phenylacetate_CoA_oxygenase_PaaI_subunitb)

#phenylacetate_CoA_oxygenase_PaaJ_subunit
#d7
m1_phenylacetate_CoA_oxygenase_PaaJ_subunit <- lm(phenylacetate_CoA_oxygenase_PaaJ_subunit ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetate_CoA_oxygenase_PaaJ_subunit)
m1_phenylacetate_CoA_oxygenase_PaaJ_subunitb <- aov(phenylacetate_CoA_oxygenase_PaaJ_subunit ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetate_CoA_oxygenase_PaaJ_subunitb)
TukeyHSD(m1_phenylacetate_CoA_oxygenase_PaaJ_subunitb)
#d35
m2_phenylacetate_CoA_oxygenase_PaaJ_subunit <- lm(phenylacetate_CoA_oxygenase_PaaJ_subunit ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetate_CoA_oxygenase_PaaJ_subunit)
m2_phenylacetate_CoA_oxygenase_PaaJ_subunitb <- aov(phenylacetate_CoA_oxygenase_PaaJ_subunit ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetate_CoA_oxygenase_PaaJ_subunitb)
TukeyHSD(m2_phenylacetate_CoA_oxygenase_PaaJ_subunitb)
#d78
m3_phenylacetate_CoA_oxygenase_PaaJ_subunit <- lm(phenylacetate_CoA_oxygenase_PaaJ_subunit ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetate_CoA_oxygenase_PaaJ_subunit)
m3_phenylacetate_CoA_oxygenase_PaaJ_subunitb <- aov(phenylacetate_CoA_oxygenase_PaaJ_subunit ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetate_CoA_oxygenase_PaaJ_subunitb)
TukeyHSD(m3_phenylacetate_CoA_oxygenase_PaaJ_subunitb)

#phenylacetate_CoA_oxygenase_reductase_PaaK_subunit
#d7
m1_phenylacetate_CoA_oxygenase_reductase_PaaK_subunit <- lm(phenylacetate_CoA_oxygenase_reductase_PaaK_subunit ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetate_CoA_oxygenase_reductase_PaaK_subunit)
m1_phenylacetate_CoA_oxygenase_reductase_PaaK_subunitb <- aov(phenylacetate_CoA_oxygenase_reductase_PaaK_subunit ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetate_CoA_oxygenase_reductase_PaaK_subunitb)
TukeyHSD(m1_phenylacetate_CoA_oxygenase_reductase_PaaK_subunitb)
#d35
m2_phenylacetate_CoA_oxygenase_reductase_PaaK_subunit <- lm(phenylacetate_CoA_oxygenase_reductase_PaaK_subunit ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetate_CoA_oxygenase_reductase_PaaK_subunit)
m2_phenylacetate_CoA_oxygenase_reductase_PaaK_subunitb <- aov(phenylacetate_CoA_oxygenase_reductase_PaaK_subunit ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetate_CoA_oxygenase_reductase_PaaK_subunitb)
TukeyHSD(m2_phenylacetate_CoA_oxygenase_reductase_PaaK_subunitb)
#d78
m3_phenylacetate_CoA_oxygenase_reductase_PaaK_subunit <- lm(phenylacetate_CoA_oxygenase_reductase_PaaK_subunit ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetate_CoA_oxygenase_reductase_PaaK_subunit)
m3_phenylacetate_CoA_oxygenase_reductase_PaaK_subunitb <- aov(phenylacetate_CoA_oxygenase_reductase_PaaK_subunit ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetate_CoA_oxygenase_reductase_PaaK_subunitb)
TukeyHSD(m3_phenylacetate_CoA_oxygenase_reductase_PaaK_subunitb)

#phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB
#d7
m1_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB <- lm(phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB)
m1_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaBb <- aov(phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaBb)
TukeyHSD(m1_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaBb)
#d35
m2_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB <- lm(phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB)
m2_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaBb <- aov(phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaBb)
TukeyHSD(m2_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaBb)
#d78
m3_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB <- lm(phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB)
m3_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaBb <- aov(phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaBb)
TukeyHSD(m3_phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaBb)

#phenylacetate_CoA_ligase
#d7
m1_phenylacetate_CoA_ligase <- lm(phenylacetate_CoA_ligase ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetate_CoA_ligase)
m1_phenylacetate_CoA_ligaseb <- aov(phenylacetate_CoA_ligase ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetate_CoA_ligaseb)
TukeyHSD(m1_phenylacetate_CoA_ligaseb)
#d7
m2_phenylacetate_CoA_ligase <- lm(phenylacetate_CoA_ligase ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetate_CoA_ligase)
m2_phenylacetate_CoA_ligaseb <- aov(phenylacetate_CoA_ligase ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetate_CoA_ligaseb)
TukeyHSD(m2_phenylacetate_CoA_ligaseb)
#d7
m3_phenylacetate_CoA_ligase <- lm(phenylacetate_CoA_ligase ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetate_CoA_ligase)
m3_phenylacetate_CoA_ligaseb <- aov(phenylacetate_CoA_ligase ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetate_CoA_ligaseb)
TukeyHSD(m3_phenylacetate_CoA_ligaseb)

#phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX
#d7
m1_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX <- lm(phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX)
m1_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaXb <- aov(phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaXb)
TukeyHSD(m1_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaXb)
#d35
m2_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX <- lm(phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX)
m2_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaXb <- aov(phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaXb)
TukeyHSD(m2_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaXb)
#d78
m3_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX <- lm(phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX)
m3_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaXb <- aov(phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaXb)
TukeyHSD(m3_phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaXb)

#phenylacetic_acid_degradation_protein_PaaD
#d7
m1_phenylacetic_acid_degradation_protein_PaaD <- lm(phenylacetic_acid_degradation_protein_PaaD ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetic_acid_degradation_protein_PaaD)
m1_phenylacetic_acid_degradation_protein_PaaDb <- aov(phenylacetic_acid_degradation_protein_PaaD ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetic_acid_degradation_protein_PaaDb)
TukeyHSD(m1_phenylacetic_acid_degradation_protein_PaaDb)
#d35
m2_phenylacetic_acid_degradation_protein_PaaD <- lm(phenylacetic_acid_degradation_protein_PaaD ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetic_acid_degradation_protein_PaaD)
m2_phenylacetic_acid_degradation_protein_PaaDb <- aov(phenylacetic_acid_degradation_protein_PaaD ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetic_acid_degradation_protein_PaaDb)
TukeyHSD(m2_phenylacetic_acid_degradation_protein_PaaDb)
#d78
m3_phenylacetic_acid_degradation_protein_PaaD <- lm(phenylacetic_acid_degradation_protein_PaaD ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetic_acid_degradation_protein_PaaD)
m3_phenylacetic_acid_degradation_protein_PaaDb <- aov(phenylacetic_acid_degradation_protein_PaaD ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetic_acid_degradation_protein_PaaDb)
TukeyHSD(m3_phenylacetic_acid_degradation_protein_PaaDb)

#phenylacetic_acid_degradation_protein_paaN
#d7
m1_phenylacetic_acid_degradation_protein_paaN <- lm(phenylacetic_acid_degradation_protein_paaN ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetic_acid_degradation_protein_paaN)
m1_phenylacetic_acid_degradation_protein_paaNb <- aov(phenylacetic_acid_degradation_protein_paaN ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetic_acid_degradation_protein_paaNb)
TukeyHSD(m1_phenylacetic_acid_degradation_protein_paaNb)
#d7
m2_phenylacetic_acid_degradation_protein_paaN <- lm(phenylacetic_acid_degradation_protein_paaN ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetic_acid_degradation_protein_paaN)
m2_phenylacetic_acid_degradation_protein_paaNb <- aov(phenylacetic_acid_degradation_protein_paaN ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetic_acid_degradation_protein_paaNb)
TukeyHSD(m2_phenylacetic_acid_degradation_protein_paaNb)
#d7
m3_phenylacetic_acid_degradation_protein_paaN <- lm(phenylacetic_acid_degradation_protein_paaN ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetic_acid_degradation_protein_paaN)
m3_phenylacetic_acid_degradation_protein_paaNb <- aov(phenylacetic_acid_degradation_protein_paaN ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetic_acid_degradation_protein_paaNb)
TukeyHSD(m3_phenylacetic_acid_degradation_protein_paaNb)

#phenylacetic_acid_degradation_protein_PaaY
#d7
m1_phenylacetic_acid_degradation_protein_PaaY <- lm(phenylacetic_acid_degradation_protein_PaaY ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetic_acid_degradation_protein_PaaY)
m1_phenylacetic_acid_degradation_protein_PaaYb <- aov(phenylacetic_acid_degradation_protein_PaaY ~ antibiotic, data = Function2_d7)
summary(m1_phenylacetic_acid_degradation_protein_PaaYb)
TukeyHSD(m1_phenylacetic_acid_degradation_protein_PaaYb)
#d35
m2_phenylacetic_acid_degradation_protein_PaaY <- lm(phenylacetic_acid_degradation_protein_PaaY ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetic_acid_degradation_protein_PaaY)
m2_phenylacetic_acid_degradation_protein_PaaYb <- aov(phenylacetic_acid_degradation_protein_PaaY ~ antibiotic, data = Function2_d35)
summary(m2_phenylacetic_acid_degradation_protein_PaaYb)
TukeyHSD(m2_phenylacetic_acid_degradation_protein_PaaYb)
#d78
m3_phenylacetic_acid_degradation_protein_PaaY <- lm(phenylacetic_acid_degradation_protein_PaaY ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetic_acid_degradation_protein_PaaY)
m3_phenylacetic_acid_degradation_protein_PaaYb <- aov(phenylacetic_acid_degradation_protein_PaaY ~ antibiotic, data = Function2_d78)
summary(m3_phenylacetic_acid_degradation_protein_PaaYb)
TukeyHSD(m3_phenylacetic_acid_degradation_protein_PaaYb)

ggplot(Function2, aes(x=day, y=tryptophan_synthase_alpha_subunit_, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Tryptophan synthase alpha subunit") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=23), axis.title.y = element_text(color="black", size=23)) + 
  theme(axis.text.x = element_text(color = "black", size = 18), axis.text.y = element_text(color = "black", size = 18))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output2/Tryptophan_syhthase_alpha", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=tryptophan_synthase_beta_subunit_, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Tryptophan synthase beta subunit") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=23), axis.title.y = element_text(color="black", size=23)) + 
  theme(axis.text.x = element_text(color = "black", size = 18), axis.text.y = element_text(color = "black", size = 18))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output2/Tryptophan_syhthase_beta", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=glutamine_amidotransferase_anthranilate_synthase_aminodeoxychorismate_synthase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Aminodeoxychorismate synthase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=23), axis.title.y = element_text(color="black", size=23)) + 
  theme(axis.text.x = element_text(color = "black", size = 18), axis.text.y = element_text(color = "black", size = 18))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output2/aminodeoxychorismate_synthase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=anthranilate_synthase_component_I, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Anthranilate synthase component I") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=23), axis.title.y = element_text(color="black", size=23)) + 
  theme(axis.text.x = element_text(color = "black", size = 18), axis.text.y = element_text(color = "black", size = 18))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output2/anthranilate_synthase_component_I", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=anthranilate_phosphoribosyltransferase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Anthranilate phosphoribosyltransferase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=22.4), axis.title.y = element_text(color="black", size=22.4)) + 
  theme(axis.text.x = element_text(color = "black", size = 18), axis.text.y = element_text(color = "black", size = 18))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output2/anthranilate_phosphoribosyltransferase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Function2, aes(x=day, y=chorismate_mutase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Chorismate mutase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=23), axis.title.y = element_text(color="black", size=23)) + 
  theme(axis.text.x = element_text(color = "black", size = 18), axis.text.y = element_text(color = "black", size = 18))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output2/chorismate_mutase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Function2, aes(x=day, y=chorismate_synthase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Chorismate synthase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=23), axis.title.y = element_text(color="black", size=23)) + 
  theme(axis.text.x = element_text(color = "black", size = 18), axis.text.y = element_text(color = "black", size = 18))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output2/chorismate_synthase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Function2, aes(x=day, y=isochorismate_synthase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Isochorismate synthase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/isochorismate_synthase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Function2, aes(x=day, y=aminodeoxychorismate_lyase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Aminodeoxychorismate lyase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/aminodeoxychorismate_lyase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Function2, aes(x=day, y=shikimate_dehydrogenase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Shikimate dehydrogenase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=23), axis.title.y = element_text(color="black", size=23)) + 
  theme(axis.text.x = element_text(color = "black", size = 18), axis.text.y = element_text(color = "black", size = 18))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output2/shikimate_dehydrogenase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Function2, aes(x=day, y=shikimate_5_dehydrogenase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Shikimate-5-dehydrogenase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/shikimate_5_dehydrogenase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Function2, aes(x=day, y=phosphoshikimate_1_carboxyvinyltransferase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("3-Phosphoshikimate-1-carboxyvinyltransferase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/phosphoshikimate_1_carboxyvinyltransferase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Function2, aes(x=day, y=membrane_bound_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("membrane bound PQQ dependent dehydrogenase glucose quinate shikimate family") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/membrane_bound_PQQ_dependent_dehydrogenase_glucose_quinate_shikimate_family", ".pdf"), height=7.1, width=9, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Function2, aes(x=day, y=deoxy_7_phosphoheptulonate_synthase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("3-deoxy-7-phosphoheptulonate_synthase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/deoxy_7_phosphoheptulonate_synthase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

#Not used for the manuscript
ggplot(Function2, aes(x=day, y=imidazole_glycerol_phosphate_synthase_glutamine_amidotransferase_subunit, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("glutamine amidotransferase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/glutamine_amidotransferase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Indolepyruvate ferredoxin oxidoreductase alpha subunit") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/indolepyruvate_ferredoxin_oxidoreductase_alpha_subunit", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=indolepyruvate_ferredoxin_oxidoreductase_beta_subunit, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Indolepyruvate ferredoxin oxidoreductase beta subunit") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/indolepyruvate_ferredoxin_oxidoreductase_beta_subunit", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=phenylacetate_CoA_oxygenase_PaaG_subunit, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Phenylacetate CoA oxygenase PaaG subunit") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/phenylacetate_CoA_oxygenase_PaaG_subunit", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=phenylacetate_CoA_oxygenase_PaaH_subunit, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Phenylacetate CoA oxygenase PaaH subunit") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/phenylacetate_CoA_oxygenase_PaaH_subunit", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=phenylacetate_CoA_oxygenase_PaaI_subunit, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Phenylacetate CoA oxygenase PaaI subunit") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/phenylacetate_CoA_oxygenase_PaaI_subunit", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=phenylacetate_CoA_oxygenase_PaaJ_subunit, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Phenylacetate CoA oxygenase PaaJ subunit") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/phenylacetate_CoA_oxygenase_PaaJ_subunit", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=phenylacetate_CoA_oxygenase_reductase_PaaK_subunit, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Phenylacetate CoA oxygenase reductase PaaK subunit") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/phenylacetate_CoA_oxygenase_reductase_PaaK_subunit", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

#Not included in manuscript 
ggplot(Function2, aes(x=day, y=acyl_CoA_thioester_hydrolase_YbgC_YbaW_family, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("acyl CoA thioester hydrolase YbgC/YbaW family") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/acyl_CoA_thioester_hydrolase_YbgC_YbaW_family", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

#Not included in manuscript 
ggplot(Function2, aes(x=day, y=alpha_L_glutamate_ligase_RimK_family, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("alpha L glutamate ligase RimK family") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/alpha_L_glutamate_ligase_RimK_family", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

#Not included in manuscript 
ggplot(Function2, aes(x=day, y=glutamate_synthase_NADPH_homotetrameric, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("glutamate synthase (NADPH) homotetrameric") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/glutamate synthase_(NADPH)_homotetrameric", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

#Not included in manuscript 
ggplot(Function2, aes(x=day, y=glutamate_synthase_NADH_NADPH_small_subunit, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("glutamate synthase NADH NADPH small subunit") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/glutamate_synthase_NADH_NADPH_small_subunit", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

#Not included in manuscript 
ggplot(Function2, aes(x=day, y=glutamate_synthase_small_subunit, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("glutamate synthase small subunit") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/glutamate_synthase_small_subunit", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

#Not included in manuscript 
ggplot(Function2, aes(x=day, y=glutamate_1_semialdehyde_2_1_aminomutase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("glutamate 1 semialdehyde 2 1 aminomutase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/glutamate_1_semialdehyde_2_1_aminomutase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

#Not included in manuscript 
ggplot(Function2, aes(x=day, y=glutamine_synthetase_type_I, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("glutamine synthetase type I") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/glutamine_synthetase_type_I", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

#Not included in manuscript 
ggplot(Function2, aes(x=day, y=glyceraldehyde_3_phosphate_dehydrogenase_type_I, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("glyceraldehyde-3-phosphate dehydrogenase type I") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/glyceraldehyde_3_phosphate_dehydrogenase_type_I", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Phenylacetate degradation probable enoyl CoA hydratase PaaB") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/phenylacetate_degradation_probable_enoyl_CoA_hydratase_PaaB", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=phenylacetate_CoA_ligase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Phenylacetate CoA ligase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/phenylacetate_CoA_ligase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Phenylacetic acid degradation operon negative regulatory protein PaaX") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/phenylacetic_acid_degradation_operon_negative_regulatory_protein_PaaX", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=phenylacetic_acid_degradation_protein_PaaD, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Phenylacetic acid degradation protein PaaD") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/phenylacetic_acid_degradation_protein_PaaD", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=phenylacetic_acid_degradation_protein_paaN, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Phenylacetic acid degradation protein PaaN") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/phenylacetic_acid_degradation_protein_paaN", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=phenylacetic_acid_degradation_protein_PaaY, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Phenylacetic acid degradation protein PaaY") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/phenylacetic_acid_degradation_protein_PaaY", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=poly_gamma_glutamate_biosynthesis_protein_PgsC, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("poly gamma glutamate biosynthesis protein PgsC") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/poly-gamma_glutamate_biosynthesis_protein_PgsC", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=poly_gamma_glutamate_synthase_PgsB, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("poly gamma glutamate synthase PgsB") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/poly_gamma_glutamate_synthase_PgsB", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=poly_gamma_glutamate_system_protein, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("poly gamma glutamate system protein") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/poly_gamma_glutamate_system_protein", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=putative_glutamate_gamma_aminobutyrate_antiporter, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("putative glutamate gamma aminobutyrate antiporter") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/putative_glutamate_gamma_aminobutyrate_antiporter", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=sodium_glutamate_symporter, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("sodium glutamate symporter") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/sodium_glutamate_symporter", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=succinylglutamate_desuccinylase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("succinylglutamate desuccinylase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/succinylglutamate_desuccinylase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=succinylglutamate_semialdehyde_dehydrogenase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("succinylglutamate semialdehyde dehydrogenase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/succinylglutamate_semialdehyde_dehydrogenase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=tryptophan_tRNA_ligase, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("tryptophan tRNA ligase") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/tryptophan_tRNA_ligase", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=Ni_specific_regulatory_protein, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Ni specific regulatory protein") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/Ni_specific_regulatory_protein", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=nif11_like_leader_peptide_domain, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("nif11 like leader peptide domain") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/nif11_like_leader_peptide_domain", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=nitrate_reductase_molybdenum_cofactor_assembly_chaperone, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("nitrate reductase molybdenum cofactor assembly chaperone") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/nitrate_reductase_molybdenum_cofactor_assembly_chaperone", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=nitrite_transporter, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("nitrite transporter") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/nitrite_transporter", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=nitrogen_regulation_protein_N_I, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("nitrogen regulation protein N(I)") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/nitrogen_regulation_protein_N_I", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=nitrogenase_cofactor_biosynthesis_protein_NifB, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("nitrogenase cofactor biosynthesis protein NifB") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/nitrogenase_cofactor_biosynthesis_protein_NifB", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=nitrogenase_iron_protein, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("nitrogenase iron protein") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/nitrogenase_iron_protein", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=nitrogenase_MoFe_cofactor_biosynthesis_protein_NifE, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("nitrogenase MoFe cofactor biosynthesis protein NifE") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/nitrogenase_MoFe_cofactor_biosynthesis_protein_NifE", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=nitrogenase_molybdenum_iron_protein_alpha_chain, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("nitrogenase molybdenum iron protein alpha chain") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/nitrogenase_molybdenum_iron_protein_alpha_chain", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=NusG_family_protein, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("NusG family protein") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/NusG_family_protein", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=protein_RecA, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Protein RecA") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/protein_RecA", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=glycosyltransferase_WecB_TagA_CpsF_family, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("glycosyltransferase WecB TagA CpsF family") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/glycosyltransferase_WecB_TagA_CpsF_family", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=Holliday_junction_DNA_helicase_RuvA, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Holliday junction DNA helicase RuvA") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/Holliday_junction_DNA_helicase_RuvA", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=Holliday_junction_DNA_helicase_RuvB, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("Holliday junction DNA helicase RuvB") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/Holliday_junction_DNA_helicase_RuvB", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=DNA_binding_regulatory_protein_YebC_PmpR_family, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("DNA_binding regulatory protein YebC PmpR family") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/DNA_binding_regulatory_protein_YebC_PmpR_family", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=repressor_LexA, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("repressor LexA") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/repressor_LexA", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=transcription_termination_factor_NusA, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("transcription termination factor NusA") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/transcription_termination_factor_NusA", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=transcription_termination_antitermination_factor_NusG, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("transcription termination antitermination factor NusG") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/transcription_termination_antitermination_factor_NusG", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(Function2, aes(x=day, y=, color=antibiotic)) + 
  geom_boxplot() + 
  theme_bw()+
  ylab("") + xlab ("Day") + 
  scale_color_manual(values = my_colors) +
  labs(color= "Antibiotic") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
#theme(text=element_text(family="Times New Roman"))
ggsave(paste0("output/", ".pdf"), height=6, width=8, device="pdf") # save a PDF 3 inches by 4 inches

