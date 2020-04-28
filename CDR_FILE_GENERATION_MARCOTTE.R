require(dplyr)
require(tidyr)


GFP <- fread("Unique_Combined_CDRs_GFP2_nucleotides_and_peptides.csv")

AB <- fread("Unique_Combined_CDRs_Ab2_nucleotides_and_peptides.csv")

#Input <- fread("Unique_Combined_CDRs_Input_nucleotides_and_peptides.csv")

#GFP <- filter(GFP, INP1 > 1)
AB <- filter(AB, INP1 > 1)
stops <- rename(stops, combined = V1)
AB <- anti_join(AB, stops, by = c("combined"))

#NANOBODIES FROM AB THAT ARE *NOT* IN GFP CONTROL

Ab_s_pep1 <- anti_join(AB, GFP, by = c("peptide1")) 
Ab_s_pep2 <- anti_join(AB, GFP, by = c("peptide2")) 
Ab_s_pep3 <- anti_join(AB, GFP, by = c("peptide3")) 
Ab_s_combined <- anti_join(AB, GFP, by = c("combined"))


#NANOBODIES PRESENT IN GFP CONTROL
TU_pep1 <- semi_join(AB, GFP, by = c("peptide1")) 
TU_pep2 <- semi_join(AB, GFP, by = c("peptide2")) 
TU_pep3 <- semi_join(AB, GFP, by = c("peptide3")) 
TU_combined <- semi_join(AB, GFP, by = c("combined"))

#NANOBODIES FROM GFP THAT ARE *NOT* IN AB

GFP_s_pep1 <- anti_join(GFP, AB, by = c("peptide1")) 
GFP_s_pep2 <- anti_join(GFP, AB, by = c("peptide2")) 
GFP_s_pep3 <- anti_join(GFP, AB, by = c("peptide3")) 
GFP_s_combined <- anti_join(GFP, AB, by = c("combined"))

#FILTER FOR INP > 1

#Ab_s_pep1 <- filter(Ab_s_pep1, INP1 > 1)
#Ab_s_pep2 <- filter(Ab_s_pep2, INP1 > 1)
#Ab_s_pep3 <- filter(Ab_s_pep3, INP1 > 1)
#Ab_s_combined <- filter(Ab_s_combined, INP1 > 1)

#TU_pep1 <- filter(TU_pep1, INP1 > 1)
#TU_pep2 <- filter(TU_pep2, INP1 > 1)
#TU_pep3 <- filter(TU_pep3, INP1 > 1)

#GFP_s_pep1 <- filter(GFP_s_pep1, INP1 > 1)
#GFP_s_pep2 <- filter(GFP_s_pep2, INP1 > 1)
#GFP_s_pep3 <- filter(GFP_s_pep3, INP1 > 1)
#GFP_s_combined <- filter(GFP_s_combined, INP1 > 1)

#WRITE RELEVANT DATA TO FILES

Ab_combined <- select(Ab_s_combined, combined)
Ab_pep1 <- select(Ab_s_pep1, peptide1)
Ab_pep2 <- select(Ab_s_pep2, peptide2)
Ab_pep3 <- select(Ab_s_pep3, peptide3)

Ab_combined <- distinct(Ab_combined)
Ab_pep1 <- distinct(Ab_pep1)
Ab_pep2 <- distinct(Ab_pep2)
Ab_pep3 <- distinct(Ab_pep3)

write.table(Ab_combined, file = "Ab_s_combined_CDRs.txt", row.names = FALSE, col.names = FALSE)
write.table(Ab_pep1, file = "Ab_s_CDR1.txt", row.names = FALSE, col.names = FALSE)
write.table(Ab_pep2, file = "Ab_s_CDR2.txt", row.names = FALSE, col.names = FALSE)
write.table(Ab_pep3, file = "Ab_s_CDR3.txt", row.names = FALSE, col.names = FALSE)

#GFP_combined <- select(GFP_s_combined, combined)
#GFP_pep1 <- select(GFP_s_pep1, peptide1)
#GFP_pep2 <- select(GFP_s_pep2, peptide2)
#GFP_pep3 <- select(GFP_s_pep3, peptide3)

#write.table(GFP_combined, file = "GFP_s_combined_CDRs.txt", row.names = FALSE, col.names = FALSE)
#write.table(GFP_pep1, file = "GFP_s_CDR1.txt", row.names = FALSE, col.names = FALSE)
#write.table(GFP_pep2, file = "GFP_s_CDR2.txt", row.names = FALSE, col.names = FALSE)
#write.table(GFP_pep3, file = "GFP_s_CDR3.txt", row.names = FALSE, col.names = FALSE)

TU_combined <- select(TU_combined, combined)
TU_pep1 <- select(TU_pep1, peptide1)
TU_pep2 <- select(TU_pep2, peptide2)
TU_pep3 <- select(TU_pep3, peptide3)

TU_combined <- distinct(TU_combined)
TU_pep1 <- distinct(TU_pep1)
TU_pep2 <- distinct(TU_pep2)
TU_pep3 <- distinct(TU_pep3)


write.table(TU_combined, file = "TU_combined_CDRs.txt", row.names = FALSE, col.names = FALSE)
write.table(TU_pep1, file = "TU_CDR1.txt", row.names = FALSE, col.names = FALSE)
write.table(TU_pep2, file = "TU_CDR2.txt", row.names = FALSE, col.names = FALSE)
write.table(TU_pep3, file = "TU_CDR3.txt", row.names = FALSE, col.names = FALSE)


#CREATE FILES FOR CLUSTERING

CDR1_Ab.TU <- bind_rows(Ab_pep1, TU_pep1)
CDR2_Ab.TU <- bind_rows(Ab_pep2, TU_pep2)
CDR3_Ab.TU <- bind_rows(Ab_pep3, TU_pep3)
Combined_CDR_Ab.TU <- bind_rows(Ab_combined, TU_combined)

CDR1_AbTU <- distinct(CDR1_Ab.TU)
CDR2_AbTU <- distinct(CDR2_Ab.TU)
CDR3_AbTU <- distinct(CDR3_Ab.TU)
Combined_CDR_AbTU <- distinct(Combined_CDR_Ab.TU)


write.table(CDR1_AbTU, file = "CDR1_AbTU.txt", row.names = FALSE, col.names = FALSE)
write.table(CDR2_AbTU, file = "CDR2_AbTU.txt", row.names = FALSE, col.names = FALSE)
write.table(CDR3_AbTU, file = "CDR3_AbTU.txt", row.names = FALSE, col.names = FALSE)
write.table(Combined_CDR_AbTU, file = "Combined_CDR_AbTU.txt", row.names = FALSE, col.names = FALSE)

