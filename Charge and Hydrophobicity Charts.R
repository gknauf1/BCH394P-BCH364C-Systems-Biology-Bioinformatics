require(data.table)
require(seqinr)
require(stringi)
require(stringr)
require(dplyr)
require(tidyr)
require(Peptides)
require(ggplot2)
require(cowplot)
theme_set(theme_cowplot())


#rename
Ab_s_CDR1 <- rename(Ab_s_CDR1, Ab = V1)
Ab_s_CDR2 <- rename(Ab_s_CDR2, Ab = V1)
Ab_s_CDR3 <- rename(Ab_s_CDR3, Ab = V1)
Ab_s_combined_CDRs <- rename(Ab_s_combined_CDRs, Ab = V1)

TU_CDR1 <- rename(TU_CDR1, TU = V1)
TU_CDR2 <- rename(TU_CDR2, TU = V1)
TU_CDR3 <- rename(TU_CDR3, TU = V1)
TU_combined_CDRs <- rename(TU_combined_CDRs, TU = V1)

#combined CDR
Ab_s_combined_CDRs <- mutate(Ab_s_combined_CDRs, Hydrophob.Ab = hydrophobicity(Ab_s_combined_CDRs$Ab))

Ab_s_combined_CDRs <- mutate(Ab_s_combined_CDRs, Charge.Ab = charge(Ab_s_combined_CDRs$Ab, pH = 7.2))

TU_combined_CDRs <- mutate(TU_combined_CDRs, Hydrophob.TU = hydrophobicity(TU_combined_CDRs$TU))

TU_combined_CDRs <- mutate(TU_combined_CDRs, Charge.TU = charge(TU_combined_CDRs$TU, pH = 7.2))

ggplot(Ab_s_combined_CDRs, aes(Hydrophob.Ab)) + geom_histogram(bins = 68, colour = "black", fill = "lightblue") + 
  geom_vline(aes(xintercept = mean(Hydrophob.Ab)), size = 1, colour = "black",linetype = "dashed") +
  scale_x_continuous(breaks =c(-2, -1.5, -1.0, -0.5, 0, 0.5) ) +
  labs(title = "Ab Specific Nanobody Combined CDR Hydrophophobicities", x = "Hydrophobicity", y = "Count")

ggplot(Ab_s_combined_CDRs, aes(Charge.Ab)) + geom_histogram(bins = 10, colour = "black", fill = "pink") + 
  geom_vline(aes(xintercept = mean(Charge.Ab)), size = 1, colour = "black",linetype = "dashed") +
  labs(title = "Ab Specific Nanobody Combined CDR Charges", x = "Charge", y = "Count")

ggplot(TU_combined_CDRs, aes(Hydrophob.TU)) + geom_histogram(bins = 68, colour = "black", fill = "lightblue") + 
  geom_vline(aes(xintercept = mean(Hydrophob.TU)), size = 1, colour = "black",linetype = "dashed") +
  scale_x_continuous(breaks =c(-2, -1.5, -1.0, -0.5, 0, 0.5) ) +
  labs(title = "TU Nanobody Combined CDR Hydrophophobicities", x = "Hydrophobicity", y ="Count")

ggplot(TU_combined_CDRs, aes(Charge.TU)) + geom_histogram(bins = 10, colour = "black", fill = "pink") + 
  geom_vline(aes(xintercept = mean(Charge.TU)), size = 1, colour = "black",linetype = "dashed") +
  labs(title = "TU Nanobody Combined CDR Charges", x = "Charge", y = "Count")

#CDR1
Ab_s_CDR1 <- mutate(Ab_s_CDR1, Hydrophob.Ab = hydrophobicity(Ab_s_CDR1$Ab))

Ab_s_CDR1 <- mutate(Ab_s_CDR1, Charge.Ab = charge(Ab_s_CDR1$Ab, pH = 7.2))

TU_CDR1 <- mutate(TU_CDR1, Hydrophob.TU = hydrophobicity(TU_CDR1$TU))

TU_CDR1 <- mutate(TU_CDR1, Charge.TU = charge(TU_CDR1$TU, pH = 7.2))

ggplot(Ab_s_CDR1, aes(Hydrophob.Ab)) + geom_histogram(bins = 68, colour = "black", fill = "lightblue") + 
  geom_vline(aes(xintercept = mean(Hydrophob.Ab)), size = 1, colour = "black",linetype = "dashed") +
  scale_x_continuous(breaks =c(-2, -1.5, -1.0, -0.5, 0, 0.5) ) +
  labs(title = "Ab Specific Nanobody CDR1 Hydrophophobicities", x = "Hydrophobicity", y = "Count")

ggplot(Ab_s_CDR1, aes(Charge.Ab)) + geom_histogram(bins = 10, colour = "black", fill = "pink") + 
  geom_vline(aes(xintercept = mean(Charge.Ab)), size = 1, colour = "black",linetype = "dashed") +
  labs(title = "Ab Specific Nanobody CDR1 Charges", x = "Charge", y = "Count")

ggplot(TU_CDR1, aes(Hydrophob.TU)) + geom_histogram(bins = 68, colour = "black", fill = "lightblue") + 
  geom_vline(aes(xintercept = mean(Hydrophob.TU)), size = 1, colour = "black",linetype = "dashed") +
  scale_x_continuous(breaks =c(-2, -1.5, -1.0, -0.5, 0, 0.5) ) +
  labs(title = "TU Nanobody CDR1 Hydrophophobicities", x = "Hydrophobicity", y ="Count")

ggplot(TU_CDR1, aes(Charge.TU)) + geom_histogram(bins = 10, colour = "black", fill = "pink") + 
  geom_vline(aes(xintercept = mean(Charge.TU)), size = 1, colour = "black",linetype = "dashed") +
  labs(title = "TU Nanobody CDR1 Charges", x = "Charge", y = "Count")


#CDR2
Ab_s_CDR2 <- mutate(Ab_s_CDR2, Hydrophob.Ab = hydrophobicity(Ab_s_CDR2$Ab))

Ab_s_CDR2 <- mutate(Ab_s_CDR2, Charge.Ab = charge(Ab_s_CDR2$Ab, pH = 7.2))

TU_CDR2 <- mutate(TU_CDR2, Hydrophob.TU = hydrophobicity(TU_CDR2$TU))

TU_CDR2 <- mutate(TU_CDR2, Charge.TU = charge(TU_CDR2$TU, pH = 7.2))

ggplot(Ab_s_CDR2, aes(Hydrophob.Ab)) + geom_histogram(bins = 68, colour = "black", fill = "lightblue") + 
  geom_vline(aes(xintercept = mean(Hydrophob.Ab)), size = 1, colour = "black",linetype = "dashed") +
  scale_x_continuous(breaks =c(-2, -1.5, -1.0, -0.5, 0, 0.5) ) +
  labs(title = "Ab Specific Nanobody CDR2 Hydrophophobicities", x = "Hydrophobicity", y = "Count")

ggplot(Ab_s_CDR2, aes(Charge.Ab)) + geom_histogram(bins = 10, colour = "black", fill = "pink") + 
  geom_vline(aes(xintercept = mean(Charge.Ab)), size = 1, colour = "black",linetype = "dashed") +
  labs(title = "Ab Specific Nanobody CDR2 Charges", x = "Charge", y = "Count")

ggplot(TU_CDR2, aes(Hydrophob.TU)) + geom_histogram(bins = 68, colour = "black", fill = "lightblue") + 
  geom_vline(aes(xintercept = mean(Hydrophob.TU)), size = 1, colour = "black",linetype = "dashed") +
  scale_x_continuous(breaks =c(-2, -1.5, -1.0, -0.5, 0, 0.5) ) +
  labs(title = "TU Nanobody CDR2 Hydrophophobicities", x = "Hydrophobicity", y ="Count")

ggplot(TU_CDR2, aes(Charge.TU)) + geom_histogram(bins = 10, colour = "black", fill = "pink") + 
  geom_vline(aes(xintercept = mean(Charge.TU)), size = 1, colour = "black",linetype = "dashed") +
  labs(title = "TU Nanobody CDR2 Charges", x = "Charge", y = "Count")


#CDR3
Ab_s_CDR3 <- mutate(Ab_s_CDR3, Hydrophob.Ab = hydrophobicity(Ab_s_CDR3$Ab))

Ab_s_CDR3 <- mutate(Ab_s_CDR3, Charge.Ab = charge(Ab_s_CDR3$Ab, pH = 7.2))

TU_CDR3 <- mutate(TU_CDR3, Hydrophob.TU = hydrophobicity(TU_CDR3$TU))

TU_CDR3 <- mutate(TU_CDR3, Charge.TU = charge(TU_CDR3$TU, pH = 7.2))

ggplot(Ab_s_CDR3, aes(Hydrophob.Ab)) + geom_histogram(bins = 68, colour = "black", fill = "lightblue") + 
  geom_vline(aes(xintercept = mean(Hydrophob.Ab)), size = 1, colour = "black",linetype = "dashed") +
  scale_x_continuous(breaks =c(-2, -1.5, -1.0, -0.5, 0, 0.5) ) +
  labs(title = "Ab Specific Nanobody CDR3 Hydrophophobicities", x = "Hydrophobicity", y = "Count")

ggplot(Ab_s_CDR3, aes(Charge.Ab)) + geom_histogram(bins = 10, colour = "black", fill = "pink") + 
  geom_vline(aes(xintercept = mean(Charge.Ab)), size = 1, colour = "black",linetype = "dashed") +
  labs(title = "Ab Specific Nanobody CDR3 Charges", x = "Charge", y = "Count")

ggplot(TU_CDR3, aes(Hydrophob.TU)) + geom_histogram(bins = 68, colour = "black", fill = "lightblue") + 
  geom_vline(aes(xintercept = mean(Hydrophob.TU)), size = 1, colour = "black",linetype = "dashed") +
  scale_x_continuous(breaks =c(-2, -1.5, -1.0, -0.5, 0, 0.5) ) +
  labs(title = "TU Nanobody CDR3 Hydrophophobicities", x = "Hydrophobicity", y ="Count")

ggplot(TU_CDR3, aes(Charge.TU)) + geom_histogram(bins = 10, colour = "black", fill = "pink") + 
  geom_vline(aes(xintercept = mean(Charge.TU)), size = 1, colour = "black",linetype = "dashed") +
  labs(title = "TU Nanobody CDR3 Charges", x = "Charge", y = "Count")

#determine mean and sds

mean(Ab_s_CDR1$Hydrophob.Ab)
mean(Ab_s_CDR1$Charge.Ab)
mean(Ab_s_CDR2$Hydrophob.Ab)
mean(Ab_s_CDR2$Charge.Ab)
mean(Ab_s_CDR3$Hydrophob.Ab)
mean(Ab_s_CDR3$Charge.Ab)
mean(Ab_s_combined_CDRs$Hydrophob.Ab)
mean(Ab_s_combined_CDRs$Charge.Ab)

mean(TU_CDR1$Hydrophob.Ab)
mean(TU_CDR1$Charge.Ab)
mean(TU_CDR2$Hydrophob.Ab)
mean(TU_CDR2$Charge.Ab)
mean(TU_CDR3$Hydrophob.Ab)
mean(TU_CDR3$Charge.Ab)
mean(TU_combined_CDRs$Hydrophob.Ab)
mean(TU_combined_CDRs$Charge.Ab)


sd(Ab_s_CDR1$Hydrophob.Ab)
sd(Ab_s_CDR1$Charge.Ab)
sd(Ab_s_CDR2$Hydrophob.Ab)
sd(Ab_s_CDR2$Charge.Ab)
sd(Ab_s_CDR3$Hydrophob.Ab)
sd(Ab_s_CDR3$Charge.Ab)
sd(Ab_s_combined_CDRs$Hydrophob.Ab)
sd(Ab_s_combined_CDRs$Charge.Ab)

sd(TU_CDR1$Hydrophob.Ab)
sd(TU_CDR1$Charge.Ab)
sd(TU_CDR2$Hydrophob.Ab)
sd(TU_CDR2$Charge.Ab)
sd(TU_CDR3$Hydrophob.Ab)
sd(TU_CDR3$Charge.Ab)
sd(TU_combined_CDRs$Hydrophob.Ab)
sd(TU_combined_CDRs$Charge.Ab)


t.test(Ab_s_CDR3$Charge.Ab, TU_CDR3$Charge.TU)
t.test(Ab_s_combined_CDRs$Charge.Ab, TU_combined_CDRs$Charge.TU)
