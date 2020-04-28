require(data.table)
require(seqinr)
require(stringi)
require(stringr)
require(dplyr)
require(tidyr)

#grep the reads with correct nanobody frame. Used regular expression here to allow mutations that don't change aa sequence
#read in file
temp1 <- fread("Input_tally.tsv",header=FALSE)

#add headers
colnames(temp1)<-c("nucleotide","INP1")

#EXACT NUCLEOTIDE
#temp1<-temp1[grepl("A(?:CGT|CGC|CGA|CGG|AGA|ACG)(?:TTA|TTG|CTT|CTC|CTA|CTG)(?:TCT|TCC|TCA|TCG|AGT|AGC)TG[T|C]GC.GC.(?:TCT|TCC|TCA|TCG|AGT|AGC)GG..{21}ATGGGGTGGTTTCGCCAGGCACCTGGCAAAGAACGTGAATTTGTTGCAGCAATTAGT.{21}TACTACGCAGATTCCGTTAAGGGACGCTTCACAATTTCGCGCGACAATGCAAAAAATACCGTGTATTTACAAATGAATTCGTTGAAGCCGGAAGACACTGCGACTTATTATTGTGCG.{27}TATTGG",nucleotide)]temp1<-temp1[grepl(".{21}ATGGG.TGGTT(T|C).{24}GA[A|G]TT[C|T]GT.GC.GC.AT.{25}TA",nucleotide)]

#ALLOW SYNOMYNOUS NUCLEOTIDE SUBSITUTIONS
temp1<-temp1[grepl("A(?:CGT|CGC|CGA|CGG|AGA|ACG)(?:TTA|TTG|CTT|CTC|CTA|CTG)(?:TCT|TCC|TCA|TCG|AGT|AGC)TG[T|C]GC.GC.(?:TCT|TCC|TCA|TCG|AGT|AGC)GG..{21}ATGGG.TGGTT[T|C](?:CGT|CGC|CGA|CGG|AGA|ACG)CA[A|G]GC.CC.GG.AA[A|G]GA[A|G](?:CGT|CGC|CGA|CGG|AGA|ACG)GA[A|G]TT[T|C]GT.GC.GC.AT[T|C|A](?:TCT|TCC|TCA|TCG|AGT|AGC).{21}TA[T|C]TA[T|C]GC.GA[T|C](?:TCT|TCC|TCA|TCG|AGT|AGC)GT.AA[A|G]GG.(?:CGT|CGC|CGA|CGG|AGA|ACG)TT[T|C]AC.AT[T|C|A](?:TCT|TCC|TCA|TCG|AGT|AGC)(?:CGT|CGC|CGA|CGG|AGA|ACG)GA[T|C]AA[T|C]GC.AA[A|G]AA[T|C]AC.GT.TA[T|C](?:TTA|TTG|CTT|CTC|CTA|CTG)CA[A|G]ATGAA[T|C](?:TCT|TCC|TCA|TCG|AGT|AGC)(?:TTA|TTG|CTT|CTC|CTA|CTG)AA[A|G]CC.GA[A|G]GA[T|C]AC.GC.AC.TA[T|C]TA[T|C]TG[T|C]GC..{27}TA[T|C]TGG",nucleotide)]

#REMOVE FIRST NUCLEOTIDE TO PUT SEQUENCE IN FRAME
temp1 <- temp1[, nuc:=gsub("^A","",nucleotide)]

#temp1<-temp1[grepl(".{27}.GC[A|G]CA[A|G]TA[A|G]TA.GT.GC",nucleotide)]

#grep CDR1, 2, 3 by trimming off the upstream and downstream sequences
temp1[,CDR1:=gsub("ATGGG.TGGTT(T|C).{24}GA[A|G]TT[C|T]GT.GC.GC.AT.{4}","_",nucleotide)
      ][,CDR2:=str_extract(CDR1,"_.{21}")
        ][,CDR2:=gsub("_","", CDR2)
          ][,CDR3:=gsub(".GC[A|G]CA[A|G]TA[A|G]TA.GT.GC","8",nucleotide)
            ][,CDR3:=str_extract(CDR3,".{27}8")
              ][,CDR3:=gsub("8","",CDR3)
                ][,CDR1:=str_extract(CDR1,".{21}_")]

#REMOVE NUCLEOTIDE COLUMN
temp1 <- select(temp1, -c(nucleotide))

```

b. translate CDR1, 2, 3, into peptide1, 2, 3
```{r}
require(data.table)
require(seqinr)
require(stringi)
require(stringr)

#TRANSLATE CDRs

temp1[,peptide1:=paste(translate(s2c(CDR1),
                                  frame=0),sep="",collapse=""),by=CDR1
       ][,peptide2:=paste(translate(s2c(CDR2),
                                    frame=0),sep="",collapse=""),by=CDR2
         ][,peptide3:=paste(translate(s2c(CDR3),
                                      frame=0, sens="R"),sep="",collapse=""),by=CDR3]

#MAKE TABLE WITH OUT NA IN PEPTIDE3 COLUMN
temp1 <- na.omit(temp1, cols ="peptide3")

#MAKE A COMBINED CDR COLUMN
temp1 <- unite(temp1, combined, peptide1:peptide3, sep ="", remove = FALSE)

fwrite(temp1,"Input_nucleotides_and_peptides.csv")

#MAKE DATA WITH ONLY UNIQUE COMBINED CDRs
unique <- temp1[!duplicated(temp1$combined), ]


fwrite(unique, "Unique_Combined_CDRs_Input_nucleotides_and_peptides.csv")

#Isolate CDR combinations and CDR1-3 as individual data
combined <- select(unique, combined)
CDR1 <- select(unique, peptide1)
CDR2 <- select(unique, peptide2)
CDR3 <- select(unique, peptide3)

#Write to txt file


write.table(combined, file = "Unique_Input_Combined_CDRs.txt", row.names = FALSE, col.names = FALSE)
write.table(CDR1, file = "Unique_Input_CDR1.txt", row.names = FALSE, col.names = FALSE)
write.table(CDR2, file = "Unique_Input_CDR2.txt", row.names = FALSE, col.names = FALSE)
write.table(CDR3, file = "Unique_Input_CDR3.txt", row.names = FALSE, col.names = FALSE)

