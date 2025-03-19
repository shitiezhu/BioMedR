
library(tidyverse)
# seqlist = read.csv("pGIPZshRNA_sequence list_mSpry12.csv")

pGIPZ <- function(shList = shList,species= "human"){

  shList <- shList %>%
  mutate(ensembl_transcript_id = gsub( "\\.\\d*", "", ensembl_transcript_id))

  if(species == "human"){
    ensembl = biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
    seq <- biomaRt::getSequence(id = shList$ensembl_transcript_id, type = "ensembl_transcript_id",
                            seqType = "cdna", mart = ensembl )
  }else if(species == "mouse"){
    ensembl = biomaRt::useMart("ensembl", dataset="mmusculus_gene_ensembl")
    seq <- biomaRt::getSequence(id = shList$ensembl_transcript_id, type = "ensembl_transcript_id",
                                seqType = "cdna", mart = ensembl )
  }

seq2 <- seq %>%
  full_join(shList, by = "ensembl_transcript_id") %>%
  rowwise() %>%
  mutate(shSense = as.character(Biostrings::reverseComplement(Biostrings::DNAString(shAntisense)))) %>%
  mutate(shSense_22 = str_extract(cdna, paste0("..", shSense, "."))) %>%
  na.omit() %>%
  rowwise() %>%
  mutate(shAntisense_22 = as.character(Biostrings::reverseComplement(Biostrings::DNAString(shSense_22)))) %>%
  mutate(FORWARD = ifelse(str_sub(shSense_22, 1,1) == "A" | str_sub(shSense_22, 1,1) == "T", sub("^.", "c", shSense_22),
                          ifelse(str_sub(shSense_22, 1,1) == "C" | str_sub(shSense_22, 1,1) == "G", sub("^.", "a", shSense_22), shSense_22)),
         REVERSE = ifelse(str_sub(shSense_22, 1,1) == "A", sub(".$", "t", shAntisense_22),
                          ifelse(str_sub(shSense_22, 1,1) == "G", sub(".$", "c", shAntisense_22),
                                 ifelse(str_sub(shSense_22, 1,1) == "T", sub(".$", "a", shAntisense_22),
                                        ifelse(str_sub(shSense_22, 1,1) == "C", sub(".$", "g", shAntisense_22),shAntisense_22)
         )))) %>%
  rowwise() %>%
  mutate(shTarget = paste0("TGCTGTTGACAGTGAGCG", FORWARD, "TAGTGAAGCCACAGATGTA", REVERSE, "TGCCTACTGCCTCGGA"),
         shTarget_F = str_sub(shTarget, 1,59),
         shTarget_R = str_sub(as.character(Biostrings::reverseComplement(Biostrings::DNAString(shTarget))), 1,59)) |>
  select(Genes,shSense, shTarget, shTarget_F, shTarget_R) |>
  pivot_longer(cols = c("shTarget_F", "shTarget_R"), names_to = "ID", values_to = "sequence") |>
  rowwise() |>
  mutate(ID = gsub("Target", Genes, ID))

write.csv(seq2,"pGIPZshRNA_sequence_ready.csv", row.names = F)
}

# library(tidyverse)
# pGIPZ(shList = seqlist, species = "mouse")


# str_extract(seq2[1,1], "..CTGTCGTTCACCTTCAACT.")
#
# BiocManager::install("Biostrings")
#
# cc <- Biostrings::reverseComplement(Biostrings::DNAString("TCAACTTCCACTTGCTGTC")) %>% as.character()
# # stringi::stri_reverse("ABC")
#
# seq[1,1]
# #
# # grep("CCAGTTACTCTAG", seq[1,1])
# # ?grep("CCAGTTACTCTAG", seq[1,1], value = T)
# #
# # s <- "PRODUCT colgate good but not goodOKAY"
# # sub(".*PRODUCT *(.*?) *OKAY.*", "\\1", s)
# #
# # x <- "helloxxxotherstuff"
# # sub("xxx.*", "", x)
# # sub(".*xxx", "bbb", x)
# # x <- c("apple", "banana", "pear")
# stringr::str_extract(seq[1,1], "..TACGGAAGAGCCGCCACGGG.")
