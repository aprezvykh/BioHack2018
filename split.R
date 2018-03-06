install.packages("dplyr")
library("dplyr")
setwd("~/biohack/FINAL/3d7/")
N5.pfam <- read.csv("3d7_filteredNQ5.fasta.pfam", header = FALSE)
N5.am <- read.table("3d7_filteredNQ5.fasta.am.csv", header = TRUE, sep = "\t")
N5.am$SEQid <- sub("<unknown description>", "", N5.am$SEQid)
N5.am$SEQid <- sub(" ", "", N5.am$SEQid)
N5.pfam$V2 <- as.character(N5.pfam$V2)

i <- intersect(N5.pfam$V2, N5.am$SEQid)

u <- unique(N5.pfam$V2)
N5.pfam$V1 <- NULL
N5.pfam$V4 <- NULL
N5.pfam$V3 <- NULL
for (f in u){
  g <- grep(f, N5.pfam$V2)
  d <- N5.pfam[g,]
  if(nrow(d)>1){
    cs <- colSums(d[,2])
    print(cs)
  }
}

