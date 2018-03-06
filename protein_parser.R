source("https://bioconductor.org/biocLite.R")
install.packages("plyr")
library("seqinr")
library("plyr")
install.packages("stringr")
library(stringr)


x <- seqinr::read.fasta("~/biohack/plasm.fasta", seqtype = "AA", as.string = TRUE)

f <- as.character(x[1])
m <- x[1]
prot_name <- as.character(attributes(m))
z <- rep("N", 5)
gregexpr(z,f, useBytes = TRUE)
coll <- paste(z, collapse = "")
length(table(strsplit(f, coll)))


