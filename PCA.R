source("https://bioconductor.org/biocLite.R")
biocLite("devtools")
install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot")

g <- grep("pfam", list.files("~/biohack/pfam/N15/"), value = TRUE)
setwd("~/biohack/pfam/N15/")
df_N15 <- data.frame()
for (f in g){
  z <- read.csv(file = f, header = FALSE)
  names(z) <- c("num", "header", "length", "family", "score")
  z$sample <- paste(f)
  df_N15 <- rbind(z, df_N15)
}

g <- grep("pfam", list.files("~/biohack/pfam/N10/"), value = TRUE)
setwd("~/biohack/pfam/N10/")
df_N10 <- data.frame()
for (f in g){
  z <- read.csv(file = f, header = FALSE)
  names(z) <- c("num", "header", "length", "family", "score")
  z$sample <- paste(f)
  df_N10 <- rbind(z, df_N10)
}


g <- grep("pfam", list.files("~/biohack/pfam/"), value = TRUE)
setwd("~/biohack/pfam/")
df_N0 <- data.frame()
for (f in g){
  z <- read.csv(file = f, header = FALSE)
  names(z) <- c("num", "header", "length", "family", "score")
  z$sample <- paste(f)
  df_N0 <- rbind(z, df_N10)
}



