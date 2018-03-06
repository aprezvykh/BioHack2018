a <- read.csv("~/test.run.csv", header = FALSE)

x <- data.frame(table(unlist(a$V4)))
