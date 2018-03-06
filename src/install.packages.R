f <- package_installed()
if (f == TRUE){
		source("https://bioconductor.org/biocLite.R")
		install.packages("gplots")
		install.packages("ggplot2")
		install.packages("pheatmap")
		biocLite("keggorthology")
		biocLite("PFAM.db")
		biocLite("dcGOR")
} else {
		print("All R packages installed!")
}


