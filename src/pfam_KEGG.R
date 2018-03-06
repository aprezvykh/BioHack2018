
#library("keggorthology")
library("PFAM.db")
library("dcGOR")
library(ggplot2)
current.directory <- getwd()
current.directory <- "~/prFinder/test_data/"
setwd(current.directory)

z <- grep("pfam", list.files(current.directory), value = TRUE)

if (length(z) == 0){
  print("No input file specified!")
  break
}


for (f in z){
      #read and convert uproc output
      a <- read.csv(f, header = FALSE)
      names(a) <- c("number", "uniprot", 
                    "Length", "Family",
                    "Score")
      domains <- a[which(a$Score > 2),]
      pfam.families <- as.character(domains$Family)
      Pfam <- dcRDataLoader('Pfam')
      #calculating enrichment for all ontologies
      enrich.output.BP <- dcEnrichment(pfam.families, domain="Pfam", ontology="GOBP")
      enrich.output.viewed.BP <- view(enrich.output.BP, top_num=30, sortBy="pvalue", details=TRUE)
      enrich.output.MF <- dcEnrichment(data, domain="Pfam", ontology="GOMF")
      enrich.output.viewed.MF <- view(enrich.output.MF, top_num=30, sortBy="pvalue", details=TRUE)
      enrich.output.CC <- dcEnrichment(data, domain="Pfam", ontology="GOCC")
      enrich.output.viewed.CC <- view(enrich.output.CC, top_num=30, sortBy="pvalue", details=TRUE)
      #finding top motif enrichment
      motifs.BP <- as.character(view(enrich.output.BP, top_num=1, sortBy="pvalue", details=T)$members)
      motifs.BP <- unlist(strsplit(motifs.BP,","))
      top.motif.BP <- Data(Pfam)[match(motifs.BP,rowNames(Pfam)),]
      motifs.MF <- as.character(view(enrich.output.MF, top_num=1, sortBy="pvalue", details=T)$members)
      motifs.MF <- unlist(strsplit(motifs.MF,","))
      top.motif.MF <- Data(Pfam)[match(motifs.MF,rowNames(Pfam)),]
      motifs.CC <- as.character(view(enrich.output.CC, top_num=1, sortBy="pvalue", details=T)$members)
      motifs.CC <- unlist(strsplit(motifs.CC,","))
      top.motif.CC <- Data(Pfam)[match(motifs.CC,rowNames(Pfam)),]
      #writing output results
      write.csv(enrich.output.viewed.BP, paste(f,"BP_out.csv", sep = "_"))
      write.csv(enrich.output.viewed.MF, paste(f,"MF_out.csv", sep = "_"))
      write.csv(enrich.output.viewed.CC, paste(f,"CC_out.csv", sep = "_"))
      write.csv(top.motif.BP, paste(f, "top_motif_BP.csv", sep = "_"))
      write.csv(top.motif.MF, paste(f, "top_motif_MF.csv", sep = "_"))
      write.csv(top.motif.CC, paste(f, "top_motif_CC.csv", sep = "_"))
      #visualizing barplots
      #BP
      barplot.BP.subset <- enrich.output.viewed.BP[1:10,]
      g <- ggplot(barplot.BP.subset) + 
              geom_bar(aes(x = term_name, y = nAnno, 
              fill = adjp), stat = "identity") + 
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(filename = "barplot_BP.png", plot = g, units = "cm", height = 15, width = 15)
      #MF
      barplot.MF.subset <- enrich.output.viewed.MF[1:10,]
      g <- ggplot(barplot.MF.subset) + 
              geom_bar(aes(x = term_name, y = nAnno, 
              fill = adjp), stat = "identity") + 
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(filename = "barplot_MF.png", plot = g, units = "cm", height = 15, width = 15)
      
      #CC
      barplot.CC.subset <- enrich.output.viewed.CC[1:10,]
      g <- ggplot(barplot.CC.subset) + 
              geom_bar(aes(x = term_name, y = nAnno, 
              fill = adjp), stat = "identity") + 
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(filename = "barplot_CC.png", plot = g, units = "cm", height = 15, width = 15)
      
      ###Graph visualization
      ###FIX SAVING
      visEnrichment(enrich.output.BP, num_top_nodes=5, layout.orientation="top_bottom", zlim=c(0,4))
      visEnrichment(enrich.output.MF, num_top_nodes=5, layout.orientation="top_bottom", zlim=c(0,4))
      visEnrichment(enrich.output.CC, num_top_nodes=5, layout.orientation="top_bottom", zlim=c(0,4))
      
      ### PROTEIN INTERACTIONS
      #g <- dcRDataLoader('onto.GOBP')
      #Anno <- dcRDataLoader('Pfam2GOBP')
      #dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths",verbose=T)
      #dnetwork <- dcDAGdomainSim(g=dag, domains=as.character(a$Family), 
      #                           method.domain="BM.average", method.term="Resnik", 
      #                           parallel=FALSE, verbose=TRUE, force = TRUE)

}


