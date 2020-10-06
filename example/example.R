library(CADEG)
library(VennDiagram)
library(Biobase)
library(limma)
library(gplots)

# step01: read 'Brapa_exp.txt', and removed the unexpressed genes
RemoveUnexpGenes(input_file="Brapa_exp.txt", cutoff=0.5, output_file="Brapa_exp.01.RemUnexpGene.txt")

# step02: calculate values of P-value, PDR and FoldChange for the comparison between two samples of SL and SP
CalPvalueLog(input_file="Brapa_exp.01.RemUnexpGene.txt", sampleA="SL", sampleB="SP", output_file="Brapa_exp.02.SL-vs-SP.txt")

# step02: calculate values of P-value, PDR and FoldChange for the comparison between two samples of SL and XL
CalPvalueLog(input_file="Brapa_exp.01.RemUnexpGene.txt", sampleA="SL", sampleB="XL", output_file="Brapa_exp.02.SL-vs-XL.txt")

# step02: calculate values of P-value, PDR and FoldChange for the comparison between two samples of SL and XP
CalPvalueLog(input_file="Brapa_exp.01.RemUnexpGene.txt", sampleA="SL", sampleB="XP", output_file="Brapa_exp.02.SL-vs-XP.txt")

# step02: select differentailly expressed genes between two samples of SL and SP, based on P-value/PDR and FoldChange values
SelDiffTwoCriterion(input_file="Brapa_exp.02.SL-vs-SP.txt", valueA=0.01, typeA="Pvalue", valueB=2, "FoldChange", "Brapa_exp.02.SL-vs-SP.Pvalue0.01-FoldChange2.txt")

# step02: select differentailly expressed genes between two samples of SL and XL, based on P-value/PDR and FoldChange values
SelDiffTwoCriterion(input_file="Brapa_exp.02.SL-vs-XL.txt", valueA=0.01, typeA="Pvalue", valueB=2, typeB="FoldChange", output_file="Brapa_exp.02.SL-vs-XL.Pvalue0.01-FoldChange2.txt")

# step02: select differentailly expressed genes between two samples of SL and XP, based on P-value/PDR and FoldChange values
SelDiffTwoCriterion(input_file="Brapa_exp.02.SL-vs-XP.txt", valueA=0.01, typeA="Pvalue", valueB=2, typeB="FoldChange", output_file="Brapa_exp.02.SL-vs-XP.Pvalue0.01-FoldChange2.txt")

# step03: draw venn graph of DEGs
VennThreeSet(file_sampleA="Brapa_exp.02.SL-vs-SP.Pvalue0.01-FoldChange2.txt", file_sampleB="Brapa_exp.02.SL-vs-XL.Pvalue0.01-FoldChange2.txt", file_sampleC="Brapa_exp.02.SL-vs-XP.Pvalue0.01-FoldChange2.txt", sampleA="SL-vs-SP", sampleB="SL-vs-XL", sampleC="SL-vs-XP", header=TRUE, output_file = "Brapa_exp.03.venn.SL-vs-SP.XL.XP.pdf", output_geneID=TRUE)

# step04: select express information for genes in 'Brapa_exp.03.venn.SL-vs-SP.XL.XP.pdf.S1.S2.S3.overlapID.txt'
SelGeneAnnotation(input_file="Brapa_exp.01.RemUnexpGene.txt", geneID_file="Brapa_exp.03.venn.SL-vs-SP.XL.XP.pdf.S1.S2.S3.overlapID.txt", output_file="Brapa_exp.03.venn.SL-vs-SP.XL.XP.pdf.S1.S2.S3.overlapID.exp.txt")

# step05: draw heatmap for genes
DrawHeatmap(input_file="Brapa_exp.03.venn.SL-vs-SP.XL.XP.pdf.S1.S2.S3.overlapID.exp.txt", to_cluster_sample=TRUE, to_cluster_gene=TRUE, output_file="Brapa_exp.05.SL-vs-SP.XL.XP.pdf.S1.S2.S3.overlapID.exp.heatmap.pdf")

# step06: add Pfam annotation to file 'Brapa_exp.03.venn.SL-vs-SP.XL.XP.pdf.S1.S2.S3.overlapID.exp.txt'
PlusPfamAnnotation(input_file="Brapa_exp.03.venn.SL-vs-SP.XL.XP.pdf.S1.S2.S3.overlapID.exp.txt", Pfam_domain_file="Brapa_Pfam.out.1e-05.stat.genedomain", evalue=1e-5, output_file="Brapa_exp.06.S1.S2.S3.overlapID.exp.Pfam.txt")

# step07: identify enrich expressed genes
# the result by the cutoff of 'P-value < 0.05' can be found in 'Brapa_exp.07.S1.S2.S3.overlapID.exp.Pfam.dis.txt.Pvalue0.05.txt'
DEGPfamAnnotationDis(DEG_Pfam_file="Brapa_exp.06.S1.S2.S3.overlapID.exp.Pfam.txt", exp_gene_file="Brapa_exp.01.RemUnexpGene.txt", Pfam_domain_file="Brapa_Pfam.out.1e-05.stat.genedomain", evalue=1e-5, output_file="Brapa_exp.07.S1.S2.S3.overlapID.exp.Pfam.dis.txt")
