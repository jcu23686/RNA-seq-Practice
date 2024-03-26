#Setting up working directory and loading R libraries

	library(ballgown); library(genefilter); library(dplyr); library(devtools); library(ggplot2); library(cowplot)

#loading phenotypes and RNA seq FPKM data generated by Stringtie

	pheno_data <- read.csv("phenodata.csv")
	bg <- ballgown(dataDir = "ballgown", samplePattern = "[D|E|F]", pData = pheno_data)
	bg_filt <- subset(bg, "rowVars(gexpr(bg)) >1", genomesubset=TRUE)
	results_genes <- stattest(bg_filt,feature="gene",covariate="time",getFC=TRUE, meas="FPKM")

#Genewise analysis
	results_genes$mean <- rowMeans(gexpr(bg_filt))
	results_genes_complete <- results_genes[complete.cases(results_genes), ]

# Make box plots for expression levels for all 14 libraries
	
	short_names=c("Noon_D08", "Night_D12", "Noon_E02", "Night_E05", "Noon_E06", "Night_E07", "Noon_E08", "Night_F01", "Noon_F03", "Night_F04", "Night_F05", "Noon_F06", "Night_F07", "Night_F09")
	fpkm = texpr(bg_filt,meas="FPKM")
	fpkm = log2(fpkm+1)
	colnames(fpkm) <- short_names
	outfile="NoonNightLib_boxblot.pdf"
	pdf(file=outfile)
	boxplot(fpkm,col=as.numeric(pheno_data$time)+1,names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 14 libraries")
	dev.off()

# MDS Plots (similar to PCA)

	sum = apply(fpkm[,], 1, sum)
	i = which(sum[] > 2)
	r=cor(fpkm[i,], use="pairwise.complete.obs", method="pearson")
	d=1-r
	mds=cmdscale(d, k=2, eig=TRUE)
	par(mfrow=c(1,1))
	outfile="Noon_Night_MDS.pdf"
	pdf(file=outfile)
	plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes) for all libraries", xlim=c(-0.25,0.25), ylim=c(-0.05,0.05))
	points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
	text(mds$points[,1], mds$points[,2], short_names)
	dev.off()

#DEG – differential expression analysis 

	indices <- match(results_genes_complete$id, texpr(bg_filt, 'all')$gene_id)
	gene_name_for_results <- texpr(bg_filt, 'all')$t_name[indices]	
	results_genes_complete <- data.frame(geneNames=gene_name_for_results, results_genes_complete)
	sig<-results_genes_complete %>% filter(pval < 0.05)
	write.csv(sig, "Noon_Night_DEGs_pval.csv")

#MA plot (Log-Fold Change versus Log-Concentration) highlighting points with p-value < 0.05

d<-results_genes_complete[!(results_genes_complete$geneNames=="."),]

ggplot(results_genes_complete, aes(log2(mean), log2(fc), colour = pval<0.05)) + scale_color_manual(values=c("#999999", "#aa0a3c")) + geom_point(aes(alpha=pval<0.05)) + scale_alpha_manual(guide="none", values = c(0.1, 1)) + geom_hline(yintercept=0) +  xlim(-10, 20) + ylim(-8, 8) + theme_bw()

ggsave("MA_plot.noon_vs_night.pdf", height=5, width=7, useDingbats=FALSE)

#MA plot (Log-Fold Change versus Log-Concentration) highlighting points with q-value < 0.05 

ggplot(results_genes_complete, aes(log2(mean), log2(fc), colour = qval<0.05)) + scale_color_manual(values=c("#999999", "#aa0a3c")) + geom_point(aes(alpha=qval<0.05)) + scale_alpha_manual(guide="none", values = c(0.1, 1)) + geom_hline(yintercept=0) + xlim(-10, 20) + ylim(-8, 8) + theme_bw()

ggsave("MA_qvalue_plot.noon_vs_night.pdf", height=5, width=7, useDingbats=FALSE)

# Quit R

	quit ()

