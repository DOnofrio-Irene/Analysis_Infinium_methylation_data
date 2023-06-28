rm(list=ls())
setwd('C:/Users/Irene/Desktop/LM Bioinformatics/DRD/Module2/Final Report-20230616/')
library(minfi)

baseDir <- ("C:/Users/Irene/Desktop/LM Bioinformatics/DRD/Module2/Final Report-20230616/Input/")
targets <- read.metharray.sheet(baseDir)
RGset <- read.metharray.exp(targets = targets) 
save(RGset,file="RGset.RData")



Red <- data.frame(getRed(RGset))
dim(Red)
head(Red)

Green <- data.frame(getGreen(RGset))
dim(Green)
head(Green)


red_probe <- Red[rownames(Red) == "64638362", ]
green_probe <- Green[rownames(Green)=="64638362",]

load("C:/Users/Irene/Desktop/LM Bioinformatics/DRD/Module2/DRD_2023_Lesson_2-20230516/Illumina450Manifest_clean.RData")
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="64638362", 'Infinium_Design_Type']


df_address <- data.frame(Sample=colnames(green_probe), Red_fluor=unlist(red_probe, use.names = FALSE ), Green_fluor=unlist(green_probe, use.names = FALSE ), Type = "II")
df_address




MSet.raw <- preprocessRaw(RGset)
MSet.raw



qc <- getQC(MSet.raw)
plotQC(qc)

controlStripPlot(RGset, controls="NEGATIVE")

detP <- detectionP(RGset)
failed <- detP>0.05


num_failed = colSums(failed)


df_failed <- data.frame(Sample=colnames(green_probe), Failed_Position=num_failed, perc_failed_probes=(means_of_columns <- colMeans(failed))*100, row.names = NULL)

df_failed

wt <- targets[targets$Group=="WT", "Basename"]
mut <- targets[targets$Group=="MUT", "Basename"]
wt <- gsub(baseDir, "", wt)
mut <- gsub(baseDir, "", mut)

wtSet <- MSet.raw[,colnames(MSet.raw) %in% wt]
mutSet <- MSet.raw[,colnames(MSet.raw) %in% mut]

wtBeta <- getBeta(wtSet)
wtM <- getM(wtSet)
mutBeta <- getBeta(mutSet)
mutM <- getM(mutSet)

mean_wtBeta <- apply(wtBeta,MARGIN=1,mean,na.rm=T)
mean_mutBeta <- apply(mutBeta,MARGIN=1,mean,na.rm=T)
mean_wtM <- apply(wtM,MARGIN=1,mean,na.rm=T)
mean_mutM <- apply(mutM,MARGIN=1,mean,na.rm=T)


d_mean_wtBeta <- density(mean_wtBeta,na.rm=T)
d_mean_mutBeta <- density(mean_mutBeta,na.rm=T)
d_mean_wtM <- density(mean_wtM,na.rm=T)
d_mean_mutM <- density(mean_mutM,na.rm=T)


par(mfrow=c(1,2))
plot(d_mean_wtBeta,main="Density of Beta Values",col="#0067E6", lwd=2.5)
lines(d_mean_mutBeta,main="Density of Beta Values",col="#E50068", lwd=2.5)

legend('topright', legend=c("WT", "MUT"), 
       fill = c("#0067E6","#E50068"), cex=1
)
plot(d_mean_wtM,main="Density of M Values",col="#0067E6", lwd=2.5)
lines(d_mean_mutM,main="Density of M Values",col="#E50068", lwd=2.5)
legend('topright', legend=c("WT", "MUT"), 
       fill = c("#0067E6","#E50068"), cex=1
)

#creation of two dataframes, containing only type I (dfI) or type II (dfII) probes
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)

#subsetting of the beta matrix according to the chemistry of the probes
beta <- getBeta(MSet.raw)
beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID,]

#calculation of the mean, sd and their corresponding density for the raw data
mean_of_beta_I <- apply(beta_I,1,mean,na.rm=T)
mean_of_beta_II <- apply(beta_II,1,mean,na.rm=T)
d_mean_of_beta_I <- density(mean_of_beta_I,na.rm=T)
d_mean_of_beta_II <- density(mean_of_beta_II,na.rm=T)

sd_of_beta_I <- apply(beta_I,1,sd,na.rm=T)
sd_of_beta_II <- apply(beta_II,1,sd,na.rm=T)
d_sd_of_beta_I <- density(sd_of_beta_I,na.rm=T)
d_sd_of_beta_II <- density(sd_of_beta_II,na.rm=T)

#application of the Functional Normalization
preprocessFunnorm_results <- preprocessFunnorm(RGset)
beta_Funnorm <- getBeta(preprocessFunnorm_results)
beta_Funnorm_I <- beta_Funnorm[rownames(beta_Funnorm) %in% dfI$IlmnID,]
beta_Funnorm_II <- beta_Funnorm[rownames(beta_Funnorm) %in% dfII$IlmnID,]


#calculation of the mean, sd and their corresponding density for the normalized data
mean_of_beta_Funnorm_I <- apply(beta_Funnorm_I,1,mean,na.rm=T)
mean_of_beta_Funnorm_II <- apply(beta_Funnorm_II,1,mean,na.rm=T)
d_mean_of_beta_Funnorm_I <- density(mean_of_beta_Funnorm_I,na.rm=T)
d_mean_of_beta_Funnorm_II <- density(mean_of_beta_Funnorm_II,na.rm=T)

sd_of_beta_Funnorm_I <- apply(beta_Funnorm_I,1,sd,na.rm=T)
sd_of_beta_Funnorm_II <- apply(beta_Funnorm_II,1,sd,na.rm=T)
d_sd_of_beta_Funnorm_I <- density(sd_of_beta_Funnorm_I,na.rm=T)
d_sd_of_beta_Funnorm_II <- density(sd_of_beta_Funnorm_II,na.rm=T)

par(mfrow=c(2,3))
targets$Group <- as.factor(targets$Group)
palette(c("#E50068", "#0067E6"))

#Raw data graphs
plot(d_mean_of_beta_I,col="#77FF00",main="raw Beta density", lwd=2.5)
lines(d_mean_of_beta_II,col="#FF7038",lwd=2.5)
legend('topright', legend=c("Type I","Type II"), fill = c("#77FF00","#FF7038"), cex=1)
plot(d_sd_of_beta_I,col="#77FF00",main="raw sd density", lwd=2.5)
lines(d_sd_of_beta_II,col="#FF7038",lwd=2.5)
legend('topright', legend=c("Type I","Type II"), fill = c("#77FF00","#FF7038"), cex=1)
boxplot(beta, col = targets$Group)
title('Boxplot of raw Beta values')


#Normalized data graphs
plot(d_mean_of_beta_Funnorm_I,col="#77FF00",main="preprocessFunnorm Beta density",lwd=2.5)
lines(d_mean_of_beta_Funnorm_II,col="#FF7038",lwd=2.5)
legend('topright', legend=c("Type I","Type II"), fill = c("#77FF00","#FF7038"), cex=1)
plot(d_sd_of_beta_Funnorm_I,col="#77FF00",main="preprocessFunnorm sd density",lwd=2.5)
lines(d_sd_of_beta_Funnorm_II,col="#FF7038",lwd=2.5)
legend('topright', legend=c("Type I","Type II"), fill = c("#77FF00","#FF7038"), cex=1)
boxplot(beta_Funnorm, col = targets$Group)
title('Boxplot of normalized Beta values')



pca_results <- prcomp(t(beta_Funnorm),scale=T)
par(mfrow=c(1,1))


library(factoextra)
fviz_eig(pca_results, addlabels = T,xlab='PC number',ylab='% of variance', barfill = "#0063A6", barcolor = "black")

targets$Group <- as.factor(targets$Group)
palette(c("#E50068", "#0067E6"))
plot(pca_results$x[,1], pca_results$x[,2],cex=1,pch=19,col=targets$Group,xlab="PC1",ylab="PC2",xlim=c(-700,700),ylim=c(-700,700), main=' PCA (Groups)')
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.4,pos=3)
legend("bottomright",legend=levels(targets$Group),col=c(1,2),pch=19, cex=1.0)


targets$Sex <- as.factor(targets$Sex)
palette(c("#C364CA", "#8BC4F9"))
plot(pca_results$x[,1], pca_results$x[,2],cex=1,pch=19,col=targets$Sex,xlab="PC1",ylab="PC2",xlim=c(-700,700),ylim=c(-700,700), main='PCA (Sex)')
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.4,pos=3)
legend("bottomright",legend=levels(targets$Sex),col=c(1:nlevels(targets$Group)),pch=19,cex=1.0)



targets$Slide <- as.factor(targets$Slide)
palette(c("#9e0059", "green"))
plot(pca_results$x[,1], pca_results$x[,2],cex=1,pch=19,col=targets$Slide,xlab="PC1",ylab="PC2",xlim=c(-700,700),ylim=c(-700,700), main=' PCA (Batch)')
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.4,pos=3)
legend("bottomright",legend=levels(targets$Slide),col=c(1,2),pch=19, cex=1.0)

library(future.apply) #for parallelization
plan(multisession)
My_MannWhitney_function <- function(x) {
  wilcox <- wilcox.test(x~ targets$Group)
  return(wilcox$p.value)
} 
pValues_wilcox <- future_apply(beta_Funnorm,1, My_MannWhitney_function)

final_wilcox<- data.frame(beta_Funnorm, pValues_wilcox)

final_wilcox <- final_wilcox[order(final_wilcox$pValues_wilcox),]

table(final_wilcox$pValues_wilcox)



#multiple testing correction
corrected_pValues_BH <- p.adjust(final_wilcox$pValues_wilcox,"BH")
corrected_pValues_Bonf <- p.adjust(final_wilcox$pValues_wilcox,"bonferroni")

#creation of a dataframe that contains the nominal p-values and the p-values after both corrections.
final_wilcox_corrected <- data.frame(final_wilcox, corrected_pValues_BH, corrected_pValues_Bonf)

#differentially methylated probes according to the threshold (both before and after the corrections)
before_correction <- nrow(final_wilcox_corrected[final_wilcox_corrected$pValues_wilcox<= 0.05,])
after_Bonferroni <- nrow(final_wilcox_corrected[final_wilcox_corrected$corrected_pValues_Bonf <= 0.05,])
after_BH <- nrow(final_wilcox_corrected[final_wilcox_corrected$corrected_pValues_BH <= 0.05,])



diff_meth_df <- data.frame(before_correction, after_Bonferroni, after_BH)
rownames(diff_meth_df) <- c("# differentially methylated probes")
diff_meth_df

par(mfrow=c(1,1))
boxplot(final_wilcox_corrected [,9:11], col = c("#7038FF", "#C7FF38", "seashell2"), names=NA)
legend("bottomright", legend=c("raw", "BH", "Bonferroni"),col=c("#7038FF", "#C7FF38", "seashell2"), lty=1:1, cex=0.8, xpd=TRUE, pch=15)

# WT group mean
WT_group <- final_wilcox_corrected[,targets$Group=="WT"]
WT_group_mean <- apply(WT_group, 1, mean)

# MUT group mean
MUT_group <- final_wilcox_corrected[,targets$Group=="MUT"]
MUT_group_mean <- apply(MUT_group, 1, mean)

# Compute delta
delta <- WT_group_mean - MUT_group_mean

#create a dataframe
toVolcPlot <- data.frame(delta, -log10(final_wilcox_corrected$pValues_wilcox))
head(toVolcPlot)

par(mfrow=c(1,1))
HighLight <- toVolcPlot[abs(toVolcPlot[,1])>0.1 & toVolcPlot[,2]>(-log10(0.05)),]
plot(toVolcPlot[,1],toVolcPlot[,2], pch=19, cex=0.6 ,xlab="delta",ylab="-log10(Pvalue)", col="black")
abline(a=-log10(0.05),b=0, col="#fb8b24")
points(HighLight[,1], HighLight[,2], pch=17, cex=0.6, col="#00f5d4")

library(qqman)
final_wilcox_corrected_df <- data.frame(rownames(final_wilcox_corrected),final_wilcox_corrected)
colnames(final_wilcox_corrected_df)[1] <- "IlmnID"

final_wilcox_corrected_annotated <- merge(final_wilcox_corrected_df, Illumina450Manifest_clean,by="IlmnID")

input_Manhattan <- data.frame(ID=final_wilcox_corrected_annotated$IlmnID,CHR=final_wilcox_corrected_annotated$CHR, MAPINFO=final_wilcox_corrected_annotated$MAPINFO, PVAL=final_wilcox_corrected_annotated$pValues_wilcox)
head(input_Manhattan)

levels(input_Manhattan$CHR)[levels(input_Manhattan$CHR) == "X"] <- "23"               
levels(input_Manhattan$CHR)[levels(input_Manhattan$CHR) == "Y"] <- "24"              
input_Manhattan$CHR <- as.numeric(as.character(input_Manhattan$CHR))


manhattan(input_Manhattan, snp="ID",chr="CHR", bp="MAPINFO", p="PVAL",suggestiveline= -log10(0.05),col=rainbow(24), annotateTop = F)


library(gplots)
library(viridis)
input_heatmap=as.matrix(final_wilcox[1:100,1:8])


group_color = c()
i = 1
for (name in colnames(beta)){
  if (name %in% wt){group_color[i]="#0067E6"}
  else{group_color[i]="#E50068"}
  i = i+1
}

heatmap.2(input_heatmap,col=viridis(100),Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=group_color,density.info="none",trace="none",scale="none",symm=F,main="Complete linkage",key.xlab='beta-val',key.title=NA,keysize=1,labRow=NA)
legend("topright", legend=levels(targets$Group),col=c('#E50068','#0067E6'),pch = 19,cex=0.7)


heatmap.2(input_heatmap,col=viridis(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'single'),dendrogram="both",key=T,ColSideColors=group_color,density.info="none",trace="none",scale="none",symm=F,main="Single linkage",key.xlab='beta-val',key.title=NA,keysize=1,labRow=NA)
legend("topright", legend=levels(targets$Group),col=c('#E50068','#0067E6'),pch = 19,cex=0.7)


heatmap.2(input_heatmap,col=viridis,Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'average'),dendrogram="both",key=T,ColSideColors=group_color
          ,density.info="none",trace="none",scale="none",symm=F,main="Average linkage",key.xlab='beta-val',key.title=NA,keysize=1,labRow=NA)
legend("topright", legend=levels(targets$Group),col=c('#E50068','#0067E6'),pch = 19,cex=0.7)
