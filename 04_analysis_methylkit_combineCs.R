#### Packages and working DIR ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(methylKit)
library(ggplot2)
library(devtools)
install_github("drveera/ggman")
library(ggman)
setwd("~/projects/methylation_arsenic/veronika_protocol/04_stats")

#### Input and file check ####
#Inputing BismarkCytosineReport because that has the strand information, cov-files miss that.

##Make file list
file.list = list("~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample1_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample2_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample3_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample4_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample5_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample6_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample7_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample8_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample9_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample10_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample11_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample12_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample13_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample14_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample15_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample16_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample17_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample18_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample19_val_1_bismark_bt2_pe.CpG_report.txt",
                 "~/projects/methylation_arsenic/veronika_protocol/03_methylation_calling/sample20_val_1_bismark_bt2_pe.CpG_report.txt")

##read the cytosine report files, mincov is 10x by default
myobj=methRead( file.list,pipeline="bismarkCytosineReport",
                sample.id=list("ctrl1_F","ctrl2_M", "ctrl3_F","ctrl4_F","ctrl5_M","ctrl6_M","ctrl7_F","ctrl8_M","ctrl9_M","ctrl10_F",
                               "test11_F","test12_F","test13_M","test14_M","test15_M","test16_F","test17_F","test18_F", "test19_M","test20_M"),
                assembly="p.major1.1",
                treatment=c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
                context = "CpG")

getData(myobj[[1]])
str(myobj@.Data)

##check meth pattern
getMethylationStats(myobj[[15]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj[[15]],plot=TRUE,both.strands=FALSE)

##percentile filtering
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

str(filtered.myobj@.Data)

##unite samples in table, destrand set to TRUE as there is strand info to combine Cs
meth=unite(filtered.myobj, destrand=TRUE, min.per.group = 8L)
head(meth)

##explore and visualise data
#getCorrelation(meth,plot=TRUE) #not working

clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
hc = clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)
plot(hc, main=NULL, xlab=NULL)
# using dendrogram objects
hcd = as.dendrogram(hc)
# alternative way to get a dendrogram
plot(hcd)

PCASamples(meth, adj.lim=c(0.9,0.2))
PCA1= PCASamples(meth, adj.lim=c(0.9,0.2), obj.return = TRUE)
summary(PCA1)

#without 14 which seems to be an outlier
PCASamples(meth2, adj.lim=c(0.9,0.2))
PCA2 = PCASamples(meth2, adj.lim=c(0.9,0.2), obj.return = TRUE)
summary(PCA2)

PCASamples(meth2, adj.lim=c(0.9,0.2), comp = c(3,4))
PCA3=PCASamples(meth2, adj.lim=c(0.9,0.2), comp = c(3,4), obj.return = TRUE)
summary(PCA3)

#### DiffMethylation with DSS ####
#remove sample 14
meth2 =reorganize(meth,sample.ids=c("ctrl1_F","ctrl2_M", "ctrl3_F","ctrl4_F","ctrl5_M","ctrl6_M","ctrl7_F","ctrl8_M","ctrl9_M","ctrl10_F","test11_F","test12_F","test13_M","test15_M","test16_F","test17_F","test18_F","test19_M","test20_M"),
                  treatment=c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1))
dataAll_meth2 <- (getData(meth2))

#Covariates doesn't work with DSS
#covariates = data.frame(sex=c("F","M","F","F","M","M","F","M","M","F","F","F","M","M","F","F","F","M","M"))

#run with all samples
all=calculateDiffMethDSS(meth, adjust="bonferroni")
myDiff10p=getMethylDiff(all,difference=10,qvalue=0.05)

#run with 14 removed
mydiff=calculateDiffMethDSS(meth2, adjust="bonferroni")
myDiff10p=getMethylDiff(mydiff,difference=10,qvalue=0.05)
head(myDiff10p)


##make results as dataframe
dataAll_myDiff <- (getData(mydiff))
dataAll_myDiff10p <- (getData(myDiff10p))

##plotting
##QQ-plot
qq(mydiff$pvalue)
estlambda(dataAll_myDiff[,"pvalue"],plot=F)

### Volcano plot visualization 
results = as.data.frame(dplyr::mutate(as.data.frame(dataAll_myDiff), significance=ifelse(dataAll_myDiff$qvalue<0.05, "<0.05", "Not Sig")), row.names=rownames(dataAll_myDiff))
volcano = ggplot2::ggplot(results, ggplot2::aes(meth.diff, -log10(pvalue))) +
  ggplot2::geom_point(ggplot2::aes(col = significance))  +
  ggplot2::scale_color_manual(values = c( "darkorange","black")) +
  #ggplot2::ggtitle("DSS")+
  ggplot2::labs(x = "Delta methylation ( test - control )")

volcano2 = volcano + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))


volcano3 = volcano2 + geom_vline(xintercept=c(-10, 10), col="red") +
  geom_hline(yintercept=7.11, col="red")

ggsave("rplot.pdf", volcano3)

dev.off()

#Manhattan plot, need to change the chrs names..
NCBI = read.table(file="NCBI_chr.txt", header = TRUE)

#Join data. Retain only rows in both sets.
dataAll_join = dplyr::inner_join(dataAll_myDiff, NCBI, by = "chr")
dataAll_join10p = dplyr::inner_join(dataAll_myDiff10p, NCBI, by = "chr")

#For checking. All rows in a that do not have a match in b. Mito has new NCBI id.
#mismatch = dplyr::anti_join(dataAll_myDiff, NCBI, by = "chr")

#making the plot
data <- tibble::rowid_to_column(dataAll_join, "ID")
manhattan<-ggman(data, snp = "ID", bp = "start", chrom = "chr_nr", pvalue = "pvalue", 
                 sigLine = 7.11,logTransform = TRUE, xlabel = "chromosome",
                 ylabel = "-log10(pvalue)", pointSize = 1, title = "")+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,panel.background= element_blank()
    ,axis.text.x=element_text(size=6)
    ,axis.text.y=element_text(size=10))

manhattan2 <- manhattan + 
  scale_color_manual(values=c("darkslategray3","darkslategray4"))

ggsave("manhattan_14removed_Cscomb_DSS_PUBS.pdf", manhattan2)

##Save results##
write.table(dataAll_join, "DSS_arsenic_Allresults_sample14removed_Cscombined.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dataAll_join10p, "DSS_arsenic_topresults_sample14removed_Cscombined.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dataAll_meth2, "methinput_arsenic_sample14removed_Cscombined.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#### Methylkit's diffmeth #####
my.diffMeth1<-calculateDiffMeth(meth2,
                                overdispersion="MN",
                                adjust = "bonferroni",
                                test = "F")

#covariates
my.diffMeth2<-calculateDiffMeth(meth2,
                                covariates=covariates,
                                overdispersion="MN",
                                adjust = "bonferroni",
                                test = "F")

my.diffMeth1_10p=getMethylDiff(my.diffMeth1,difference=10,qvalue=0.05)
my.diffMeth2_10p=getMethylDiff(my.diffMeth2,difference=10,qvalue=0.05)

head(my.diffMeth1_10p)
head(my.diffMeth2_10p)

dataAll_myDiff_mk1 <- (getData(my.diffMeth1))
dataAll_myDiff_mk2 <- (getData(my.diffMeth2))

estlambda(dataAll_myDiff_mk1[,"pvalue"],plot=F)
estlambda(dataAll_myDiff_mk2[,"pvalue"],plot=F)

### Volcano plot visualization 
results_mk1 = as.data.frame(dplyr::mutate(as.data.frame(dataAll_myDiff_mk1), significance=ifelse(dataAll_myDiff_mk1$qvalue<0.05, "<0.05", "Not Sig")), row.names=rownames(dataAll_myDiff_mk1))
volcano_mk1 = ggplot2::ggplot(results_mk1, ggplot2::aes(meth.diff, -log10(pvalue))) +
  ggplot2::geom_point(ggplot2::aes(col = significance))  +
  ggplot2::scale_color_manual(values = c( "darkorange","black")) +
  ggplot2::ggtitle("Methylkit")+
  ggplot2::labs(x = "Delta methylation ( control vs test )")

volcano_mk1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))


#### Sex and batch effect ####
#Make new data object because no NAs allowed
sampleAnnotation=data.frame(lane_id=c("1","1","1","1","1","2","2","2","2","2","1","1","1","1","2","2","2","2","2"),
                            sex=c("F","M","F","F","M","M","F","M","M","F","F","F","M","M","F","F","F","M","M"))

sampleAnnotation_no16_3=data.frame(lane_id=c("1","1","1","1","2","2","2","2","2","1","1","1","1","2","2","2","2"),
                            sex=c("F","M","F","M","M","F","M","M","F","F","F","M","M","F","F","M","M"))


meth_as_all=unite(filtered.myobj, destrand=TRUE)
meth2_as =reorganize(meth_as_all,sample.ids=c("ctrl1_F","ctrl2_M", "ctrl3_F","ctrl4_F","ctrl5_M","ctrl6_M","ctrl7_F","ctrl8_M","ctrl9_M","ctrl10_F","test11_F","test12_F","test13_M","test15_M","test16_F","test17_F","test18_F","test19_M","test20_M"),
                  treatment=c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1))

#Remove the two outlier inds
sampleAnnotation_no16_3=data.frame(lane_id=c("1","1","1","1","2","2","2","2","2","1","1","1","1","2","2","2","2"),
                                   sex=c("F","M","F","M","M","F","M","M","F","F","F","M","M","F","F","M","M"))

meth2_as_no16_3 =reorganize(meth_as_all,sample.ids=c("ctrl1_F","ctrl2_M","ctrl4_F","ctrl5_M","ctrl6_M","ctrl7_F","ctrl8_M","ctrl9_M","ctrl10_F","test11_F","test12_F","test13_M","test15_M","test17_F","test18_F","test19_M","test20_M"),
                     treatment=c(0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1))

##Association test
as=assocComp(mBase=meth2_as,sampleAnnotation)
as=assocComp(mBase=meth2_as_no16_3,sampleAnnotation_no16_3)
as

#remove pc3 and do diffmeth analysis without it, no NAs allowed
newObj=removeComp(meth2_as,comp=3)
mydiff_as=calculateDiffMethDSS(meth2_as, adjust="bonferroni")
myDiff10p_as=getMethylDiff(mydiff_as,difference=10,qvalue=0.05)
head(myDiff10p_as)
dataAll_myDiff_as <- (getData(mydiff_as))
dataAll_join = dplyr::inner_join(dataAll_myDiff_as, NCBI, by = "chr")

dataAll_myDiff10p_as <- (getData(myDiff10p_as))
dataAll_join10p_as = dplyr::inner_join(dataAll_myDiff10p_as, NCBI, by = "chr")

write.table(dataAll_join, "DSS_arsenic_Allresults_sample14removed_Cscombined_PC3removed.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dataAll_join10p_as, "DSS_arsenic_topresults_sample14removed_Cscombined_PC3removed.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Draw PCA figures
PCASamples(meth2_as_no16, adj.lim=c(0.9,0.2), comp = c(2,3))
PCASamples(meth2_as, adj.lim=c(0.9,0.2), comp = c(2,3))


####Coverage check####
Cov <- dataAll_meth2 %>% 
  #rowwise will make sure the sum operation will occur on each row
  rowwise() %>% 
  #then a simple sum(..., na.rm=TRUE) is enough to result in what you need
  mutate(covsum = sum(coverage1, coverage2, coverage3, coverage4, coverage5, coverage6, coverage7, coverage8, coverage9, coverage10, coverage11, coverage12, coverage13, coverage14, coverage15, coverage16, coverage17, coverage18, coverage19, na.rm=TRUE))

Cov<- dplyr::select(Cov,chr, start, covsum)
Cov <- Cov %>% 
  tidyr::unite(chr_start, c("chr", "start"))

dataAll_meth2 <- dataAll_meth2 %>% 
  tidyr::unite(chr_start, c("chr", "start"))

DSS<- dplyr::select(dataAll_myDiff,chr, start, qvalue)
DSS <- DSS %>% 
  tidyr::unite(chr_start, c("chr", "start"))

DSS_cov <- dplyr::inner_join(DSS, dataAll_meth2, by = "chr_start")
DSS_cov$log10 = log10(DSS_cov$qvalue)


plot_all <- ggplot(DSS_cov, aes(x = coverage1, y = log10)) +
  geom_point()
                             