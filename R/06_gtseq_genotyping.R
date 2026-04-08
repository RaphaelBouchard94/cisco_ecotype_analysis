rm(list=ls())
##################
#GENOTYPING
##################

#set working directory, change working directory as needed
source("GTscore/GTscore.R")
library(magrittr)
library(tidyverse)

###############################################
### Generate data using GTscore.R functions ###
###############################################

#SHOW GENOTYPING USING TEST DATASET FOR SPEED
singleSNP_locusTable<-read.delim("GTscore/LocusTable_singleSNPs.txt",header=TRUE)
singleSNP_alleleReads<-read.delim("GTscore/AlleleReads_singleSNPs.txt",header=TRUE,row.names=1)

#generate singleSNP genotypes
#polyGenResults_singleSNP<-polyGen(singleSNP_locusTable,singleSNP_alleleReads)

#write.table(polyGenResults_singleSNP,"GTscore_results/polyGenResults_singleSNP.txt",quote=FALSE,sep="\t")

#singleSNP_summary<-summarizeGTscore(singleSNP_alleleReads,singleSNP_locusTable, polyGenResults_singleSNP)

#write.table(singleSNP_summary,"GTscore_results/singleSNP_summary.txt",quote=FALSE,sep="\t",row.names=FALSE)

#summarize sample genotype rate
#calculate genotype rate per sample for single SNP data

#sample_genotypeRate_singleSNP<-sampleGenoRate(polyGenResults_singleSNP)

#write.table(sample_genotypeRate_singleSNP,"GTscore_results/sample_genotypeRate_singleSNP.txt",quote=FALSE,sep="\t")

#singleSNP_sampleSummary<-summarizeSamples(polyGenResults_singleSNP,singleSNP_alleleReads)

#write.table(singleSNP_sampleSummary,"GTscore_results/singleSNP_sampleSummary.txt",quote=FALSE,sep="\t")


####################
## DATA SUMMARIES ##
####################

polyGenResults_singleSNP <- read.table("data/gtscore_analysis/polyGenResults_singleSNP.txt",header=T)

##

singleSNP_summary <- read.table("data/gtscore_analysis/singleSNP_summary.txt",header=T)

##

singleSNP_sampleSummary <- read.table("data/gtscore_analysis/singleSNP_sampleSummary.txt",header=T)

singleSNP_sampleSummary %<>% mutate(sample_id =str_replace_all(str_replace_all(sample, "\\.","_"),"X02_input_files_",""))

##

sample_genotypeRate_singleSNP <- read.table("data/gtscore_analysis/sample_genotypeRate_singleSNP.txt", header=T)

sample_genotypeRate_singleSNP %<>% mutate(sample_id = str_replace_all(str_replace_all(sample, "\\.","_"),"X02_input_files_",""))

##

GTscore_individualSummary<-read.delim("GTscore/GTscore_individualSummary.txt",header=TRUE,stringsAsFactors=FALSE)

GTscore_individualSummary %<>% mutate(sample_id = str_replace_all(str_replace_all(Sample, "-","_"),"02_input_files/",""))


#combine individual summary data with sample genotype rate
GTscore_individualSummary<-merge(GTscore_individualSummary,sample_genotypeRate_singleSNP,by.x="sample_id",by.y="sample_id")

GTscore_individualSummary<-merge(GTscore_individualSummary,singleSNP_sampleSummary,by.x="sample_id",by.y="sample_id")

#-------------#
#SNP coverage
#-------------#

ggplot(GTscore_individualSummary)+
  geom_density(aes(x = Total.Reads), color = "red")+
  geom_density(aes(x = Primer.Probe.Reads), color = "blue")+
  theme_bw(base_size = 16)

ggsave("figures/99_primer_probe.png")

mean(GTscore_individualSummary$Primer.Probe.Proportion)

sd(GTscore_individualSummary$Primer.Probe.Proportion)

ggplot(singleSNP_summary, aes(x = reorder(Locus_ID,AvgReadDepth), y = log(AvgReadDepth)))+
  geom_point()+
  theme_classic()

ggplot(singleSNP_summary, aes(x= AvgReadDepth))+
  geom_density(fill = "lightgrey")+
  geom_vline(aes(xintercept=mean(AvgReadDepth, na.rm = T)),color = "darkred", linetype = "dashed")+
  labs(title = "Average read depth = 37X")+
  theme_bw()

#ggsave("figures/read_depth_per_locus.png")

mean(singleSNP_summary$AvgReadDepth, na.rm = T)
#Average read depth 36.71
sd(singleSNP_summary$AvgReadDepth, na.rm = T)
#sd 22

length(which(singleSNP_summary$AvgReadDepth < 10))

#What is the effect of read depth on genotyping rate

singleSNP_summary %>% ggplot(aes(x =AvgReadDepth ,y=GenotypeRate))+
  geom_point()+
  geom_vline(xintercept = 10)

#ggsave("figures/avgreaddepth_genotyperate.jpeg")


poorQualitySNPs <- singleSNP_summary %>% filter(AvgReadDepth < 10 | GenotypeRate < 0.75 | conScore> 0.5) %>% pull(Locus_ID)

#A total of 67 SNPs which are categorized as poor quality based on the threshold above


#---------------------------#
# Identify problematic individuals 
# to remove from the dataset
#---------------------------#

mean(GTscore_individualSummary %>% filter(str_detect(sample_id,"m")) %>% pull(GenotypeRate.x))

GTscore_individualSummary %<>%  
  separate(sample.x, c("sample", "barcode"), sep = "_") %>% 
  mutate(sample_clean = str_replace_all(sample,"\\.","_")) %>% 
  dplyr::select(-sample.y,-GenotypeRate.y)


ggplot(GTscore_individualSummary,aes(x = conScore))+
  geom_density()

threshold_conscore <- quantile(GTscore_individualSummary$conScore,probs = 0.975,na.rm =T)
#0.25

ggplot(GTscore_individualSummary,aes(x = conScore))+
  geom_density()+
  geom_vline(xintercept = threshold_conscore)

length(which(GTscore_individualSummary$conScore > threshold_conscore))

#39 individuals with conScore > 0.25 (95% conf interval)

#What is the deal those 39 individuals

GTscore_individualSummary %<>% mutate(species = if_else(str_detect(sample_id,"coar"),"cisco","whitefish"))

GTscore_individualSummary %>% ggplot(aes(x = Heterozygosity, y = conScore, col = GenotypeRate.x,shape = species))+
  geom_point(alpha = 0.6,size = 2)+
  scale_color_viridis_c()+
  geom_hline(yintercept = threshold_conscore)+
  theme_bw()

GTscore_individualSummary %>% filter(conScore > threshold_conscore) %>% 
  mutate(high_het = if_else(Heterozygosity < 0.2,"no","yes")) %>% 
  group_by(high_het) %>% 
  count()

#It seems like individuals that have a high conScore also have EITHER HIGHER OR LOWER heterozygosity.
#These individuals are most likely whitefish and true contaminated samples. 
#The other individuals that have high conScore
#and high heterozygosity are most likely contaminated.

#What about the genotyping rate?

GTscore_individualSummary %>% ggplot(aes(x = GenotypeRate.x, y = Heterozygosity, col = conScore))+
  geom_point(alpha = 0.6,size = 2)+
  scale_color_viridis_c()+
  theme_bw()

#ggsave("figures/het_genotyperate.jpeg")

GTscore_individualSummary %>% ggplot(aes(x = Total.Reads, y = Heterozygosity, col = conScore))+
  geom_point(alpha = 0.6,size = 2)+
  scale_color_viridis_c()+
  theme_bw()


GTscore_individualSummary %>% ggplot(aes(x = Total.Reads, y = Heterozygosity, col = species))+
  geom_point(alpha = 0.6,size = 2)+
  theme_bw()

#Test if cisco samples identified in whitefish run are truly ciscos

cocl <- read.table("data/gtscore_analysis/whitefish_samples",header = T)

GTscore_individualSummary %>% 
  ggplot(aes(x = Total.Reads, y = Heterozygosity, col = if_else(sample_id %in% cocl$x, "red","black")))+
  geom_point(alpha = 0.6,size = 2)+
  theme_bw()

#Those samples were not cisco afterall (471 whitefish at least)

#ggsave("figures/het_totalreads.jpeg")

contaminatedSample <- GTscore_individualSummary %>% filter(Heterozygosity < 0.15) %>% dplyr::select(sample_id)

contaminatedSample %<>% mutate(FISHES_ID = str_replace_all(sample_id,"_S.*L001_R1_001","")) %>% pull(FISHES_ID)

write.table(contaminatedSample, "data/gtscore_analysis/contaminated_samples.txt", quote = F, row.names = F, col.names = F)

# 461 samples are removed based on this threshold

#------------------#
## EXPORT RUBIAS ##
#------------------#

data <- read.csv("data/gtscore_analysis/cisco_to_sequence.csv")

data %<>% mutate(fishes_id = str_replace_all(fishes_id, "-","_"))

data %<>% distinct()

sampleMetaData_final <- data.frame(sample_type=if_else(data$stock_clean == "Source","reference","mixture"),
                                   repunit=NA,
                                   collection=str_replace_all(data$river.name," ",""),
                                   indiv=data$fishes_id)


sampleMetaData_final %>% group_by(sample_type) %>% dplyr::count()

polyGenResults_singleSNP <- read.table("data/gtscore_analysis/polyGenResults_singleSNP.txt",header=T)

# Rename the columns, matching either coar or cocl, and allowing for underscores
colnames(polyGenResults_singleSNP) <- gsub("\\.", "_", colnames(polyGenResults_singleSNP))

##

exportRubias(polyGenResults_singleSNP,
             singleSNP_locusTable,
             sampleMetaData_final,
             sampleBlacklist=contaminatedSample,
             locusBlacklist=poorQualitySNPs,
             filename="data/gtscore_analysis/polyGenResults_singleSNP_rubias.txt")

exportRubias(polyGenResults_singleSNP,
             singleSNP_locusTable,
             sampleMetaData_final,
             sampleBlacklist=contaminatedSample,
             locusBlacklist=poorQualitySNPs,
             filename="data/gtscore_analysis/olyGenResults_singleSNP_rubias.txt")

####################
## EXPORT GENEPOP ##
####################

source_mixed_id <- read.table("data/gtscore_analysis/polyGenResults_singleSNP_rubias.txt",header = T)

source_id <- source_mixed_id %>% filter(sample_type == "reference") %>% pull(indiv)

#Export source samples only

exportGenepop(polyGenResults_singleSNP,
              singleSNP_locusTable, 
              locusBlacklist=poorQualitySNPs,
              sampleWhitelist = source_id,
              sampleBlacklist = contaminatedSample,
              filename="data/gtscore_analysis/polyGenResults_singleSNP_genepop_source.txt")

exportGenepop(polyGenResults_singleSNP,
              singleSNP_locusTable, 
              locusBlacklist=poorQualitySNPs,
              sampleWhitelist = source_id,
              sampleBlacklist = contaminatedSample,
              filename="data/gtscore_analysis/polyGenResults_singleSNP_genepop_source.gen")

