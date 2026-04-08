#########################
####### LIBRARY #########
#########################
library(tidyverse)
library(magrittr)
library(paletteer)
library(patchwork)

####################################
######  Admixture analysis  ########
####################################
# Load metadata
metadata <- read.table("data/id_lowcov_fulldata.tsv")
metadata2 <- read.csv("data/JB_samples_metadata.csv")
metadata %<>% left_join(metadata2, by = c("V1" = "sample_id"))
metadata %<>% distinct()

# Import ngsadmix proportion of admixture
files <- list.files(path = "data/ngsadmix/all_pop/", 
                    pattern = "qopt",
                    full.names = TRUE)

# Define color ONLY for NOT_summer dominant ancestry
not_summer_color <- "#E69F00"
# Define colors for dominant ancestry in each population
pop_colors <- c(
  "NOT_summer" = "#E69F00",
  "NOT_fall" = "#0072B2",
  "RUP" = "#56B4E9",
  "EAS" = "lightpink4", 
  "SAB" = "lightpink4",
  "LAG" = "#0072B2")
# Define population order for faceting
pop_order <- c("LAG","SAB","EAS","RUP", "NOT_fall", "NOT_summer")
metadata2$pop <- factor(metadata2$pop, levels = pop_order)

kfig <- vector('list', length(files) + 1)
for (i in 2:(length(files))) {
  df <- read.table(paste("data/ngsadmix/all_pop/all_maf0.01_pctind0.8_maxdepth10_pruned_singletons_", i, "_1.qopt", sep = ""))
  
  df$id <- metadata$V1
  
  # Join with metadata FIRST
  df <- df %>% 
    left_join(metadata, by = c("id" = "V1"))
  
  # Calculate max prob for each individual (before reshaping)
  # Get column names that start with V
  v_cols <- df %>% select(starts_with("V")) %>% names()
  
  df <- df %>%
    mutate(max_prob = do.call(pmax, select(., all_of(v_cols))))
  
  # Create ordering: within each pop, order by descending max_prob
  id_order <- df %>%
    arrange(pop, desc(max_prob)) %>%
    pull(id)
  # Then reshape
  dfclean <- df %>% 
    gather('ancestry_cluster', 'prob', starts_with("V")) %>%
    group_by(id) %>% 
    mutate(likely_assignment = ancestry_cluster[which.max(prob)],
           assignment_prob = max(prob)) %>% 
    ungroup() %>%
    mutate(id = factor(id, levels = id_order))
  
  kfig[[i]] <-
    ggplot(dfclean, aes(id, prob, fill = ancestry_cluster)) +
    geom_col(width = 1) +
    facet_grid(~ pop, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = paletteer_d("ggsci::category20_d3")) +
    labs(y = paste("K =", i), x = "") +
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 18),
          axis.ticks.x = element_blank(),
          strip.text = if (i > 2) element_blank() else element_text(size = 18),
          legend.position = "none")
}

wrap_plots(kfig[2:4], ncol = 1)

####################################
#Check best cluster

log_lik <- read.table("data/ngsadmix/all_pop/best_like_final.tsv")

colnames(log_lik) <- c("K","iter","like")

evanno = as.data.frame(
  cbind(1:max(log_lik$K),
        tapply(log_lik$like, log_lik$K, FUN= function(x) mean(abs(x))/sd(abs(x))))
)

log_lik %>%
  group_by(K) %>%
  dplyr::summarise(
    sd = sd(like, na.rm = TRUE),
    like = mean(like)
  ) -> llogs

#Second order change
llogs$llike = NA 
for(i in 2:nrow(llogs)){
  llogs$llike[i] = llogs$like[i] - llogs$like[i-1]
}

#Third order
llogs$lllike = NA 
for(i in 2:nrow(llogs)){
  llogs$lllike[i] = llogs$llike[i+1] - llogs$llike[i]
}

p1 = ggplot(data = llogs,
            aes(x = as.factor(K), y = like/1e8, ymin = like/1e8-sd/1e8, ymax = like/1e8+sd/1e8)) +
  geom_point(size = 1) +
  #geom_line() +
  geom_errorbar(width = 1) +
  coord_cartesian(xlim = c(0,10)) +
  labs(x = "K",
       y = "L(K)")+
  theme_classic()

p2 = ggplot(data = llogs,
            aes(x = as.factor(K), y = llike, ymin = llike-sd, ymax = llike+sd)) +
  geom_point() +
  #geom_line() +
  geom_errorbar(width = 0.5) +
  coord_cartesian(xlim = c(0,10)) +
  labs(x = "K",
       y = "L'(K)")+
  theme_classic()

p3 = ggplot(data = llogs,
            aes(x = as.factor(K), y = abs(lllike), ymin = abs(lllike)-sd, ymax = abs(lllike)+sd)) +
  geom_point() +  #geom_line() +
  geom_errorbar(width = 0.5) +
  labs(x = "K",
       y = "|L''(K)|")+
  theme_classic()

p4 = ggplot(data = llogs,
            aes(x = as.factor(K), y = abs(lllike)/sd,group = 1))+
  geom_line() +
  geom_point() +
  labs(x = "K",
       y = expression(Delta*K))+
  scale_y_log10() +
  theme_classic()

p1/p2/p3/p4

ggsave("figures/06_all_pop_bestK.jpeg")


####################################

for (i in 1:length(files)){
  new_df <-  read.table(files[i])
  assign(paste("k",substr(files[i],74,75),sep=""),new_df)
}

#K2

k_2$id <- metadata$V1

k2 <- k_2 %>% 
  gather('pop', 'prob', V1:V2) %>% 
  dplyr::group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assignment_prob = max(prob)) %>%  # Fixed typo: assingment -> assignment
  ungroup() %>%
  # Sort first by major cluster, then by assignment probability
  dplyr::arrange(likely_assignment, desc(assignment_prob)) %>%
  mutate(id = forcats::fct_inorder(as.factor(id)))

k2 %<>% left_join(metadata2,by = c("id" = "sample_id")) %>% distinct()


k2 <- k2 %>%
  arrange(pop.y, likely_assignment, desc(assignment_prob)) %>%
  mutate(id = forcats::fct_inorder(as.factor(id)))

a <-ggplot(k2 %>% filter(!is.na(pop.y)), aes(id, prob, fill = pop.x)) +
  geom_col(width = 2) +
  scale_fill_manual(values=c("#0072B2","#E69F00"))+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x",
             labeller =  as_labeller(c("LAG"="LAG fall-run", 
                                       "SAB"="SAB fall-run", 
                                       "EAS"="EAS fall-run", 
                                       "RUP"="RUP fall-run", 
                                       "NOT_fall"="NOT fall-run", 
                                       "NOT_summer"="Kuukamek")))+
  ylab("K2")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size=18),
        legend.position = "none")

a

####################################
#K3

k_3$id <- metadata$V1

k3 <- k_3 %>% 
  gather('pop', 'prob', V1:V3) %>% 
  dplyr::group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assignment_prob = max(prob)) %>%  # Fixed typo: assingment -> assignment
  ungroup() %>%
  # Sort first by major cluster, then by assignment probability
  dplyr::arrange(likely_assignment, desc(assignment_prob)) %>%
  mutate(id = forcats::fct_inorder(as.factor(id)))

k3 %<>% left_join(metadata2,by = c("id" = "sample_id")) %>% distinct()

k3 <- k3 %>%
  arrange(pop.y, likely_assignment, desc(assignment_prob)) %>%
  mutate(id = forcats::fct_inorder(as.factor(id)))

b<-ggplot(k3 %>% filter(!is.na(pop.y)), aes(id, prob, fill = pop.x)) +
  geom_col(width = 2) +
  scale_fill_manual(values=c("skyblue1","#0072B2","#E69F00"))+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
  ylab("K3")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

b

####################################
#K4

k_4$id <- metadata$V1

k4 <- k_4 %>% 
  gather('pop', 'prob', V1:V4) %>% 
  dplyr::group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assignment_prob = max(prob)) %>%  # Fixed typo: assingment -> assignment
  ungroup() %>%
  # Sort first by major cluster, then by assignment probability
  dplyr::arrange(likely_assignment, desc(assignment_prob)) %>%
  mutate(id = forcats::fct_inorder(as.factor(id)))

k4 %<>% left_join(metadata2,by = c("id" = "sample_id")) %>% distinct()

k4 <- k4 %>%
  arrange(pop.y, likely_assignment, desc(assignment_prob)) %>%
  mutate(id = forcats::fct_inorder(as.factor(id)))


c<-ggplot(k4 %>% filter(!is.na(pop.y)), aes(id, prob, fill = pop.x)) +
  geom_col(width = 2) +
  scale_fill_manual(values=c("#8B8989","#E69F00","#0072B2","skyblue1"))+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
  ylab("K4")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

c

####################################
#K4

k_5$id <- metadata$V1

k5 <- k_5 %>% 
  gather('pop', 'prob', V1:V5) %>% 
  dplyr::group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assignment_prob = max(prob)) %>%  # Fixed typo: assingment -> assignment
  ungroup() %>%
  # Sort first by major cluster, then by assignment probability
  dplyr::arrange(likely_assignment, desc(assignment_prob)) %>%
  mutate(id = forcats::fct_inorder(as.factor(id)))

k5 %<>% left_join(metadata2,by = c("id" = "sample_id")) %>% distinct()

k5 <- k5 %>%
  arrange(pop.y, likely_assignment, desc(assignment_prob)) %>%
  mutate(id = forcats::fct_inorder(as.factor(id)))


d <-ggplot(k5 %>% filter(!is.na(pop.y)), aes(id, prob, fill = pop.x)) +
  geom_col(width = 2) +
  scale_fill_manual(values=c("#8B7B8B","#8B8989","#E69F00","#0072B2","skyblue1"))+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
  ylab("K4")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

d

#################################

a/b/c

ggsave("figures/06_admixture_allpop_lcwgs_pruned.png",height=15,width =27, unit = "cm", dpi = 600)


#################################
#Assign genetic cluster to each individual with best k = 4

individual_assignments <- k_4 %>%
  mutate(id = metadata$V1) %>%
  rowwise() %>%
  mutate(genetic_cluster = names(.)[which.max(c_across(V1:V4))],
         assignment_probability = max(c_across(V1:V4))) %>%
  ungroup() %>%
  left_join(metadata2, by = c("id" = "sample_id")) %>%
  dplyr::select(individual_id = id,
                pop = pop,
                genetic_cluster,
                assignment_probability) %>%
  distinct()

write.table(individual_assignments,"data/ngsadmix/admixture_assign_k4_allpop.txt",quote=F,row.names = F)








