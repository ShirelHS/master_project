# this script should compare the classification results of all datasets
setwd("C:/Users/ybnkb/Desktop/master degree/Hebrew Uni/rabani lab/new_pipeline_dataset/comparison")
rm(list = ls())

# install.packages("reshape")
library(reshape)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

# read classifications of few datasets 
medina_classified <- readRDS("../medina_dataset/exp_df_classified.RDS")
meyer_classified <- readRDS("../meier_dataset/exp_df_classified.RDS")
pauli_classified <- readRDS("../pauli_dataset/exp_df_classified.RDS")
zhao_classified <- readRDS("../zhao_dataset/exp_df_classified.RDS")

class_genes_inner <- zhao_classified %>%
  select(gene_id, gene_name, classification) %>%
  rename(c("class_zhao" = "classification")) %>%
  inner_join(meyer_classified %>% 
              select(gene_id, classification) %>% 
              rename(c("class_meyer" = "classification")), 
            by = c("gene_id" = "gene_id")) %>%
  inner_join(pauli_classified %>% 
              select(gene_id, classification) %>% 
              rename(c("class_pauli" = "classification")), 
            by = c("gene_id" = "gene_id")) %>%
  inner_join(medina_classified %>% 
              select(gene_id, classification) %>% 
              rename(c("class_medina" = "classification")), 
            by = c("gene_id" = "gene_id"))


class_genes_full <- zhao_classified %>%
  select(gene_id, gene_name, classification) %>%
  rename(c("class_zhao" = "classification")) %>%
  full_join(meyer_classified %>% 
               select(gene_id, classification) %>% 
               rename(c("class_meyer" = "classification")), 
             by = c("gene_id" = "gene_id")) %>%
  full_join(pauli_classified %>% 
               select(gene_id, classification) %>% 
               rename(c("class_pauli" = "classification")), 
             by = c("gene_id" = "gene_id")) %>%
  full_join(medina_classified %>% 
               select(gene_id, classification) %>% 
               rename(c("class_medina" = "classification")), 
             by = c("gene_id" = "gene_id"))
  

saveRDS(class_genes_inner, file = "classification_4_DS_inner.RDS")
saveRDS(class_genes_full, file = "classification_4_DS_full.RDS")

# for_heatmap_freq <- class_genes_inner %>% 
#   select(-gene_id) %>% 
#   pivot_longer(!gene_name) %>%
#   group_by(name, value) %>%
#   summarise(n = n()) %>%
#   data.frame()

same_combs <- class_genes_inner %>% 
  select(-gene_id) %>% 
  pivot_longer(cols = c(-gene_name),
               names_to = 'dataset',
               values_to = 'class') %>%
  group_by(dataset, class) %>%
  summarise(n = n()) %>%
  ungroup %>%
  mutate(class_1 = paste0(dataset, '_', class),
         class_2 = class_1) %>%
  select(n, class_1, class_2) %>%
  as.data.frame()

diff_combs <- class_genes_inner %>%
  count(class_zhao, class_meyer) %>%
  mutate(class_1 = paste0('class_zhao_', class_zhao)) %>%
  mutate(class_2 = paste0('class_meyer_', class_meyer)) %>%
  select(-class_zhao,-class_meyer) %>%
  rbind(class_genes_inner %>% 
          count(class_zhao, class_pauli) %>%
          mutate(class_1 = paste0('class_zhao_', class_zhao)) %>%
          mutate(class_2 = paste0('class_pauli_', class_pauli)) %>%
          select(-class_zhao,-class_pauli)) %>%
  rbind(class_genes_inner %>% 
          count(class_zhao, class_medina) %>%
          mutate(class_1 = paste0('class_zhao_', class_zhao)) %>%
          mutate(class_2 = paste0('class_medina_', class_medina)) %>%
          select(-class_zhao,-class_medina)) %>%
  rbind(class_genes_inner %>% 
          count(class_pauli, class_meyer) %>%
          mutate(class_1 = paste0('class_pauli_', class_pauli)) %>%
          mutate(class_2 = paste0('class_meyer_', class_meyer)) %>%
          select(-class_meyer,-class_pauli)) %>%
  rbind(class_genes_inner %>% 
        count(class_meyer, class_medina) %>%
        mutate(class_1 = paste0('class_meyer_', class_meyer)) %>%
        mutate(class_2 = paste0('class_medina_', class_medina)) %>%
        select(-class_meyer,-class_medina)) %>%
  rbind(class_genes_inner %>% 
          count(class_pauli, class_medina) %>%
          mutate(class_1 = paste0('class_pauli_', class_pauli)) %>%
          mutate(class_2 = paste0('class_medina_', class_medina)) %>%
          select(-class_pauli,-class_medina))


all_combs <- rbind(same_combs, diff_combs) %>%
  mutate(dataset_1 = gsub("_M|_MZ|_Z|_NONE|_NA", "", class_1)) %>%
  mutate(dataset_2 = gsub("_M|_MZ|_Z|_NONE|_NA", "", class_2)) 

all_combs_summ <- all_combs %>%
  mutate(dataset_1 = gsub("_M|_MZ|_Z|_NONE|_NA", "", class_1)) %>%
  mutate(dataset_2 = gsub("_M|_MZ|_Z|_NONE|_NA", "", class_2)) %>%
  group_by(dataset_1, dataset_2) %>%
  summarise(summ = sum(n))

all_combs_pc <- all_combs %>%
  full_join(all_combs_summ) %>%
  mutate(percent = round(n / summ * 100, 2))

ggplot(all_combs_pc, aes(class_1, class_2, fill = percent)) +
  geom_tile() +
  geom_text(aes(label = percent), color = 'white') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Common Classification Between each Two Datasets",
       subtitle = "percentage representation") +
  theme(text = element_text(size = 16)) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18, face = "bold"))



# roc curve each 2 datasets -----------------------------------------------

library("pROC")

# medina vs. meyer, while medina is the 'predicted' and meyer as the 'real'
 
combined_medina_meyer <- medina_classified %>% 
  rename(prediction_medina = classification) %>%
  inner_join(meyer_classified %>% 
               rename(real = classification) %>% 
               select(gene_id,real), 
             by = c("gene_id" = "gene_id"))

roc_m <- roc(response = factor(ifelse(combined_medina_meyer$real == "M", "M", "non-M"), ordered = TRUE), 
             predictor = factor(ifelse(combined_medina_meyer$prediction_medina == "M", "M", "non-M"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_z <- roc(response = factor(ifelse(combined_medina_meyer$real == "Z", "Z", "non-Z"), ordered = TRUE), 
             predictor = factor(ifelse(combined_medina_meyer$prediction_medina == "Z", "Z", "non-Z"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_mz <- roc(response = factor(ifelse(combined_medina_meyer$real == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
             predictor = factor(ifelse(combined_medina_meyer$prediction_medina == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

par(pty = "s")
plot(roc_m, col = "red", legacy.axes = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", 
     main = "ROC - medina against meyer (as the real class)", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.2, lwd = 3) 
plot(roc_mz, col = "green", add = TRUE, lwd = 3)
plot(roc_z, col = "blue", add = TRUE, lwd = 3)
legend("bottomright", c("maternal", "maternal-zygotic", "zygotic"), lty=1, 
       col = c("red", "green", "blue"), bty="n", inset=c(0,0.10), lwd = 3)
text(40,19.5, paste0("AUC: ", round(roc_m$auc,1), "%"), col = "red")
text(40,15.5, paste0("AUC: ", round(roc_mz$auc,1), "%"), col = "green")
text(40,11.5, paste0("AUC: ", round(roc_z$auc,1), "%"), col = "blue")



# zhao vs. medina, while zhao is the 'predicted' and medina as the 'real'

combined_zhao_medina <- zhao_classified %>% 
  rename(prediction_zhao = classification) %>%
  inner_join(medina_classified %>% 
               rename(real = classification) %>% 
               select(gene_id,real), 
             by = c("gene_id" = "gene_id"))

roc_m <- roc(response = factor(ifelse(combined_zhao_medina$real == "M", "M", "non-M"), ordered = TRUE), 
             predictor = factor(ifelse(combined_zhao_medina$prediction_zhao == "M", "M", "non-M"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_z <- roc(response = factor(ifelse(combined_zhao_medina$real == "Z", "Z", "non-Z"), ordered = TRUE), 
             predictor = factor(ifelse(combined_zhao_medina$prediction_zhao == "Z", "Z", "non-Z"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_mz <- roc(response = factor(ifelse(combined_zhao_medina$real == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
              predictor = factor(ifelse(combined_zhao_medina$prediction_zhao == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
              percent = TRUE, auc = TRUE)

par(pty = "s")
plot(roc_m, col = "red", legacy.axes = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", 
     main = "ROC - Zhao against Medina (as the real class)", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.2, lwd = 3) 
plot(roc_mz, col = "green", add = TRUE, lwd = 3)
plot(roc_z, col = "blue", add = TRUE, lwd = 3)
legend("bottomright", c("maternal", "maternal-zygotic", "zygotic"), lty=1, 
       col = c("red", "green", "blue"), bty="n", inset=c(0,0.10), lwd = 3)
text(40,19.5, paste0("AUC: ", round(roc_m$auc,1), "%"), col = "red")
text(40,15.5, paste0("AUC: ", round(roc_mz$auc,1), "%"), col = "green")
text(40,11.5, paste0("AUC: ", round(roc_z$auc,1), "%"), col = "blue")


# zhao vs. meyer, while zhao is the 'predicted' and meyer as the 'real'

combined_zhao_meyer <- zhao_classified %>% 
  rename(prediction_zhao = classification) %>%
  inner_join(meyer_classified %>% 
               rename(real = classification) %>% 
               select(gene_id,real), 
             by = c("gene_id" = "gene_id"))

roc_m <- roc(response = factor(ifelse(combined_zhao_meyer$real == "M", "M", "non-M"), ordered = TRUE), 
             predictor = factor(ifelse(combined_zhao_meyer$prediction_zhao == "M", "M", "non-M"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_z <- roc(response = factor(ifelse(combined_zhao_meyer$real == "Z", "Z", "non-Z"), ordered = TRUE), 
             predictor = factor(ifelse(combined_zhao_meyer$prediction_zhao == "Z", "Z", "non-Z"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_mz <- roc(response = factor(ifelse(combined_zhao_meyer$real == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
              predictor = factor(ifelse(combined_zhao_meyer$prediction_zhao == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
              percent = TRUE, auc = TRUE)

par(pty = "s")
plot(roc_m, col = "red", legacy.axes = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", 
     main = "ROC - Zhao against Meyer (as the real class)", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.2, lwd = 3) 
plot(roc_mz, col = "green", add = TRUE, lwd = 3)
plot(roc_z, col = "blue", add = TRUE, lwd = 3)
legend("bottomright", c("maternal", "maternal-zygotic", "zygotic"), lty=1, 
       col = c("red", "green", "blue"), bty="n", inset=c(0,0.10), lwd = 3)
text(40,19.5, paste0("AUC: ", round(roc_m$auc,1), "%"), col = "red")
text(40,15.5, paste0("AUC: ", round(roc_mz$auc,1), "%"), col = "green")
text(40,11.5, paste0("AUC: ", round(roc_z$auc,1), "%"), col = "blue")


# zhao vs. pauli, while zhao is the 'predicted' and pauli as the 'real'

combined_zhao_pauli <- zhao_classified %>% 
  rename(prediction_zhao = classification) %>%
  inner_join(pauli_classified %>% 
               rename(real = classification) %>% 
               select(gene_id,real), 
             by = c("gene_id" = "gene_id"))

roc_m <- roc(response = factor(ifelse(combined_zhao_pauli$real == "M", "M", "non-M"), ordered = TRUE), 
             predictor = factor(ifelse(combined_zhao_pauli$prediction_zhao == "M", "M", "non-M"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_z <- roc(response = factor(ifelse(combined_zhao_pauli$real == "Z", "Z", "non-Z"), ordered = TRUE), 
             predictor = factor(ifelse(combined_zhao_pauli$prediction_zhao == "Z", "Z", "non-Z"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_mz <- roc(response = factor(ifelse(combined_zhao_pauli$real == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
              predictor = factor(ifelse(combined_zhao_pauli$prediction_zhao == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
              percent = TRUE, auc = TRUE)

par(pty = "s")
plot(roc_m, col = "red", legacy.axes = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", 
     main = "ROC - Zhao against Pauli (as the real class)", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.2, lwd = 3) 
plot(roc_mz, col = "green", add = TRUE, lwd = 3)
plot(roc_z, col = "blue", add = TRUE, lwd = 3)
legend("bottomright", c("maternal", "maternal-zygotic", "zygotic"), lty=1, 
       col = c("red", "green", "blue"), bty="n", inset=c(0,0.10), lwd = 3)
text(40,19.5, paste0("AUC: ", round(roc_m$auc,1), "%"), col = "red")
text(40,15.5, paste0("AUC: ", round(roc_mz$auc,1), "%"), col = "green")
text(40,11.5, paste0("AUC: ", round(roc_z$auc,1), "%"), col = "blue")


# plot roc againt Lior's classifications ---------------------------------------------------------------

# zhao against Liors class
liors_classes <- read.delim("class.5.txt",header = FALSE)
colnames(liors_classes) <- c("gene_name", "lior_class")
zhao_against_lior <- zhao_classified %>%
  inner_join(liors_classes %>% mutate(lior_class = toupper(lior_class)), 
             by = c("gene_name" = "gene_name"))

roc_m <- roc(response = factor(ifelse(zhao_against_lior$lior_class == "M", "M", "non-M"), ordered = TRUE), 
             predictor = factor(ifelse(zhao_against_lior$classification == "M", "M", "non-M"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_z <- roc(response = factor(ifelse(zhao_against_lior$lior_class == "Z", "Z", "non-Z"), ordered = TRUE), 
             predictor = factor(ifelse(zhao_against_lior$classification == "Z", "Z", "non-Z"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

# roc_mz <- roc(response = factor(ifelse(zhao_against_lior$lior_class == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
#               predictor = factor(ifelse(zhao_against_lior$classification == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
#               percent = TRUE, auc = TRUE)

par(pty = "s")
plot(roc_m, col = "red", legacy.axes = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", 
     main = "ROC - Zhao against Liors classes (as the real class)", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.2, lwd = 3) 
plot(roc_z, col = "blue", add = TRUE, lwd = 3)
legend("bottomright", c("maternal", "zygotic"), lty=1, 
       col = c("red", "blue"), bty="n", inset=c(0,0.10), lwd = 3)
text(30,15.5, paste0("AUC: ", round(roc_m$auc,1), "%"), col = "red")
text(30,11.5, paste0("AUC: ", round(roc_z$auc,1), "%"), col = "blue")


# pauli against Liors class
pauli_against_lior <- pauli_classified %>%
  inner_join(liors_classes %>% mutate(lior_class = toupper(lior_class)), 
             by = c("gene_name" = "gene_name"))

roc_m <- roc(response = factor(ifelse(pauli_against_lior$lior_class == "M", "M", "non-M"), ordered = TRUE), 
             predictor = factor(ifelse(pauli_against_lior$classification == "M", "M", "non-M"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_z <- roc(response = factor(ifelse(pauli_against_lior$lior_class == "Z", "Z", "non-Z"), ordered = TRUE), 
             predictor = factor(ifelse(pauli_against_lior$classification == "Z", "Z", "non-Z"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

# roc_mz <- roc(response = factor(ifelse(pauli_against_lior$lior_class == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
#               predictor = factor(ifelse(pauli_against_lior$classification == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
#               percent = TRUE, auc = TRUE)

par(pty = "s")
plot(roc_m, col = "red", legacy.axes = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", 
     main = "ROC - Pauli against Liors classes (as the real class)", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.2, lwd = 3) 
plot(roc_z, col = "blue", add = TRUE, lwd = 3)
legend("bottomright", c("maternal", "zygotic"), lty=1, 
       col = c("red", "blue"), bty="n", inset=c(0,0.10), lwd = 3)
text(30,15.5, paste0("AUC: ", round(roc_m$auc,1), "%"), col = "red")
text(30,11.5, paste0("AUC: ", round(roc_z$auc,1), "%"), col = "blue")

# meyer against Liors class
meyer_against_lior <- meyer_classified %>%
  inner_join(liors_classes %>% mutate(lior_class = toupper(lior_class)), 
             by = c("gene_name" = "gene_name"))

roc_m <- roc(response = factor(ifelse(meyer_against_lior$lior_class == "M", "M", "non-M"), ordered = TRUE), 
             predictor = factor(ifelse(meyer_against_lior$classification == "M", "M", "non-M"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_z <- roc(response = factor(ifelse(meyer_against_lior$lior_class == "Z", "Z", "non-Z"), ordered = TRUE), 
             predictor = factor(ifelse(meyer_against_lior$classification == "Z", "Z", "non-Z"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

par(pty = "s")
plot(roc_m, col = "red", legacy.axes = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", 
     main = "ROC - Meyer against Liors classes (as the real class)", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.2, lwd = 3) 
plot(roc_z, col = "blue", add = TRUE, lwd = 3)
legend("bottomright", c("maternal", "zygotic"), lty=1, 
       col = c("red", "blue"), bty="n", inset=c(0,0.10), lwd = 3)
text(30,15.5, paste0("AUC: ", round(roc_m$auc,1), "%"), col = "red")
text(30,11.5, paste0("AUC: ", round(roc_z$auc,1), "%"), col = "blue")

# meyer against Liors class
medina_against_lior <- medina_classified %>%
  inner_join(liors_classes %>% mutate(lior_class = toupper(lior_class)), 
             by = c("gene_name" = "gene_name"))

roc_m <- roc(response = factor(ifelse(medina_against_lior$lior_class == "M", "M", "non-M"), ordered = TRUE), 
             predictor = factor(ifelse(medina_against_lior$classification == "M", "M", "non-M"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_z <- roc(response = factor(ifelse(medina_against_lior$lior_class == "Z", "Z", "non-Z"), ordered = TRUE), 
             predictor = factor(ifelse(medina_against_lior$classification == "Z", "Z", "non-Z"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

par(pty = "s")
plot(roc_m, col = "red", legacy.axes = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", 
     main = "ROC - Medina against Liors classes (as the real class)", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.2, lwd = 3) 
plot(roc_z, col = "blue", add = TRUE, lwd = 3)
legend("bottomright", c("maternal", "zygotic"), lty=1, 
       col = c("red", "blue"), bty="n", inset=c(0,0.10), lwd = 3)
text(30,15.5, paste0("AUC: ", round(roc_m$auc,1), "%"), col = "red")
text(30,11.5, paste0("AUC: ", round(roc_z$auc,1), "%"), col = "blue")



# roc curve against subset of known classificated genes -------------------


# plot roc against a subset of known genes
# a known classificated genes 
# Lior gave this dataset to me
maternal <- c("buc", "dazl", "eomesa", "cldng", "cth1", "gyg1a", "idhba", "yipf5", "cldnd", "dnd1", "efhd2", "hmgb3b", "nanos3", "org", "pabpn1l")
zygotic <- c("tbxtl", "ntl", "eve1", "sp5l", "apoc2", "sox2", "bmp4", "her1", "gadd45bb", "bmp2b", "foxa3", "her5", "klf2b", "sox32", "ywhae2")

idsM = c('dazl', 'zar1l', 'zar1', 'wee2', 'tdrd5', 'tdrd1', 'tdrd7a', 'slbp2', 'buc2l', 'btg4', 'fbxo43', 'mos', 'bmp15',
  'dab2', 'h1m','cldnd', 'prkcda' ,'lrrc3', 'rubcn', 'stx11a', 'rheb', 'st3gal1', 'trmt9b')
idsMZ = c('pou5f3', 'sox11b', 'magoh', 'ddx52', 'ddx27' ,'ddx18', 'ddx21', 'ddx39aa', 'ddx19', 'ddx3b',
          'ddx5', 'ddx39ab', 'hist1h', 'hist2h', 'h3f3')
idsZ = c('tbxta', 'tbx16', 'ved', 'vox', 'sox19a', 'sox3', 'sox32', 'tfec', 'klf17', 'klf2b', 'foxd5', 'chrd', 'lft2',
         'gadd45bb', 'gadd45ga', 'gadd45ba', 'histh1l', 'znf1145', 'mxtx2', 'skilb' ,'slc26a1', 'gata3', 'hnf4a')
# idsX = c('nanog', 'cth1', 'tbpl2', 'adcy9', 'bach2a')

maternal_list <- tolower(read.delim("expression.M.txt", header = FALSE)[,1])
zygotic_list <- tolower(read.delim("expression.Z.txt", header = FALSE)[,1])

M_subset <- data.frame(gene_name = unique(c(maternal, idsM, maternal_list)), class_real = "M")
Z_subset <- data.frame(gene_name = unique(c(zygotic, idsZ, zygotic_list)), class_real = "Z")
MZ_subset <- data.frame(gene_name = unique(idsMZ), class_real = "MZ")

all_subset_2gether <- rbind(M_subset, Z_subset, MZ_subset)


# the graph against zhao
combined_zhao_subset <- zhao_classified %>% 
  rename(prediction_zhao = classification) %>%
  inner_join(all_subset_2gether,
             by = c("gene_name" = "gene_name"))

roc_m <- roc(response = factor(ifelse(combined_zhao_subset$class_real == "M", "M", "non-M"), ordered = TRUE), 
             predictor = factor(ifelse(combined_zhao_subset$prediction_zhao == "M", "M", "non-M"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_z <- roc(response = factor(ifelse(combined_zhao_subset$class_real == "Z", "Z", "non-Z"), ordered = TRUE), 
             predictor = factor(ifelse(combined_zhao_subset$prediction_zhao == "Z", "Z", "non-Z"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_mz <- roc(response = factor(ifelse(combined_zhao_subset$class_real == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
              predictor = factor(ifelse(combined_zhao_subset$prediction_zhao == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
              percent = TRUE, auc = TRUE)

par(pty = "s")
plot(roc_m, col = "red", legacy.axes = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", 
     main = "ROC - Zhao Against Subset of Known Classificated genes", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.2, lwd = 3) 
plot(roc_mz, col = "green", add = TRUE, lwd = 3)
plot(roc_z, col = "blue", add = TRUE, lwd = 3)
legend("bottomright", c("maternal", "maternal-zygotic", "zygotic"), lty=1, 
       col = c("red", "green", "blue"), bty="n", inset=c(0,0.10), lwd = 3)
text(40,19.5, paste0("AUC: ", round(roc_m$auc,1), "%"), col = "red")
text(40,15.5, paste0("AUC: ", round(roc_mz$auc,1), "%"), col = "green")
text(40,11.5, paste0("AUC: ", round(roc_z$auc,1), "%"), col = "blue")


# the graph against pauli
combined_pauli_subset <- pauli_classified %>% 
  rename(prediction_pauli = classification) %>%
  inner_join(all_subset_2gether,
             by = c("gene_name" = "gene_name"))

roc_m <- roc(response = factor(ifelse(combined_pauli_subset$class_real == "M", "M", "non-M"), ordered = TRUE), 
             predictor = factor(ifelse(combined_pauli_subset$prediction_pauli == "M", "M", "non-M"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_z <- roc(response = factor(ifelse(combined_pauli_subset$class_real == "Z", "Z", "non-Z"), ordered = TRUE), 
             predictor = factor(ifelse(combined_pauli_subset$prediction_pauli == "Z", "Z", "non-Z"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_mz <- roc(response = factor(ifelse(combined_pauli_subset$class_real == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
              predictor = factor(ifelse(combined_pauli_subset$prediction_pauli == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
              percent = TRUE, auc = TRUE)

par(pty = "s")
plot(roc_m, col = "red", legacy.axes = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", 
     main = "ROC - Pauli Against Subset of Known Classificated genes", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.2, lwd = 3) 
plot(roc_mz, col = "green", add = TRUE, lwd = 3)
plot(roc_z, col = "blue", add = TRUE, lwd = 3)
legend("bottomright", c("maternal", "maternal-zygotic", "zygotic"), lty=1, 
       col = c("red", "green", "blue"), bty="n", inset=c(0,0.10), lwd = 3)
text(40,19.5, paste0("AUC: ", round(roc_m$auc,1), "%"), col = "red")
text(40,15.5, paste0("AUC: ", round(roc_mz$auc,1), "%"), col = "green")
text(40,11.5, paste0("AUC: ", round(roc_z$auc,1), "%"), col = "blue")


# the graph against meyer
combined_meyer_subset <- meyer_classified %>% 
  rename(prediction_meyer = classification) %>%
  inner_join(all_subset_2gether,
             by = c("gene_name" = "gene_name"))

roc_m <- roc(response = factor(ifelse(combined_meyer_subset$class_real == "M", "M", "non-M"), ordered = TRUE), 
             predictor = factor(ifelse(combined_meyer_subset$prediction_meyer == "M", "M", "non-M"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_z <- roc(response = factor(ifelse(combined_meyer_subset$class_real == "Z", "Z", "non-Z"), ordered = TRUE), 
             predictor = factor(ifelse(combined_meyer_subset$prediction_meyer == "Z", "Z", "non-Z"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_mz <- roc(response = factor(ifelse(combined_meyer_subset$class_real == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
              predictor = factor(ifelse(combined_meyer_subset$prediction_meyer == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
              percent = TRUE, auc = TRUE)

par(pty = "s")
plot(roc_m, col = "red", legacy.axes = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", 
     main = "ROC - Meyer Against Subset of Known Classificated genes", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.2, lwd = 3) 
plot(roc_mz, col = "green", add = TRUE, lwd = 3)
plot(roc_z, col = "blue", add = TRUE, lwd = 3)
legend("bottomright", c("maternal", "maternal-zygotic", "zygotic"), lty=1, 
       col = c("red", "green", "blue"), bty="n", inset=c(0,0.10), lwd = 3)
text(40,19.5, paste0("AUC: ", round(roc_m$auc,1), "%"), col = "red")
text(40,15.5, paste0("AUC: ", round(roc_mz$auc,1), "%"), col = "green")
text(40,11.5, paste0("AUC: ", round(roc_z$auc,1), "%"), col = "blue")


# the graph against meyer
combined_medina_subset <- medina_classified %>% 
  rename(prediction_medina = classification) %>%
  inner_join(all_subset_2gether,
             by = c("gene_name" = "gene_name"))

roc_m <- roc(response = factor(ifelse(combined_medina_subset$class_real == "M", "M", "non-M"), ordered = TRUE), 
             predictor = factor(ifelse(combined_medina_subset$prediction_medina == "M", "M", "non-M"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_z <- roc(response = factor(ifelse(combined_medina_subset$class_real == "Z", "Z", "non-Z"), ordered = TRUE), 
             predictor = factor(ifelse(combined_medina_subset$prediction_medina == "Z", "Z", "non-Z"), ordered = TRUE), 
             percent = TRUE, auc = TRUE)

roc_mz <- roc(response = factor(ifelse(combined_medina_subset$class_real == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
              predictor = factor(ifelse(combined_medina_subset$prediction_medina == "MZ", "MZ", "non-MZ"), ordered = TRUE), 
              percent = TRUE, auc = TRUE)

par(pty = "s")
plot(roc_m, col = "red", legacy.axes = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", 
     main = "ROC - Medina Against Subset of Known Classificated genes", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.2, lwd = 3) 
plot(roc_mz, col = "green", add = TRUE, lwd = 3)
plot(roc_z, col = "blue", add = TRUE, lwd = 3)
legend("bottomright", c("maternal", "maternal-zygotic", "zygotic"), lty=1, 
       col = c("red", "green", "blue"), bty="n", inset=c(0,0.10), lwd = 3)
text(40,19.5, paste0("AUC: ", round(roc_m$auc,1), "%"), col = "red")
text(40,15.5, paste0("AUC: ", round(roc_mz$auc,1), "%"), col = "green")
text(40,11.5, paste0("AUC: ", round(roc_z$auc,1), "%"), col = "blue")




















# --------------------------------------------------------------------------------------------

# calculate pearson correlation of classifications between zhao to the other datasets.
# no such thing!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# --------------------------------------------------------------------------------------------

# zhao_n <- class_genes_full %>%
#   group_by(class_zhao) %>%
#   summarise(n = n()) %>%
#   rename(class = class_zhao) %>%
#   rename(n_zhao = n) %>%
#   as.data.frame()
# 
# meyer_n <- class_genes_full %>%
#   group_by(class_meyer) %>%
#   summarise(n = n()) %>%
#   rename(class = class_meyer) %>%
#   rename(n_meyer = n) %>%
#   as.data.frame()
# 
# zhao_vs_meyer_comp <- zhao_n %>%
#   left_join(meyer_n) %>%
#   select(n_zhao, n_meyer) %>%
#   mutate(summed = rowSums(.)) %>%
#   add_column(class = zhao_n$class)
# 
# sum_col_zhao <- sum(zhao_vs_meyer_comp$n_zhao)
# sum_col_meyer <- sum(zhao_vs_meyer_comp$n_meyer)
# 
