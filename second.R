# this script is:
# do manual classifications to Maternal/Zygotic/Maternal-Zygotic: 
# first take the three criteria for classification
# second plot the distribution of these three criteria
# third classify genes

# the data needed:
# log_exp_tbl.csv file

setwd(getwd()) # set the directory to work in
rm(list = ls()) # delete/clear global environment

# ---------------------------------------------------------------------------------
# 1.
# this part takes the MZT hour of this organism
# read the filtered normalized data from csv file
# cut the data to mature and precursor 

file_name <- "normalized_exp.csv"

read_relev_data <- function(file_name) {
  
  normalized_exp <- read.csv(file = file_name, header = TRUE, sep = ",")[,-1]
  return(normalized_exp)
  
}

get_mzt <- function() {
  
  mzt <- readline(prompt = "Enter hour of the MZT: ")
  return(mzt)
}


normalized_exp <- read_relev_data("normalized_exp.csv")

# for zebrafish:
# MZT <- 3.5 
MZT <- as.numeric(get_mzt())
  


# take the mature transcripts
mature <- cbind(normalized_exp[,1:2], normalized_exp[,grep("transcript", colnames(normalized_exp))])

# take the pre-mRNA
precursor <- cbind(normalized_exp[,1:2], normalized_exp[,grep("precursor", colnames(normalized_exp))])


# ----------------------------------------------------------------------------------
# 2.
# this part is taking the three criteria 

# ----------------
# first criterion:
# take maximum hour for mature transcript before the mzt:
gene_hours <- as.numeric(gsub("transcript_|_\\d+","",colnames(mature)[-c(1,2)]))
hours_bf_mzt <- gene_hours[gene_hours < MZT]
before_mzt_col_name <- mature[,1:(length(hours_bf_mzt)+2)] 
before_mzt_col_name$max_ <- unlist(apply(before_mzt_col_name[,-c(1,2)], 1, max))


# -----------------
# second criterion:
# take the ratio of maximum expression of pre-mRNA after MZT to before MZT
pre_hours <- as.numeric(gsub("precursor_|_\\d+","",colnames(precursor)[-c(1,2)]))
pre_hour_names_aft_mzt <- pre_hours[pre_hours > MZT]
pre_hour_names_bf_mzt <- pre_hours[pre_hours < MZT]

pre_hour_bf_mzt <- precursor[,1:(length(pre_hour_names_bf_mzt)+2)]
pre_hour_aft_mzt <- precursor[,c(1,2,(length(pre_hour_names_bf_mzt)+3):(ncol(precursor)))]
pre_hour_bf_mzt$max_ <- unlist(apply(pre_hour_bf_mzt[,-c(1,2)], 1, max))
pre_hour_aft_mzt$max_ <- unlist(apply(pre_hour_aft_mzt[,-c(1,2)], 1, max))

# becuase the maximum expression of precursor before the mzt is the deliminator, it shouldn't be zero 
if (any(pre_hour_bf_mzt$max_ == 0))
  pre_hour_bf_mzt$max_[which(pre_hour_bf_mzt$max_ ==0)] <- 0.1^5

ratio_pre <- cbind(pre_hour_bf_mzt[,1:2], ratio_precursor = pre_hour_aft_mzt$max_ - pre_hour_bf_mzt$max_)


# ----------------
# third criterion:
# ratio of mature log max hour after mzt to mature log max hour before mzt
mature_hour_aft_mzt <- mature[,c(1,2,(length(hours_bf_mzt)+3):(ncol(mature)))]
mature_hour_aft_mzt$max_ <- unlist(apply(mature_hour_aft_mzt[,-c(1,2)], 1, max))

# becuase the maximum expression of mature before the mzt is the deliminator, it shouldn't be zero 
if (any(before_mzt_col_name$max_ == 0))
  before_mzt_col_name$max_[which(before_mzt_col_name$max_ ==0)] <- 0.1^5

ratio_mature <- cbind(mature_hour_aft_mzt[,1:2], ratio_mature = mature_hour_aft_mzt$max_ - before_mzt_col_name$max_)


# ------------------------------------------------------------------------------------------
# 3.
# combine all criteria to one data frame
# determine threshold for the 3 criteria

three_df <- cbind(ratio_mature[,1:2],
                  mature_bf_mzt = before_mzt_col_name$max_ , 
                  ratio_precursor = ratio_pre$ratio_precursor, 
                  ratio_mature = ratio_mature$ratio_mature)


threshold_mat_bf_mzt <- (mean(three_df$mature_bf_mzt) - 1*sd(three_df$mature_bf_mzt))
threshold_ratio_precursor <- log2(1.25)
threshold_ratio_mature <- log2(0.9)

# -------------------------------------------------------------------------------------------
# 4.
# create histograms of each one from the three criteria in PDF file

pdf("three_criteria_graph.pdf")

# maximum hour for mature transcript before the mzt from section 1.
ggplot(data = three_df, aes(x = mature_bf_mzt)) + 
  geom_histogram(binwidth = 0.1, color="plum4", fill="plum4")  + 
  coord_cartesian(ylim = c(0,700)) + 
  geom_vline(aes(xintercept = threshold_mat_bf_mzt), col = "orange", show.legend = F) + 
  labs(title = "Maximum Exp for Mature Transcripts Before MZT", x = "log(FPKM)") + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 18), legend.text = element_text(size = 18))

# ratio precursor after the mzt to before max exp from section 2.
ggplot(data = three_df, aes(x = ratio_precursor)) + 
  geom_histogram(binwidth = 0.03, color="plum4", fill="plum4")  + 
  coord_cartesian(ylim = c(0,700)) + 
  geom_vline(aes(xintercept = threshold_ratio_precursor), col = "orange", show.legend = F) + 
  labs(title = "Ratio of Precursor After MZT to Before MZT", x = "log(FPKM) ratio") + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 18), legend.text = element_text(size = 18))

# ratio of mature max exp after mzt to before mzt from section 3.
ggplot(data = three_df, aes(x = ratio_mature)) + 
  geom_histogram(binwidth = 0.1, color="plum4", fill="plum4")  + 
  coord_cartesian(xlim = c(- 10, 10), ylim = c(0,700)) + 
  geom_vline(aes(xintercept = threshold_ratio_mature), col = "orange", show.legend = F) + 
  labs(title = "Ratio of Mature Transcript After MZT to Before MZT", x = "log(FPKM) ratio") + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 18), legend.text = element_text(size = 18))

dev.off()


# ----------------------------------------------------------------------------------------
# 5.
# caracterize the genes by the three criteria and three thresholds
# classify by its characterization

three_df$class_mat_bf <- "-"
three_df$class_mat_bf[which(three_df$mature_bf_mzt >= threshold_mat_bf_mzt)] <- "+"

three_df$class_ratio_pre <- "-"
three_df$class_ratio_pre[which(three_df$ratio_precursor >= threshold_ratio_precursor)] <- "+"

three_df$class_ratio_mat <- "-"
three_df$class_ratio_mat[which(three_df$ratio_mature >= threshold_ratio_mature)] <- "+"

three_df$classification <- "NONE"
# for maternal class is should:
# have high expression before the MZT (+)
# have low or no change of precursor expression (-)
# have negative change of mature expression, because it suppose to degrade after the MZT (-)
three_df$classification[(three_df$class_mat_bf == "+") & 
                          (three_df$class_ratio_pre == "-") & 
                          (three_df$class_ratio_mat == "-")] <- "M"


# for zygotic class is should:
# have low to no expression before the MZT (-)
# have high change of precursor expression (+)
# have high positive change of mature expression, because it suppose to express only after the MZT (+)
three_df$classification[(three_df$class_mat_bf == "-") &
                          (three_df$class_ratio_pre == "+") &
                          (three_df$class_ratio_mat == "+")] <- "Z"


# for maternal-zygotic class is should:
# have high expression before the MZT (+) 
# have high change of precursor expression (+)
# have high positive (or not) change of mature expression, because it suppose to express only after the MZT (+/-)
three_df$classification[(three_df$class_mat_bf == "+") &
                          (three_df$class_ratio_pre == "+") &
                          (three_df$class_ratio_mat == "+")] <- "MZ"

three_df$classification[(three_df$class_mat_bf == "+") &
                          (three_df$class_ratio_pre == "+") &
                          (three_df$class_ratio_mat == "-")] <- "MZ"

saveRDS(three_df, file = "class_criteria.RDS")

# -----------------------------------------------------------------------------------------
# 6.
# plot a graph of classification division


ggplot(three_df, aes(classification)) + 
  geom_bar() + 
  labs(title = "manual classification - counts of each category") + 
  theme(text = element_text(size = 18))

# save the expression data frame with the classification
exp_classified <- normalized_exp %>%
  inner_join(three_df %>%
               select(gene_id, gene_name, classification), by = c("gene_id", "gene_name"))

saveRDS(exp_classified, file = "exp_df_classified.RDS")

write_csv(exp_classified, "exp_df_classified.csv")

