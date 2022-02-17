# this script is:
# calculate linear regression model for the maternals only:
# first filter for the classified maternal genes
# second calculate mathematical models for the degradation
# third choose the better model with likelihood ratio test

# the data needed:
# class_criteria.RDS
# log_exp_tbl.csv file

setwd(getwd())
rm(list = ls())

# load data ----------------------------------------------------------------------------------
# 1.
# this part take the classifications and the exp levels and 
# filter for mature transcript of maternal genes 

class_criteria <- readRDS(file = "class_criteria.RDS")

normalized_exp <- read.csv(file = "normalized_exp.csv", header = TRUE, sep = ",")[,-1]

maternal_only <- normalized_exp %>%
  inner_join(class_criteria %>% 
               select(gene_id, classification)) %>%
  filter(classification == "M") %>%
  select(gene_id, gene_name, starts_with("transcript_"))


# model 0 data --------------------------------------------------------------------------------- 
# 2.
# prepering data for M0 - constant model

# take only exp levels for linear fit
m0_matrix <- maternal_only %>% 
  select(starts_with("transcript_"))

# move the name and id to row names
gene_id_name <- unlist(sapply(1:nrow(maternal_only), function(x) paste0(maternal_only$gene_id[x], "_", maternal_only$gene_name[x])))
row.names(m0_matrix) <- gene_id_name

# change colnames to hours
colnames(m0_matrix) <- gsub("transcript_|_\\d+", "", colnames(m0_matrix))

m_not_numeric <- apply(m0_matrix,1, stack)

# this function is turning the hours, for stacked gene, into numeric
hours_2_numeric <- function(gene){
  gene$ind <- as.numeric(levels(gene$ind))[gene$ind]
  return(gene)
}
m_data <- lapply(m_not_numeric, hours_2_numeric)

# change hours to 0 and create list of stacks for zeros
m0_zeros <- m0_matrix
colnames(m0_zeros) <- rep(0, times = ncol(m0_zeros))
m0_data_zeros <- apply(m0_zeros,1, stack)
m0_data <- lapply(m0_data_zeros, hours_2_numeric)

# model 1 data -----------------------------------------------------------------------
# 3.
# prepering data for M1 - simple linear model

# expression per hour w/o changes for early degradation onset
m1_data <- m_data

# model 2 data -----------------------------------------------------------------------
# 4.
# find the best onset for each hour

# this function is getting a gene with expression and hours
# it calculates linear models for each start degradation time 
# it computes r square for each model to determine the best t1
# it returns the gene with hours manipulated according to t1
find_t1 <- function(gene, to_hour) {
  
  best_t1 <- choose_t1(gene, to_hour)
  gene$ind[gene$ind < best_t1] <- best_t1
  return(gene)
  
}

# calculate all linear models for changing t1 (start checking from 0)
# take a gene with expression and hours
# returns the t1 for which the r square was the maximum
choose_t1 <- function(gene, to_hour) {
  
  r_square <- 0
  t1 <- 0
  for (h in seq(0,to_hour, 0.25)) { # gene$ind[nrow(gene)]
    gene_h <- gene
    gene_h$ind[gene_h$ind < h] <- h
    lm_h <- calculate_lm(gene_h)
    if (is.null(lm_h)) 
      next
    if (r_square < summary(lm_h[[1]])$r.squared) {
      r_square <- summary(lm_h[[1]])$r.squared
      t1 <- h
    }
  }
  
  return(t1)
  
}

latest_t1 <- unique(m1_data[[1]]$ind)[(length(unique(m1_data[[1]]$ind)))-2]

# for each gene- calculate the t1
m2_data <- lapply(m1_data, function(x) find_t1(x, latest_t1))


# model 3 data -----------------------------------------------------------------------
# 5.
# prepering data for M2 - complecated linear model
# first find the best onset for each hour
# second manipualate hours according to t1

# calculate the linear model for a gene with expression per hours
calculate_lm <- function(gene) {
  tmp_lm <- lm(values ~ ind, gene)
  
  # if the slope is greater than 0 it returns NULL
  if ((is.na(tmp_lm)) || (is.null(tmp_lm)) || 
      (!is.na(tmp_lm$coefficients[["ind"]]) && tmp_lm$coefficients[["ind"]] > 0)) {
    return(NULL)
  }
  return(list(tmp_lm,tmp_lm$fitted.values))
}

# get a gene with expression and hours
# calculate linear models for each start degradation time 
# compute r square for each model to determine the best t1
# return the gene with hours manipulated according to t1
find_t1_t2 <- function(gene, to_hour) {
  
  best_t <- choose_t1_t2(gene, last_hour = to_hour)
  gene$ind[gene$ind < best_t[[1]]] <- best_t[[1]]
  gene$ind[gene$ind > best_t[[2]]] <- best_t[[2]]
  
  return(gene)
  
}

# calculate all linear models for changing t1 (start checking from 0)
# take a gene with expression and hours
# returns the t1 for which the r square was the maximum
choose_t1_t2 <- function(gene, last_hour) {
  
  r_square <- 0
  t1 <- 0
  t2 <- gene$ind[(length(gene$ind))]
  
  # the sequence of samples from the first sample to the last_hour sample
  gene_hours_uniq <- unique(gene$ind)[1:which(unique(gene$ind) == last_hour)]
  
  for (h1 in gene_hours_uniq) { # gene$ind[nrow(gene)]
    gene_h <- gene
    gene_h$ind[gene_h$ind < h1] <- h1
    
    # the sequence of samples from h2 sample to the last sample (12 in Zhao et al. dataset)
    genes_h2_seq <- unique(gene$ind)[(which(unique(gene$ind) == h1)+2):length(unique(gene$ind))]
    # if this is the last sample, don't take na as an sample 
    if (h1 == gene_hours_uniq[length(gene_hours_uniq)]) {
      genes_h2_seq <- unique(gene$ind)[length(unique(gene$ind))]
    }
    
    # calculate t2 for the last few samples
    for (h2 in genes_h2_seq) {
      gene_h <- gene
      gene_h$ind[gene_h$ind < h1] <- h1
      gene_h$ind[gene_h$ind > h2] <- h2
      
      lm_h <- calculate_lm(gene_h)
      if (is.null(lm_h)) 
        next
      if (r_square < summary(lm_h[[1]])$r.squared) {
        r_square <- summary(lm_h[[1]])$r.squared
        t1 <- h1
        t2 <- h2
      }
      
    }
    
  }
  return(list(t1, t2))
}

# limit onset time to 3 hours before the last timepoint
latest_sample_hour <- unique(m1_data[[1]]$ind)[(length(unique(m1_data[[1]]$ind)))-2]

# for each gene- calculate the t1
m3_data <- lapply(m1_data, function(x) find_t1_t2(x, latest_sample_hour))


# save data ----------------------------------------------------------------------
# 5. 
# save data for each gene

# save data for models
save(m0_data, file = "data_M0.rdata")
save(m1_data, file = "data_M1.rdata")
save(m2_data, file = "data_M2.rdata")
save(m3_data, file = "data_M3.rdata")

# compute linear models ----------------------------------------------------------
# 6.
# compute nested linear models to represent the degradation rates
# option: flatten last few samples to permenant expression level 

model_M0 <- lapply(m0_data, calculate_lm)
model_M1 <- lapply(m1_data, calculate_lm)
model_M2 <- lapply(m2_data, calculate_lm)
model_M3 <- lapply(m3_data, calculate_lm)

# remove from the list genes that may have a slope > 0 from models lists and from data lists
# take the names of the genes that has a model null because of slope > 0
null_model_1 <- names(m2_data[unlist(lapply(1:length(model_M1), function(i) if(is.null(model_M1[[i]])) return(i)))])
null_model_2 <- names(m2_data[unlist(lapply(1:length(model_M2), function(i) if(is.null(model_M2[[i]])) return(i)))])
null_model_3 <- names(m3_data[unlist(lapply(1:length(model_M3), function(i) if(is.null(model_M3[[i]])) return(i)))])

# intersect the names and remove from the models list and the data
gene_names_to_drop <- union(null_model_1, null_model_2)
model_M0 <- model_M0[names(model_M0) %in% gene_names_to_drop == FALSE]
model_M1 <- model_M1[names(model_M1) %in% gene_names_to_drop == FALSE]
model_M2 <- model_M2[names(model_M2) %in% gene_names_to_drop == FALSE]
model_M3 <- model_M3[names(model_M3) %in% gene_names_to_drop == FALSE]

m0_data <- m0_data[names(m0_data) %in% gene_names_to_drop == FALSE]
m1_data <- m1_data[names(m1_data) %in% gene_names_to_drop == FALSE]
m2_data <- m2_data[names(m2_data) %in% gene_names_to_drop == FALSE]
m3_data <- m3_data[names(m3_data) %in% gene_names_to_drop == FALSE]

save(model_M0, file = "model_M0.rdata")
save(model_M1, file = "model_M1.rdata")
save(model_M2, file = "model_M2.rdata")
save(model_M3, file = "model_M3.rdata")

# save data for models
save(m0_data, file = "data_M0_wo_slope_gt_0.rdata")
save(m1_data, file = "data_M1_wo_slope_gt_0.rdata")
save(m2_data, file = "data_M2_wo_slope_gt_0.rdata")
save(m3_data, file = "data_M3_wo_slope_gt_0.rdata")


# wilks test ---------------------------------------------------------------------
# 7.
# model examining with Wilks test (log likelihood ratio with chi square test) for the 3 nested models

# calculate std for datasets with more than 1 repetition
# std for specific gene and given hour
std_specific_hour <- function(gene, h) {
  indices <- which(gene$ind %in% h)
  return(sqrt((sum((abs(gene$values[indices] - mean(gene$values[indices])))^2) / length(indices))))
  
}

# likelihood calculation - MSE
# calculate the std for each hour and compute the likelihood (MSE)
calculate_MSE <- function(gene, fitted_vals) {
  
  std_whole_gene <- unlist(sapply(gene$ind, function(x) std_specific_hour(gene,x)))
  std_whole_gene[std_whole_gene == 0] <- (0.1)^10
  mse_val <- sum(((abs(gene$values - fitted_vals))^2) / std_whole_gene^2)
  return(mse_val)
  
}

# the statistics
calculate_LR <- function(gene, model_x, model_y) {
  # model_x[[2]] are the fitted values
  return((calculate_MSE(gene, model_x[[2]]) - calculate_MSE(gene, model_y[[2]])))
  
}


# chi square test ----------------------------------------------------------------
# 8.
# model examining with Wilks test (log likelihood ratio with chi square test) for the 3 nested models


# choose the right model from M0 and M1
# and return the value of calculate_LR with the best model 
send_right_model <- function(gene, first_match, model0, model1, model2) {
  
  tmp_model <- if (first_match == "M0") model0 else model1
  return(calculate_LR(gene, tmp_model, model2))
}


# input p value and names of two models
# output the name of the accepted model
choose_better_model <- function(p_val, m_x,m_y) {
  first_match_gene <- if (p_val < 0.05) m_y else m_x
  return(first_match_gene)
}


# returns the degrees of freedom for two models [0,2]
degrees_of_freedom <- function(m_x, m_y) {
  if (((endsWith(m_y,"2")) & (endsWith(m_x, "0"))) | ((endsWith(m_y,"3")) & (endsWith(m_x, "2"))))
    dof <- 2
  else if(((endsWith(m_y,"3")) & (endsWith(m_x, "0"))))
    dof <- 4
  else if(((endsWith(m_y,"3")) & (endsWith(m_x, "1"))))
    dof <- 3
  else 
    dof <- 1
  return(dof)
}


# get the gene, the two models while the x is the more primitive one
# get also the names of the models m_x, m_y
chi_test_1_comp <- function(gene, model_x, model_y, m_x, m_y) {
  #  dof <- degrees_of_freedom(m_x,m_y)
  #  print(model_x)
  p_val <- pchisq(calculate_LR(gene, model_x, model_y),degrees_of_freedom(m_x, m_y), lower.tail = FALSE)
  return(p_val)
}


# chi square test for the second comparison - 
# M2 against M0 or M1
chi_test_2_comp <- function(gene, p_val,model_x, model_y, model_z) {
  
  first_match <- choose_better_model(p_val, "M0", "M1")
  p_val <- pchisq(send_right_model(gene,first_match, model_x, model_y, model_z),
                  degrees_of_freedom(first_match, "M2"), lower.tail = FALSE)
  return(p_val)
}

model_by_name <- function(second_comp, model0,model1, model2, model3) {
  
  tmp_model <- if (second_comp == "M0") model0 else if (second_comp == "M1") model1 else model2
  return(tmp_model)
  
}

# chi square test for the second comparison -
# M3 against M0 or M1 or M2
chi_test_3_comp <- function(gene, second_comp,model_x, model_y, model_z, model_q) {
  
  previous_comp <- model_by_name(second_comp, model_x, model_y, model_z, model_q)
  p_val <- pchisq(calculate_LR(gene, previous_comp, model_q),
                  degrees_of_freedom(second_comp, "M3"), lower.tail = FALSE)
  return(p_val)
}

# compare models --------------------------------------------------------------------
# 9.
# this section run the wilks test and the chi square test
# order results in dataframe

# first comparison - model M0 vs. M1
# for each gene, calculate the chi square distribution
m0_vs_m1 <- unlist(lapply(names(m1_data), function(x) { chi_test_1_comp(m1_data[[x]], model_M0[[x]], model_M1[[x]], "M0", "M1")}))

# compare M2 with the best model for each gene
m_num_vs_m2 <- unlist(lapply(1:length(names(m1_data)), function(x) chi_test_2_comp(m1_data[[x]], m0_vs_m1[[x]], model_M0[[x]], model_M1[[x]], model_M2[[x]])))

most_fitted_model <- data.frame(gene = names(m1_data), comp1 = m0_vs_m1, comp2 = m_num_vs_m2)

most_fitted_model <- most_fitted_model %>%
  mutate(first_comp_mod = derivedFactor("M1" = (comp1 < 0.05),"M0" = (comp1 >= 0.05),.method = "first",.default = NA)) %>%
  mutate(second_comp_mod = derivedFactor("M1" = (comp1 < 0.05),"M0" = (comp1 >= 0.05),.method = "first",.default = NA))

most_fitted_model <- most_fitted_model %>%
  mutate(second_comp_mod = case_when((comp2 < 0.05) ~ "M2",
                                     (comp2 >= 0.05) ~ as.character(first_comp_mod)))

# compare M3 with the best model for each gene 
m_num_vs_m3 <- unlist(lapply(1:length(names(m1_data)), function(x) chi_test_3_comp(m1_data[[x]], most_fitted_model$second_comp_mod[x], model_M0[[x]], model_M1[[x]], model_M2[[x]], model_M3[[x]])))

most_fitted_model$comp3 <- m_num_vs_m3 
most_fitted_model <- most_fitted_model %>% 
  mutate(third_comp_mod = second_comp_mod) %>%
  mutate(third_comp_mod = case_when((comp3 < 0.05) ~ "M3",
                                    (comp3 >= 0.05) ~ as.character(second_comp_mod)))



# parameters table ---------------------------------------------------------------------------
# 10.
# I build the table of the parameters of the models


parameter_4_each_gene <- function(genes, fitted_models, models_list, parameter_type) {
  
  param_4_genes <- sapply(1:length(genes), function(x) switch_parameter(genes[x], fitted_models[x], models_list, parameter_type))
  
}

switch_parameter <- function(gene, model, models_list, parameter_type) {
  
  return( switch(parameter_type,
                 x0 = {find_x(gene, model, models_list, "t0")},
                 x1 = {find_x(gene, model, models_list, "t1")},
                 B = {find_B(gene, model, models_list)},
                 t0 = {find_t0(gene, model, models_list)},
                 t1 = {find_t1(gene, model, models_list)},
                 gene = {return(gene)},
                 mod = {return(model)},
                 {return(NULL)}
  ))
}

# onset time
find_t0 <- function(gene, gene_fitted_model, models_list) {
  
  model_number <- as.numeric(stri_sub(gene_fitted_model, -1))
  t0 <- models_list[[model_number+1]][[gene]][[1]]$model$ind[1]
  return(t0)
}

# offset time
find_t1 <- function(gene, gene_fitted_model, models_list) {
  
  model_number <- as.numeric(stri_sub(gene_fitted_model, -1))
  t1 <- models_list[[model_number+1]][[gene]][[1]]$model$ind[length(models_list[[model_number+1]][[gene]][[1]]$model$ind)]
  return(t1)
}

# expression state (constant)
find_x <- function(gene, gene_fitted_model, models_list, t0_or_1) {
  if ((gene_fitted_model == "M0" | gene_fitted_model == "M1") & t0_or_1 == "t0") {
    model_number <- as.numeric(stri_sub(gene_fitted_model, -1))
    return(models_list[[model_number+1]][[gene]][[1]]$coefficients[1])
  }
  
  # in order to find x0 i should find t0
  t <- switch(t0_or_1,
              t0 = find_t0(gene, gene_fitted_model, models_list),
              t1 = find_t1(gene, gene_fitted_model, models_list))
  
  
  if (gene_fitted_model == "M2" | gene_fitted_model == "M3") {
    # take the right model
    model_number <- as.numeric(stri_sub(gene_fitted_model, -1))
    model_details <- models_list[[model_number+1]][[gene]][[1]]$coefficients
    m <- model_details[2]
    b <- model_details[1]
    
    return(m*t + b)
    
  } else {
    return(NA)
  }
}

# degradation rate
find_B <- function(gene, gene_fitted_model, models_list) {
  model_number <- as.numeric(stri_sub(gene_fitted_model, -1))
  return(models_list[[model_number+1]][[gene]][[1]]$coefficients[2])
}

parameter_list <- c("gene","mod", "x0", "x1", "B", "t0", "t1")
models_list <- list(model_M0, model_M1, model_M2, model_M3)

list_parameters_per_gene <- as.data.frame(sapply(parameter_list, function(p) parameter_4_each_gene(most_fitted_model$gene,
                                                                                                   most_fitted_model$third_comp_mod, 
                                                                                                   models_list, p)))

saveRDS(list_parameters_per_gene, "list_parameters_lm.RDS")

# plot model parameter distribution ----------------------------------------------------------
# 11.
# plot an histogram for each parameter of the models: X0, X1, t0, t1, B

plot_parameter_hist <- function(param_vec, param_type) {
  
  ggplot(param_vec, aes(x)) +
    geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth = 0.1) +
    labs(y = "% of genes", x = param_type, title = paste("Histogram of", param_type)) +
    theme(text = element_text(size = 18))
  
}

plot_params_list <- function(list_parameters_per_gene) {
  
  B <- data.frame(x = as.numeric(unlist(list_parameters_per_gene$B)))
  plot_B <- plot_parameter_hist(B, "B") 
  
  t0 <- data.frame(x = as.numeric(unlist(list_parameters_per_gene$t0)), mod = unlist(list_parameters_per_gene$mod))
  t0 <- t0 %>%
    filter(mod == "M2" | mod == "M3") %>%
    select(x)
  plot_t0 <- plot_parameter_hist(t0, "t0") 
  
  t1 <- data.frame(x = as.numeric(unlist(list_parameters_per_gene$t1)), mod = unlist(list_parameters_per_gene$mod))
  t1 <- t1 %>%
    filter(mod == "M3") %>%
    select(x)
  plot_t1 <- plot_parameter_hist(t1, "t1") 
  
  x0 <- data.frame(x = as.numeric(unlist(list_parameters_per_gene$x0)), mod = unlist(list_parameters_per_gene$mod))
  x0 <- x0 %>%
    filter(mod == "M3" | mod == "M2" | mod == "M1" | mod == "M0") %>%
    select(x)
  plot_x0 <- plot_parameter_hist(x0, "x0") 
  
  x1 <- data.frame(x = as.numeric(unlist(list_parameters_per_gene$x1)), mod = unlist(list_parameters_per_gene$mod))
  x1 <- x1 %>%
    filter(mod == "M3" | mod == "M2") %>%
    select(x)
  plot_x1 <- plot_parameter_hist(x1, "x1") 
  
  pdf(file = "parameters_graphs.pdf")
  par(mfrow=c(2,2))
  
  plot(plot_B)
  plot(plot_t0)
  plot(plot_t1)
  plot(plot_x0)
  plot(plot_x1)
  
  dev.off()
  
}

# --------------------------------------------------------------------------------------------


# take the first hour for genes starts to degrade after the fertilization
onset_deg <- unlist(sapply(m2_data, function(x) return(x$ind[1])))

# for models M0 fit "constant"
# for models M1 or M2 with onset before MZT fit "early model"
# for models M2 with onset after the MZT fit "late model"
most_fitted_model <- most_fitted_model %>%
  mutate(late_onset = onset_deg) %>%
  mutate(model = case_when((second_comp_mod == "M0") ~ "constant",
                           (second_comp_mod == "M1") ~ "early onset",
                           (second_comp_mod == "M2" & late_onset > mzt) ~ "late onset",
                           (second_comp_mod == "M2") ~ "early onset")) 


# save the models fitted for each gene
to_save <- most_fitted_model[,c("gene","model")]
saveRDS(to_save, file = "lm_fitted_calc.RDS")


# plot distributions --------------------------------------------------------------
# 10.
# plot graphs - 
# for counts of each linear group 
# for late onset distribution
# for half life distribution

# plot the percentage of count of each model 
ggplot(most_fitted_model, aes(model)) + geom_bar(position = "stack", aes(y = (..count..)/sum(..count..))) +
  labs(title = "Linear Models", subtitle = "according to chi square test", x = "model", y = "percentage") + 
  theme(text = element_text(size = 18))


# take only the late model genes 
late_mod_genes <- most_fitted_model %>%
  filter(model == "late onset")

png(filename = "late_onset_dist.png")

ggplot(late_mod_genes, aes(late_onset)) +
  geom_histogram(binwidth = 0.25) +
  labs(title = "Late Onset Time Distribution", x = "hpf", subtitle = paste0("starts from ", mzt, " hpf")) +
  theme(text = element_text(size = 18))

dev.off()

# take only the early model genes 
early_mod_genes <- most_fitted_model %>%
  filter(second_comp_mod == "M2") %>%
  select(late_onset) %>%
  mutate(early_dist = late_onset)

zero_onset <- most_fitted_model %>%
  filter(second_comp_mod == "M1") %>%
  mutate(early_dist = 0)

# ??????????????????????????????????????????????????????????????????????????????
# need to bind rows for early_mod_genes and zero_onset


png(filename = "early_onset_dist.png")

ggplot(early_mod_genes, aes(early_onset)) +
  geom_histogram(binwidth = 0.25) +
  labs(title = "Early Onset Time Distribution", x = "hpf", subtitle = paste0("ends at ", mzt, " hpf")) +
  theme(text = element_text(size = 18))

dev.off()


# half life calculation --------------------------------------------------------------
# 11.
# calculate half life for each gene and plot distribution

# take gene names for each model 
m0_names <- most_fitted_model %>%
  filter(model == "constant") %>%
  select(gene)

m1_names <- most_fitted_model %>%
  filter(model == "early onset") %>%
  select(gene)

m2_names <- most_fitted_model %>%
  filter(model == "late onset") %>%
  select(gene)

# choose for each gene its fitted model
m0_mod <- model_M0[m0_names$gene]
m1_mod <- model_M1[m1_names$gene]
m2_mod <- model_M2[m2_names$gene]

# get the intercept for model 0 for each gene
df_mod0 <- data.frame(intercept = unlist(lapply(m0_mod, function(x) x[[1]][[1]][1])), 
                      slope = 0) 
row.names(df_mod0) <- gsub("\\.\\(Intercept\\)","",row.names(df_mod0))

# get the intercept for model 1 for each gene
df_mod1 <- data.frame(intercept = unlist(lapply(m1_mod, function(x) x[[1]][[1]][1])),
                      slope = unlist(lapply(m1_mod, function(x) x[[1]][[1]][2])))
row.names(df_mod1) <- gsub("\\.\\(Intercept\\)","",row.names(df_mod1))

df_mod2 <- data.frame(intercept = unlist(lapply(m2_mod, function(x) x[[1]][[1]][1])),
                      slope = unlist(lapply(m2_mod, function(x) x[[1]][[1]][2])))
row.names(df_mod2) <- gsub("\\.\\(Intercept\\)","",row.names(df_mod2))

# bind these three data frames of the 3 models to one data frame
df_mods <- rbind(df_mod0, df_mod1, df_mod2)

# calculate half life out of slope
df_mods$half_life <- log(2)/(-df_mods$slope)
df_mods$half_life[df_mods$half_life > 20] <- 20

# plot half life distribution
png(filename = "half_life_dist.png")

ggplot(df_mods, aes(half_life)) +
  geom_histogram(bins = 200) +  
  xlim(0,max(df_mods$half_life[!is.na(df_mods$half_life)])) +
  labs(title = "half life distribution", subtitle = "model 1 and 2") +
  theme(text = element_text(size = 18))

dev.off()