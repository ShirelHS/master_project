# this script is:
# read the input files - the fpkm files contianing gene, pre-mRNA FPKM, mature-mRNA FPKM
# determine a threshold for filtering
# filer unexpressed genes in the results file
# save into new csv file

setwd(getwd())
rm(list = ls())

# load all libraries here.
# install if needed
requiredPackages = c('reshape','dplyr','gtools', 'tibble', 'stats',
                     'ggplot2', 'gridExtra', 'grid', 'data.table',
                     'mosaic', 'ggpubr', 'stringi', 'MLmetrics', 'readr')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
rm(requiredPackages,p)

# ------------------------------------------------------------------------------------
# 1.
# in this part read result files 
# assume results files in the same directory as getwd()
# create one expression table

# get the results file path (same as getwd())
path <- getwd()

#' Title
#' this function read files of FPKM of different samples or repetitions from the same directory
#' user should concentrate all fpkm files at the same directory w/o other files
#' which is called 'results"
#' the file names should be at the format of < seq_0h_SRR3228717.some_suffix >
#' the file format:
#'
#' @param path 
#'
#' @return list of all files content
#' @export
#'
#' @examples 
read_fpkm_files <- function(path) {
  
  # read all (isoforms) results files in a directory into a list of dataframes
  file_list <- list.files(path = paste0(path, "/results"), full.names=TRUE)
  file_list <- mixedsort(sort(file_list)) # sort files by names (in case the names starts with numbers), so 2 is after 1 (and not 10 is after 1)
  
  # read all the files in the directory into a list
  all_files <- lapply(file_list, read.table, sep = '\t',header = F, stringsAsFactors = T)
  
  return(all_files)
  
}


#' Title
#' this function ask for the hours from the user including repetitions.
#' change the files names in the list
#'
#'
#' @param files_content 
#'
#' @return al list of the: hours list from the user and the files FPKM content list.
#' @export
#'
#' @examples
fix_hours <- function(files_content) {
  
  # get the hours post fertilization of the different samples 
  hours <- readline(prompt = "Enter hours post fertilization in ascending numbers with commas saperating (including each repetition): ")
  hours <- as.numeric(unlist(strsplit(hours, ",")))
  
  # save hours for external memory
  saveRDS(hours, file = "hours.RDS")
  
  # for repetitions add the number of repetition using the function "repetition_num"
  occur_hours <- data.frame(table(hours))
  sample_name <- unlist(apply(occur_hours, 1,function(h) repetition_num(h[1], h[2]))) 
  
  
  if (length(sample_name) == (length(files_content))) {
    for (i in 1:length(sample_name)) {
      if (length(files_content[[i]]) == 3) {
        names(files_content[[i]]) <- c(paste0("gene_", sample_name[i]), 
                                       paste0("transcript_", sample_name[i]), 
                                       paste0("precursor_", sample_name[i]))
      }
      else {
        names(files_content[[i]]) <- c(paste0("gene_", sample_name[i]), 
                                       paste0("name_", sample_name[i]), 
                                       paste0("transcript_", sample_name[i]), 
                                       paste0("precursor_", sample_name[i])) 
      }
      files_content[[i]] <- files_content[[i]][order(files_content[[i]][,1]),]
    }
  } else {
    print("wrong number of hours post fertilization")
  }
  
  return(list(hours, files_content))
}


#' Title
#' this function returns the name of the repetition of the sample
#' for example: 1_1, 1_2, 3_2
#' the first number is the hour and the second the the repetition number
#'
#' @param hour 
#' @param num2_repeat 
#'
#' @return the name of the repetition of the sample
#' @export
#'
#' @examples
repetition_num <- function(hour, num2_repeat) {
  
  length_rep <- seq(1,num2_repeat)
  return(unlist(sapply(length_rep, function(x) paste(hour, x, sep = "_"))))
}

# ask about repetitions

#' Title
#' this function checks for repetitions in the datasets
#' it returns boolean flag of repetition or not
#'
#' @param hours 
#'
#' @return boolean flag of repetitions
#' @export
#'
#' @examples
is_include_repetitions <- function(hours) {
  
  occur_hours <- data.frame(table(hours))
  if (any(occur_hours$Freq > 1))
    return(TRUE)
  return(FALSE)
  
}


#' Title
#' this function makes a data frame of all expression levels
#' from the different samples of mature and precursor of each gene
#'
#' @param all_files 
#'
#' @return full_df_ordered - a data frame of all expression levels (FPKM)
#' @export
#'
#' @examples
order_files2_tbl <- function(all_files) {
  
  if (ncol(all_files[[1]]) == 4) {
    
    # read samples (hours) from the user
    df_exp <- data.frame(Reduce(cbind, all_files)) %>%
      select(-starts_with("gene"), -starts_with("name") ) %>%
      add_column(name = all_files[[1]][,2], .before = 1) %>%
      add_column(gene = all_files[[1]][,1], .before = 1) 
    
    df_ordered <- df_exp[,-c(1,2)]
    df_ordered <- df_ordered[, mixedsort(sort(colnames(df_ordered)))]
    full_df_ordered <- cbind(df_exp[,1:2],df_ordered)
    
  } else {
    
    # read samples (hours) from the user
    df_exp <- data.frame(Reduce(cbind, all_files)) %>%
      select(-starts_with("gene")) %>%
      add_column(gene = all_files[[1]][,1], .before = 1) 
    
    df_ordered <- df_exp[,-c(1)]
    df_ordered <- df_ordered[, mixedsort(sort(colnames(df_ordered)))]
    full_df_ordered <- cbind(gene = df_exp[,1],df_ordered)
    
  }
  
  
  return(full_df_ordered)
}


#' Title
#'
#' @param result_tbl
#'
#' @return
#' @export
#'
#' @examples
add_gene_name2_rsem <- function(result_tbl) {
  # read gff file 
  gff_path <- readline(prompt = "Enter path of gff file: ")
 
  if (endsWith(gff_path, 'gff3') == TRUE) {
    id_name_df <- gff_gene_names(gff_path)
  }
  else if (endsWith(gff_path, 'gtf') == TRUE) {
    id_name_df <- gtf_gene_names(gff_path)
  }
  
  # add colum of name to results file
  with_names <- result_tbl %>%
    inner_join(id_name_df, by = c("gene" = "id")) %>%
    relocate(name, .after = gene) 
  
  return(with_names)
}

gff_gene_names <- function(gff_path) {
  
  gff <- read.delim(file = gff_path, skip = 1001, header = F)
  genes_only <- gff %>% 
    filter(V3 == "gene") %>%
    select(V9)
  
  splited_info_genes <- strsplit(genes_only$V9, ";")
  id_name <- (unlist(lapply(splited_info_genes, function(x) x[1:2])))
  id <- gsub("ID=","",id_name[grep(pattern = "ID=", x = id_name)])
  name <- gsub("Name=","",id_name[grep(pattern = "Name=", x = id_name)])
  id_name_df <- data.frame(id,name)
  return(id_name_df)
}

gtf_gene_names <- function(gtf_path) {
  
  gff <- read.delim(file = gtf_path, header = F)[9]

  genes_only <- unique(lapply(gff$V9, function(x) return(str_split(x, ";")[[1]][c(1,3)])))
  
  details_striped <- lapply(genes_only, function(x) 
    return(list(gsub("gene_id ","", x[[1]][1]),gsub(" gene_name ","", x[[2]][1])))) 
  id_name_df <- do.call(rbind.data.frame, details_striped)
  colnames(id_name_df) <- c("id", "name")
  return(id_name_df)
}

#' Title
#'
#' @param result_tbl
#'
#' @return
#' @export
#'
#' @examples
add_gene_name2_rsem_gtf <- function(result_tbl) {
  # read gff file 
  gff_path <- readline(prompt = "Enter path of gff file: ")
  gff <- read.delim(file = gff_path, skip = 1001, header = F)
  genes_only <- gff %>% 
    filter(V3 == "gene") %>%
    select(V9)
  
  splited_info_genes <- strsplit(genes_only$V9, ";")
  id_name <- (unlist(lapply(splited_info_genes, function(x) x[1:2])))
  id <- gsub("ID=","",id_name[grep(pattern = "ID=", x = id_name)])
  name <- gsub("Name=","",id_name[grep(pattern = "Name=", x = id_name)])
  id_name_df <- data.frame(id,name)
  
  # add colum of name to results file
  with_names <- result_tbl %>%
    inner_join(id_name_df, by = c("gene" = "id")) %>%
    select(-name) %>%
    add_column(name = name, .before = 2) 
  
  return(with_names)
}

#' Title
#' ask from the user for gene names list
#' split them and return the vector
#'
#' @return vector of gene names
#' @export
#'
#' @examples slc35a5,dap,dcaf8,upf3b,ofd1,daxx,lcmt2,scfd1,nr4a1,dstyk
gene_names_from_user <- function() {
  
  # input example: slc35a5,dap,dcaf8,XB5959486 [provisional,upf3b,ofd1,daxx,lcmt2,scfd1,nr4a1,dstyk,unk,ncor2,bysl,dicer1,hsd17b8,fkbp8,ptpra,phgdh,man2b1,lmcd1,ppil2,pygm,mnat1,pygl,dcaf13,ankrd22
  genes_to_present <- readline(prompt = "Enter list of gene names saperated with commas only (example: buc,eomesa,slbp2): ")
  splited_express <- data.frame(gene_name = unlist(strsplit(genes_to_present, ",")))
  
  return(splited_express)
  
}


#' Title
#' write data frame into a csv file
#'
#' @param tble 
#' @param file_name 
#'
#' @return
#' @export
#'
#' @examples
write_csv <- function(tble, file_name) {
  
  # save the expression table in csv file
  write.csv(tble, file = paste0(file_name,".csv"))
  
}

# 0,0,2,2,4,4,6,6,6,6,8,8,8,8,12,12
#' Title
#' check all the functions above
#' like a main, but for check
#'
#' @return full data frame of expression and gene names 
#' @export
#'
#' @examples
check_main1 <- function() {
  #path <- getwd()
  
  file_list <- read_fpkm_files(path = path)
  hours_n_files <- fix_hours(file_list)
  hs <- hours_n_files[[1]]
  file_contetnt <- hours_n_files[[2]]
  is_repetition <- is_include_repetitions(hs)
  full_exp_df <- order_files2_tbl(file_contetnt)
  
  if (ncol(file_contetnt[[1]]) == 3) {
    full_exp_df <- add_gene_name2_rsem(full_exp_df)
  }
  
  colnames(full_exp_df)[1] <- "gene_id"
  colnames(full_exp_df)[2] <- "gene_name"
  if (length(which(full_exp_df$gene_id == "gene_id")) > 0) {
    full_exp_df <- full_exp_df[-which(full_exp_df$gene_id == "gene_id"),]
  }
  genes_to_plot <- gene_names_from_user()
  write_csv(full_exp_df,"exp_df")
  return(list(full_exp_df, genes_to_plot))
  
}   


# pre-mRNA and mRNA dist plot ------------------------------------------------------------------
# 2.
# this part is ploting an histogram for precursor and mature 
# and also scatter plot with density

extract_precursor_mature <- function(exp_of_all) {
  
  # generating the precursor data frame
  mature <- cbind(gene = exp_of_all$gene_id, name = exp_of_all$gene_name,
                  exp_of_all[,startsWith(colnames(exp_of_all), "transcript")])
  precursor <- cbind(gene = exp_of_all$gene_id, name = exp_of_all$gene_name,
                     exp_of_all[,startsWith(colnames(exp_of_all), "precursor")])
  
  return(list(mature, precursor))
  
}


#' Title
#' this function get the expression df and reorgenized according to mature and precursor
#'
#' @param exp_df 
#'
#' @return unified - a data frame of rows of mature and rows of precursor
#' @export
#'
#' @examples
unify_precursor_mature <- function(exp_df) {
  
  list_mat_pre <- extract_precursor_mature(exp_df)
  mature <- list_mat_pre[[1]]
  precursor <- list_mat_pre[[2]]
  
  colnames(mature) <- gsub("transcript_", "", colnames(mature))
  colnames(precursor) <- gsub("precursor_", "", colnames(precursor))
  
  mature$type <- rep("mature", nrow(mature))
  precursor$type <- rep("precursor", nrow(precursor))
  
  unified <- rbind(mature, precursor)
  
  return(unified)
}


#' Title
#'
#' @param exp_df 
#'
#' @return
#' @export
#'
#' @examples
take_max_exp_for_plot <- function(exp_df) {
  
  reorganized_df <- unify_precursor_mature(exp_df)
  reorganized_df$max <- as.numeric(unlist(apply(reorganized_df %>% select(-"gene",-"name",-"type"), 1, max)))
  reorganized_df$max <- reorganized_df$max + 0.00000001
  reorganized_df$log <- log2(reorganized_df$max)
  reorganized_df$log[reorganized_df$log < -10 | reorganized_df$log == -Inf] <- -10
  
  # order only the mature transcripts and take the threshold for 50% highly expressed genes
  mature_to_filter_by <- reorganized_df %>%
    filter(type == "mature") %>%
    arrange(desc(max), .by_group = FALSE) 
  threshold_to_filter <-  (mature_to_filter_by %>% slice_head(n = nrow(mature_to_filter_by)*0.5))[nrow(mature_to_filter_by)*0.5, c("log")]
  
  pdf(file = "exp_hist_pre_mat.pdf")
  par(mfrow=c(2,2))
  
  # here call the plot function
  plot(plot_exp_hist(reorganized_df, threshold_to_filter))
  
  dev.off()
}


plot_exp_hist <- function(df, threshold_num) {
  
  ggplot(data = df, aes(x = df$log, fill = type)) +
    geom_histogram(binwidth = 0.1, alpha = 0.4, position = "identity") +
    coord_cartesian(ylim = c(0,500)) +
    labs(title = "Histogram of Max Exp mRNA and pre-mRNA", x = "log2 FPKM", y = "count") +
    theme(text = element_text(size = 18)) +
    geom_vline(xintercept = threshold_num) +
    geom_text(aes(x = threshold_num, label = paste0("\n", round(threshold_num, 2)), y = 1), 
              colour = "black", angle = 90, text = element_text(size = 11))
  
}


pre_mat_max <- function(exp_df) {
  
  exp_df$log_max_transcript <- log2(as.numeric(unlist(apply(exp_df %>% select(starts_with("transcript_")), 1, max))) + 0.00000001)
  exp_df$log_max_precursor <- log2(as.numeric(unlist(apply(exp_df %>% select(starts_with("precursor_")), 1, max))) + 0.00000001)
  
  exp_df$log_max_transcript[exp_df$log_max_transcript < -10] <- -10
  exp_df$log_max_precursor[exp_df$log_max_precursor < -10] <- -10
  
  pdf(file = "exp_scatter_pre_mat.pdf")
  par(mfrow=c(2,2))
  
  plot(plot_exp_dens(exp_df))
  
  dev.off()
}

plot_exp_dens <- function(df) {
  
  ggplot(data = df, aes(x = log_max_transcript, y = log_max_precursor)) + 
    geom_hex(bins = 80) +
    scale_fill_continuous(limits = c(0, 150), type = "viridis") + 
    theme_bw() +
    geom_abline(slope = 1, intercept = 0, colour = 'red', size = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x = 0, label.y = 8) +
    labs(x = "mRNA (log2 FPKM)", y = "pre-mRNA (log2 FPKM)", title = "Max Expression Density Plot") +
    theme(text = element_text(size = 18))
  
}

check_main2 <- function() {
  
  listed_exp_genes <- check_main1()
  exp_to_filter <- listed_exp_genes[[1]]
  gene_names <- listed_exp_genes[[2]]
  
  take_max_exp_for_plot(exp_to_filter)
  pre_mat_max(exp_to_filter)
  
  return(list(exp_to_filter, gene_names))
  
  
}

# ---------------------------------------------------------------------------------------------
# 3.
# this part should filter the genes taking only the most expressed genes
# 50 percent of the top percent genes

#' Title
#' this function take the full data frame of expression levels 
#' filter the genes to the top 50% by the maximum expression of the mature transcript 
#' of each gene 
#'
#' @param full_df_ordered 
#'
#' @return filtered_df_ordered - filtered dataset for top 50% expressed genes
#' @export
#'
#' @examples
filter_data <- function(full_df_ordered) {
  
  # take only the mature transcripts 
  mature <- cbind(gene = full_df_ordered$gene_id, name = full_df_ordered$gene_name, full_df_ordered[,startsWith(colnames(full_df_ordered), "transcript")])
  
  # add a maximum column 
  mature$max <- unlist(apply(mature[,-c(1,2)], 1, max))
  
  # order by the maximum expression (for filtering)
  ordered_mature <- mature[order(mature$max, decreasing = T),]
  
  ordered_mature$log_on_max <- log2(as.numeric(ordered_mature$max))
  ordered_mature$log_on_max[ordered_mature$log_on_max < -10] <- -10
  
  # remove rows of log2(fpkm)=-10 expression
  only_expressed <- ordered_mature %>%
    filter(!(log_on_max == -10))
  
  # filter for top 50% expressed genes
  exp_top_50pc <- ordered_mature[1:(0.5*nrow(only_expressed)),] 
  
  # the filtered expression table 
  filtered_df_ordered <- exp_top_50pc %>%
    select(gene) %>%
    inner_join(full_df_ordered, by = c("gene" = "gene_id"))
  
  return(filtered_df_ordered)
}

# 0,0,2,2,4,4,6,6,6,6,8,8,8,8,12,12

#' Title
#' this function operates the check_main1 
#' and then checks the filter of genes  
#'
#' @return filtered data and gene names
#' @export
#'
#' @examples
check_main3 <- function() {
  
  listed_exp_genes <- check_main2()
  exp_to_filter <- listed_exp_genes[[1]]
  gene_names <- listed_exp_genes[[2]]
  
  filtered <- filter_data(exp_to_filter)
  
  return(list(filtered, gene_names))
  
}


# ------------------------------------------------------------------------------------------
# 4.
# normalize the expression 
# minimum threshold

#' Title
#' normalize the expression with ln
#'
#' @param filtered_dataset 
#'
#' @return log_expression - a data frame of expression after ln
#' @export
#'
#' @examples
normalize_exp <- function(filtered_dataset) {
  
  log_expression <- data.frame(filtered_dataset)
  as_numeric_exp <- as.data.frame(lapply(log_expression[,-c(1,2)], function(x) as.numeric(as.character(x))))
  log_expression[,-c(1,2)] <- as_numeric_exp + 0.1^7
  log_expression[,-c(1,2)] <- log2(log_expression[,-c(1,2)])
  
  return(log_expression)
}


#' Title
#' generate a minimal threshold for ln(fpkm)>-2.5
#'
#' @param log_exp 
#'
#' @return filtered_min_thre - a data frame with minimal threshold as -2.5   
#' @export
#'
#' @examples
minimal_exp <- function(log_exp) {
  
  only_exp <- log_exp[,-c(1,2)]
  only_exp[only_exp < -2.5] <- -2.5
  
  filtered_min_thre <- cbind(gene_id = log_exp$gene, gene_name = log_exp$gene_name, only_exp) 
  filtered_min_thre <- filtered_min_thre %>%
    mutate(gene_id = gsub("gene:","",gene_id))
  
  return(filtered_min_thre)
  
}

#' Title
#' this function operates the check_main2
#' then normalize the expression table  
#'
#' @return
#' @export
#'
#' @examples
check_main4 <- function() {
  
  list_genes_filter <- check_main3()
  filtered_data <- list_genes_filter[[1]]
  gene_names <- list_genes_filter[[2]]
  normalized <- normalize_exp(filtered_data)
  min_exp <- minimal_exp(normalized)
  write_csv(min_exp, "normalized_exp")
  
  return(list(min_exp, gene_names))
}


# --------------------------------------------------------------------------------------------
# 5.
# plot expression graphs for chosen genes

# map_id_name <- read.table("map_gene_id_name.txt", header = TRUE)
# genes_to_display <- filtered_min_thre %>%
#   left_join(map_id_name, by = c("gene_id" = "id")) %>%
#   inner_join(splited_express, by = c("name" = "gene_name"))


#' Title
#' this function gets a list of genes given by the user and 
#'
#' @param exp_tbl 
#' @param genes_to_plot 
#'
#' @return
#' @export
#'
#' @examples
get_specific_genes_exp <- function(exp_tbl, genes_to_plot) {
  
  genes_to_display <- exp_tbl %>%
    inner_join(genes_to_plot, by = c("gene_name" = "gene_name"))
  
  return(genes_to_display)
  
}


# function that plots the expression of precursor and mature transcript of a gene
# return graph object

#' Title
#' this function plots expression graph for gene given in a data frame with expression and precursor
#'
#' @param fpkm_df - a data frame of expression of gene (- mature transcript and precursor) 
#'
#' @return plot of expression of mature and precursor
#' @export
#'
#' @examples
plot_summed_iso <- function(fpkm_df) {
  
  current_gene_name <- fpkm_df$gene_name
  transcript_name <- "pre-mRNA"
  
  # put the pre-mRNA first in the list of isoforms so it always will get the same color in the plot
  pre_mRNA_row <- fpkm_df[,grep("precursor", colnames(fpkm_df))]
  isoform_row <- fpkm_df[,grep("transcript", colnames(fpkm_df))]
  
  # change the names to hours and repetitions only
  colnames(pre_mRNA_row) <- gsub("precursor_","",colnames(pre_mRNA_row))
  colnames(isoform_row) <- gsub("transcript_","",colnames(isoform_row))
  
  data_pre_first <- rbind.data.frame(pre_mRNA_row, isoform_row)
  
  # if there are repetitions to the same gene:
  hours <- as.numeric(gsub(pattern = "_\\d+", replacement = "", x = colnames(data_pre_first)))
  
  stacked <- stack(data_pre_first[,1:(length(colnames(data_pre_first)))])
  stacked$ind <- rep(hours,each = 2)
  stacked$part <- rep(c("precursor","transcript"), length(hours))
  
  
  graph <- ggplot(data = stacked, aes(x = ind, y = values, color = part, group = part)) + 
    geom_point() + 
    stat_summary(aes(group = part), geom = "line") +
    labs(title = current_gene_name, x = "hours post fertilization", y = "log2(FPKM)") +
    coord_cartesian(ylim = c(-5, 15), xlim = c(0, max(stacked$ind))) +
    theme(text = element_text(size = 18), axis.text = element_text(size = 18), 
          legend.text = element_text(size = 18))
  
}


#' Title
#' plot genes into a one long pdf file
#'
#' @param genes_to_display 
#'
#' @return no return
#' @export
#'
#' @examples
plot_2_pdf <- function(genes_to_display) {
  
  # call function to plot summed isoforms:
  # call function plot_fpkm and create a pdf with all maternal transcripts' expression 
  pdf(file = "expression_graphs.pdf")
  par(mfrow=c(2,2))
  
  for (i in 1:nrow(genes_to_display)) {
    par(mar = c(5,4,4,5)+.1)
    my_plot <- plot_summed_iso(genes_to_display[i,])
    print(my_plot)
    
  }
  
  dev.off()
  
}

check_main5 <- function() {
  
  # 0,0,2,2,4,4,6,6,6,6,8,8,8,8,12,12
  # buc,eomesa,dazl,hmgb3b,btg4,tbxta,celf1,slbp2
  # first get the data from the previuos check_main (3,2,1...)
  list_exp_genes <- check_main4()
  expression_tbl <- list_exp_genes[[1]]
  gene_names <- list_exp_genes[[2]]
  g_names_exp <- get_specific_genes_exp(expression_tbl, gene_names)
  plot_2_pdf(g_names_exp)
  
  return(expression_tbl)
  
}


# ---------------------------------------------------------------------------------------------
# 6.
# plot histograms for mature and precursors


#' Title
#' this function plots histogram for a vector of expression
#'
#' @param one_hour_exp - exp vector
#' @param current_hour - name of hour to present
#' @param sub_ - sub title to present in the plot
#' @param threshold - the threshold for minimum expression (after ln, the values might be very low)
#'
#' @return a plot 
#' @export
#'
#' @examples
hist_ggplot <- function(one_hour_exp, current_hour, sub_, threshold) {
  
  df_hour_exp <- data.frame(one_hour_exp)
  title_ <- paste("Distribution of FPKM at hour", current_hour)
  above_threshold <- sum(one_hour_exp > threshold)
  below_threshold <- sum(one_hour_exp < threshold)
  sub_title <- paste(paste(above_threshold, sub_,"above",threshold) , paste(below_threshold, sub_,"below", threshold), sep = "\n")
  
  ggplot(data = df_hour_exp, aes(x = one_hour_exp)) + geom_histogram(bins = 50) +
    labs(title = title_, subtitle = sub_title, x = "fpkm levels (log2 scale)") + 
    coord_cartesian(ylim = c(0,3000)) +
    geom_vline(xintercept = threshold, colour = 'red') +
    theme(plot.subtitle = element_text(size = 8))
}


#' Title
#' this function generates histogram plots for the different samples (hours) 
#'
#' @param exp_df - full expression data frame 
#' @param hours - sample hour
#' @param trans_type - mature or precursor
#' @param pdf_name - name of file to save in it
#'
#' @return no return
#' @export
#'
#' @examples
generate_hist_plot <- function(exp_df, hours, trans_type, pdf_name) {
  
  count = 1
  my_plots <- list()
  for (i in 1:length(hours)) {
    my_plots[[count]] <- hist_ggplot(exp_df[,i+2],hours[i],trans_type, threshold = -2.5)
    count <- count + 1
  }
  
  ml <- marrangeGrob(my_plots, nrow=2, ncol=2)
  ## non-interactive use, multipage pdf
  ggsave(pdf_name, ml, width = 6.5, height = 6.5)
  
}


check_main6 <- function() {
  
  exp_df <- check_main5()
  hs <- readRDS("hours.RDS")
  list_mat_pre <- extract_precursor_mature(exp_df)
  mature <- list_mat_pre[[1]]
  precursor <- list_mat_pre[[2]]
  generate_hist_plot(mature, hs, "mature transcript", "mRNA_hist.pdf")
  generate_hist_plot(precursor, hs, "pre-mRNA", "pre_mRNA_hist.pdf")
  
  return(list(exp_df, hs))
}


# ---------------------------------------------------------------------------------------------
# 7.
# plot scatter plot for expression levels

# function that creates scatter plots with density and abline
scatter_graph <- function(cols, h1, h2, ab_line, max_count) {
  
  # take only the two colonms into a df in order to remove the genes that doesn't express at all in both hours. 
  #  my_data <- my_data[!(my_data[,1] == ln(0.00001) & my_data[,2] == ln(0.00001)),c(1,2)] # ?????? but why??????
  title_ <- paste("log2(FPKM) at hours ", h1, " and ",h2)
  
  return(ggplot(data = cols, aes(x = cols[,1], y = cols[,2])) + geom_hex(bins = 80) +
           scale_fill_continuous(limits = c(0, max_count), type = "viridis") + theme_bw() + 
           coord_cartesian(ylim = c(-10, 18), xlim = c(-10, 18)) + 
           geom_abline(slope = 1, intercept = 0, colour = 'red', size = 0.3) +
           geom_vline(xintercept = ab_line, colour = 'mediumorchid2') +
           geom_hline(yintercept = ab_line, colour = 'mediumorchid2') +
           labs(title = title_, x = paste(h1, " hpf"), y = paste(h2, " hpf")) +
           theme(text = element_text(size = 15)))
}


generate_scatter_plots <- function(mature, hours, pdf_name) {
  # not sure if needed:
  # take the hours:
  # create list of graphs of each two following hours
  graph_list <- list()
  count = 1
  ab_lin <- -2.5
  max_coun <- 350
  
  for (col in 3:(length(hours)+1)) {
    graph_list[[count]] <- scatter_graph(mature[,c(col,col+1)], hours[col-2], hours[col - 1], ab_lin, max_coun) 
    count <- count + 1
    
  }
  
  # save the graphs into a file
  ml <- marrangeGrob(graph_list, nrow = 1, ncol = 1)
  ## non-interactive use, scatter_plot.pdf
  ggsave(pdf_name, ml, width = 5.5, height = 5.5)
  
}


check_main7 <- function() {
  
  # zhao - 0,0,2,2,4,4,6,6,6,6,8,8,8,8,12,12
  # meyer - 2.5,2.5,2.5,3.3,3.3,3.3,4.5,4.5,4.5,5.3,5.3,5.3,10,10,10
  # buc,eomesa,dazl,hmgb3b,btg4,tbxta,celf1,slbp2,eve1,sp5l,apoc2,actb2,ybx1,apoeb,dynll1
  exp_df <- check_main1()[[1]]
  exp_df_normal <- normalize_exp(exp_df)
  exp_only <- exp_df_normal[, -c(1,2)]
  exp_only[exp_only < -10] <- -10
  with_gene_names <- cbind(exp_df_normal[,c(1,2)], exp_only)
  
  mature <- extract_precursor_mature(with_gene_names)[[1]]
  
  hours <- gsub("precursor_|transcript_", "",names(mature)[-c(1,2)])
  
  generate_scatter_plots(mature, hours, "scatter_plots.pdf")
  
}




