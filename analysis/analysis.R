
load_column <- function(sample_name, column_name) {
  results_path <- file.path(base_dir, 'results', sample_name,
    c('paper', 'paper_forward', 'paper_reverse'),
    'express', 'results.xprs')
  res <- lapply(results_path, read.table, header = TRUE, stringsAsFactors = FALSE)
  res <- lapply(res, function(x) x[, c('target_id', column_name)])
  res <- Map(function(df, new_name) {
    data.table::setnames(df, column_name, new_name)
    df
  }, res, paste0(c('default', 'forward', 'reverse'), '_', column_name))

  Reduce(function(x, y) dplyr::inner_join(x, y, by = 'target_id'), res)
}

base_dir <- '..'
current_sample <- 'ENCLB037ZZZ'
current_sample_fpkm <- load_column(current_sample, 'fpkm')

ggplot(current_sample_fpkm, aes(log(forward_fpkm + 1), log(reverse_fpkm + 1))) +
  geom_point(alpha = 0.4)

load_all_samples <- function(which_mode, column_name) {
  which_path <- switch(which_mode,
    default = 'paper',
    forward = 'paper_forward',
    reverse = 'paper_reverse')

  sample_names <- c('ENCLB037ZZZ', 'ENCLB038ZZZ', 'ENCLB055ZZZ', 'ENCLB056ZZZ')
  results_path <- file.path(base_dir, 'results', sample_names, which_path,
    'express', 'results.xprs')

  res <- lapply(results_path, read.table, header = TRUE, stringsAsFactors = FALSE)
  res <- lapply(res, function(x) x[, c('target_id', column_name)])
  res <- Map(function(df, new_name) {
    data.table::setnames(df, column_name, new_name)
    df
  }, res, sample_names)

  res <- Reduce(function(x, y) dplyr::inner_join(x, y, by = 'target_id'), res)
  res <- as.data.frame(res)
  rownames(res) <- res$target_id
  res$target_id <- NULL

  res
}

fasta_ordering <- system("grep '>' ../index/gencode.v16.pc_transcripts.fa | sed 's/^>//g'",
  intern = TRUE)

reverse_tpm <- load_all_samples('reverse', 'tpm')
reverse_tpm <- reverse_tpm[fasta_ordering, ]

# this is what was submitted to http://rafalab.rc.fas.harvard.edu/rnaseqbenchmark
write.table(reverse_tpm, file = '../results/reverse_merged.tpm',
  quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)


forward_tpm <- load_all_samples('forward', 'tpm')

forward_plot <- ggplot(forward_tpm, aes(log(ENCLB037ZZZ + 1), log(ENCLB038ZZZ + 1))) +
  geom_point(alpha = 0.4) +
  ggtitle('running express with forward strand')

reverse_plot <- ggplot(reverse_tpm, aes(log(ENCLB037ZZZ + 1), log(ENCLB038ZZZ + 1))) +
  geom_point(alpha = 0.4) +
  ggtitle('running express with reverse strand')

plot_grid(forward_plot, reverse_plot)

reverse_counts <- load_all_samples('reverse', 'est_counts')
forward_counts <- load_all_samples('forward', 'est_counts')

# NOTE: very few reads map if you do not specify the correct strand (1.4-1.8M versus 60-75M reads)

sapply(forward_counts, sum)
sapply(reverse_counts, sum)


###
# figuring out what unit they are using
###

tmp <- read.table('http://rafalab.rc.fas.harvard.edu/rnaseqcomp/encodeexample.txt',
  header = TRUE)

# appears to be tpm
sapply(tmp, sum)
