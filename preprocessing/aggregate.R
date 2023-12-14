args <- commandArgs(trailingOnly = TRUE)

rsem_dir <- args[1]

files <- list.files(rsem_dir, pattern = '*.genes.results', full.names = T)

message('Aggregating samples')

files_name <- gsub('.genes.results', '', basename(files))
sample <- unique(unlist(lapply(strsplit(files_name, '-'), `[[`, 1)))

## Reading files
tab <- lapply(files, function(x) {
  read.table(x, sep = '\t', header = T, row.names = 1)[, c('TPM', 'expected_count'), drop = F]
})

## TPM
tpm <- do.call(cbind, lapply(tab, "[", , 'TPM', drop = F))
colnames(tpm) <- files_name
write.csv(tpm, file.path(rsem_dir, paste0('cm_tpm_', sample, '.csv')), quote = F)

## COUNTS
counts <- do.call(cbind, lapply(tab, "[", , 'expected_count', drop = F))
colnames(counts) <- files_name
write.csv(round(counts), file.path(rsem_dir, paste0('cm_counts_', sample, '.csv')), quote = F)
