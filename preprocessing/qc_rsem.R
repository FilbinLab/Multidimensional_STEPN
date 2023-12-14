args <- commandArgs(trailingOnly = TRUE)

rsem_dir <- args[1]

files <- list.files(rsem_dir, pattern = '*.cnt', full.names = T, recursive = T)

message('QC from samples')

files_name <- gsub('.cnt', '', basename(files))
sample <- unique(unlist(lapply(strsplit(files_name, '-'), `[[`, 1)))

qc <- lapply(files, function(x) {
  tmp <- read.table(x, header = F, nrows = 2, sep = '\t')
  c(sapply(strsplit(tmp[1,], ' '), `[`, 2), sapply(strsplit(tmp[1,], ' '), `[`, 4), sapply(strsplit(tmp[2,], ' '), `[`, 1), sapply(strsplit(tmp[2,], ' '), `[`, 2))
})

dt <- data.frame(name = gsub('-', '.', files_name), 
                 sample = sample, 
                 plate = sapply(strsplit(files_name, '-'), `[`, 2), 
                 well = sapply(strsplit(files_name, '-'), `[`, 3), 
                 nAligned = sapply(qc, `[`, 1), 
                 nTotal = sapply(qc, `[`, 2), 
                 nUniq = sapply(qc, `[`, 3), 
                 nMulti = sapply(qc, `[`, 4))

write.csv(dt, file.path(rsem_dir, paste0('qc_', sample, '.csv')), quote = F, row.names = F)
