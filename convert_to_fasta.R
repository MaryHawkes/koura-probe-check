#!/usr/bin/env Rscript

library(tidyverse)

# read the file
raw_probes <- readxl::read_xlsx("data/Order_3381_candidate_probes (1).xlsx",
                                col_names = TRUE) %>% 
    filter(!is.na(`target binding region`)) %>%
    mutate(fasta_name = paste(Gene, Probe, sep = "_probe"))

hmmer_ref <- "data/Trinity.fasta"

# don't ever do this. write a function instead
for (i in 1:nrow(raw_probes)){
    # make the fasta
    my_sq <- toupper(raw_probes[[i, "target binding region"]])
    my_name <- paste0(">", raw_probes[[i, "fasta_name"]])
    my_fasta <- list(my_name, my_sq)
    # write it to a file
    my_filename <- paste0(raw_probes[[i, "fasta_name"]], ".fasta")
    my_filepath <- paste0("output/probes/", my_filename)
    write_lines(my_fasta, my_filepath)    
    # run hmmer
    my_resultpath <- paste0("output/hmmer/", raw_probes[[i, "fasta_name"]], "_hits.txt")
    system2("nhmmer", args = c("--tblout", my_resultpath, my_filepath, hmmer_ref))
}

# read results back in
cn <- c("target name",
        "accession",
        "query name",
        "accession",
        "hmmfrom",
        "hmm to",
        "alifrom",
        "ali to",
        "envfrom",
        "env to",
        "sq len",
        "strand",
        "E-value",
        "score",
        "bias",
        "description of target1",
        "description of target2")

results_files <- list.files("output/hmmer",
                            pattern = "hits.txt",
                            full.names = TRUE,
                            recursive = FALSE)
names(results_files) <- sub("_hits.txt", 
                            "",
                            basename(results_files))
results_list <- lapply(results_files, 
                       read_table2, 
                       col_names = cn, 
                       comment = "#")
all_hmmer_results <- bind_rows(results_list, 
                               .id = "probe") %>%
    filter(!is.na(`target name`)) %>%
    mutate(should_hit = gsub("_probe.*$", "", probe))

# look for off targets
off_targets <- filter(all_hmmer_results, should_hit != `target name`)
write_csv(off_targets, "output/off_targets.csv")
