#setwd("/data/dragon226/awatson/dej_rip/bin/")

source('./05_visualise_results/fxns/functions.R')

source("master_eps_df.R")

plot_dbscan <- function(csv, dist, dunbrack_length, chain, cdrx, cdr_len){

  if (!file.exists(paste0(csv))) {
    print(paste0("XXX - ", csv, " does not exist."))
    return()
  }
  
  print(paste0("XXX - ", csv))
  df <- read_csv(csv, show_col_types=F)
  df <- select_k(df, chain)
  df <- join_dunbrack(df, dunbrack_length, chain, cdrx, cdr_len)
  print(unique(df$DBS))

  samps <- df %>% arrange(Coord1, Coord2) %>%
    group_by(DBS, Origin) %>%
    slice(sample(1:n(), 50, replace = TRUE)) %>%
    ungroup()
  
  # Export samples with colors to use in pymol
  write_csv(samps, file = paste0("../output/aligned/", "CDR", cdrx,
                                 "_Length", cdr_len, ".csv"))

  ov1 <- plot_overlay(df, cdrx, cdr_len)
  ov2 <- plot_overlay_bare(df, cdrx, cdr_len)
  ov3 <- plot_overlay_white(df, cdrx, cdr_len)
  dbs <- plot_dbs(df, dist, chain, cdrx, cdr_len)
  minidbs <- plot_minidbs(df, dist, chain, cdrx, cdr_len)
  logo <- plot_logo(df, cdrx, cdr_len)
  tabs <- plot_tables(df, dist, chain, cdrx, cdr_len)
  
}

pmap(master_df, plot_dbscan)
