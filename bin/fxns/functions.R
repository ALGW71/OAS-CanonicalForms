library(tidyverse)
library(gridExtra)
library(ggseqlogo)
library(grid)

# load CDR seqs
cdr_seqs <- read_rds(file = "../resources/rds/cdr_seqs.rds")

group.colors <- c("0" = "black", "1" = "blue", "2" = "red", "3" ="green",
                  "4" = "yellow", "5" = "magenta", "6" = "purple")


select_k <- function(df, chain){

  ct <- df %>% pivot_longer(4:ncol(.), names_to="K", values_to="DBS") %>%
                mutate(kNN = as.numeric(str_extract(pattern = ".$", string = K))) %>%
                mutate(DBS == as.numeric(DBS)) %>%
                group_by(K, kNN) %>% summarise(DBS = max(DBS), .groups="drop")

  if (nrow(ct %>% filter(kNN == DBS)) > 0){
    ct <- ct %>% filter(kNN == DBS) %>% arrange(desc(K)) %>% slice(1) %>% pull(K)

  }else if (nrow(ct %>% filter((kNN+1) == DBS)) > 0) {
    ct <- ct %>% filter((kNN+1) == DBS) %>% arrange(desc(K)) %>% slice(1) %>% pull(K)

  }else if (nrow(ct %>% filter((kNN-1) == DBS)) > 0) {
    ct <- ct %>% filter((kNN-1) == DBS) %>% arrange(desc(K)) %>% slice(1) %>% pull(K)

  }else{
    ct <- ct %>% arrange(desc(DBS)) %>% slice(1) %>% pull(K)

  }

  print(ct)
  df <- df %>% dplyr::select(all_of(c("Coord1", "Coord2","PDB", !!ct))) %>% rename("DBS"=!!ct)
  
  df <- df %>% 
    arrange(desc(PDB)) %>%
    mutate(path = PDB, .before=2) %>%
    separate(PDB, into = c("PDB", "chain"), sep = "_chain_", remove = F, fill = "right")  %>%
    mutate(chain = str_replace(chain, pattern = "_Ig(.*)", replacement = "")) %>%
    unite(PDB, chain, sep = "_", col = "join", remove = F)  %>%
    select("Coord1", "Coord2","PDB", "DBS", "chain", "join", "path")

  if(chain == "H"){
    fails <- read_csv("../resources/csvs/anarci_fails_heavy.csv", show_col_types = F)
  }else if(chain == "L"){
    fails <- read_csv("../resources/csvs/anarci_fails_light.csv", show_col_types = F)
  }
  
  df <- df %>% filter(!PDB %in% fails$ID)
  return(df)
}

join_dunbrack <- function(df, dunbrack_length, chain, cdrx, cdr_len){
  
  cdr_id <- paste0("cdr", str_extract(pattern = ".$", string = cdrx), "_len")
  sabdab_dunbrack <- read_csv(paste0("../resources/csvs/",
                                     cdrx,
                                     "_cdrs.csv"),
                              show_col_types = F) %>%
    # ==========================================================================
    # This is an issue...
    filter(get(cdr_id) == !!cdr_len) %>%
    # ==========================================================================
    mutate(PDB = tolower(pdb), .before=1) %>%
    select(-pdb) %>%
    arrange(desc(PDB)) %>%
    unite(PDB, dunbrack_chain_id, sep = "_", col = "join", remove = F)

  vis_mds <- df %>%
    mutate(Origin = case_when(join %in% sabdab_dunbrack$join ~ "SAbDab-Dubrack",
                              T ~ "ImmuneBuilder")) %>%
    left_join(sabdab_dunbrack %>% select(-PDB, -chain), by = "join") %>%
    filter(!(is.na(cluster_nocutoff)&Origin=="SAbDab-Dubrack")) %>%
    mutate(join = join %>% str_remove(pattern = "_NA")) %>%
    filter(((Origin == "SAbDab-Dubrack") & (resolution <= 3.5) & 
              (cdr_len == dunbrack_length)) | 
             Origin == "ImmuneBuilder")

    #print(paste0("Row difference = ", nrow(df) - nrow(vis_mds))  )
  return(vis_mds)

}

get_dunbrack <- function(cdrx, cdr_len){
  cdr_id <- paste0("cdr", str_extract(pattern = ".$", string = cdrx), "_len")
  sabdab_dunbrack <- read_csv(paste0("../resources/csvs/",
                                     cdrx,
                                     "_cdrs.csv"),
                              show_col_types = F) %>%
    filter(get(cdr_id) == !!cdr_len)

    return(sabdab_dunbrack)
}

plot_overlay <- function(dbs, cdrx, cdr_len){
  
  pt <- dbs %>%
    ggplot(aes(x = Coord1, y = Coord2)) +
    geom_point(data = dbs %>% filter(Origin == "ImmuneBuilder"),
               aes(x = Coord1, y = Coord2),
               color = "black",
               alpha = 0.25,
               size = 0.5) +
    geom_point(data = dbs %>% filter(Origin == "SAbDab-Dubrack"),
               aes(x = Coord1, y = Coord2,
                   color = cluster),
               alpha = 0.75,
               size = 1.5) +
    theme_classic() +
    theme(legend.position = "right",
          aspect.ratio = 5/5,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          legend.text = element_text(size = 20,  color = "black"),
          text = element_text(size = 20,  color = "black")) +
    guides(color=guide_legend(ncol=1, byrow=TRUE)) +
    ggtitle(paste0("Predictions: ", 
                   nrow(dbs %>% filter(Origin == "ImmuneBuilder") %>%
                          select(join) %>% distinct), "\n",
                   "Plotted PDBs: ", 
                   nrow(dbs %>% filter(!Origin == "ImmuneBuilder") %>% 
                          select(join) %>% distinct),  "\n",
                   "PyIg in train data: ",
                   nrow(get_dunbrack(cdrx, cdr_len))))

    # Define the directory and file paths
    output_dir <- paste0("../output/cdrs/", cdrx, "_", cdr_len)
    # Check if the directory exists, if not create it
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  
    ggsave(plot = pt, paste0("../output/cdrs/", 
    cdrx, "_", cdr_len, "/",
    cdrx, "_", cdr_len, "_", "overlay.svg"), 
         device = "svg", 
         height = 7, width = 6)
  return() 
}

plot_overlay_bare <- function(dbs, cdrx, cdr_len){
   
  pt <- dbs %>%
    ggplot(aes(x = Coord1, y = Coord2)) +
    geom_point(data = dbs %>% filter(Origin == "ImmuneBuilder"),
               aes(x = Coord1, y = Coord2),
               color = "black",
               alpha = 0.25,
               size = 0.5) +
    geom_point(data = dbs %>% filter(Origin == "SAbDab-Dubrack"),
               aes(x = Coord1, y = Coord2,
                   color = cluster),
               alpha = 0.75,
               size = 1.5) +
    theme_classic() +
    theme(legend.position = "right",
          aspect.ratio = 5/5,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          legend.text = element_text(size = 20,  color = "black"),
          text = element_text(size = 20,  color = "black")) +
    guides(color=guide_legend(ncol=1, byrow=TRUE))

    # Define the directory and file paths
    output_dir <- paste0("../output/cdrs/", cdrx, "_", cdr_len)
    # Check if the directory exists, if not create it
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
 
    ggsave(plot = pt, paste0("../output/cdrs/", 
    cdrx, "_", cdr_len, "/",
    cdrx, "_", cdr_len, "_", "overlay_bare.svg"), 
         device = "svg", 
         height = 5, width = 6)
  return() 
  
}

plot_overlay_white <- function(dbs, cdrx, cdr_len){
   
  pt <- dbs %>%
    ggplot(aes(x = Coord1, y = Coord2)) +
    geom_point(data = dbs %>% filter(Origin == "ImmuneBuilder"),
               aes(x = Coord1, y = Coord2),
               color = "white",
               alpha = 0.25,
               size = 0.5) +
    geom_point(data = dbs %>% filter(Origin == "SAbDab-Dubrack"),
               aes(x = Coord1, y = Coord2,
                   color = cluster),
               alpha = 0.75,
               size = 1.5) +
    theme_classic() +
    theme(legend.position = "right",
          aspect.ratio = 5/5,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          legend.text = element_text(size = 20,  color = "black"),
          text = element_text(size = 20,  color = "black")) +
    guides(color=guide_legend(ncol=1, byrow=TRUE))

    # Define the directory and file paths
    output_dir <- paste0("../output/cdrs/", cdrx, "_", cdr_len)
    # Check if the directory exists, if not create it
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
 
    ggsave(plot = pt, paste0("../output/cdrs/", 
    cdrx, "_", cdr_len, "/",
    cdrx, "_", cdr_len, "_", "overlay_white.svg"), 
         device = "svg", 
         height = 5, width = 6)
  return() 
  
}

plot_dbs <-  function(df, dist, chain, cdrx, cdr_len){
  
  dbs <- df
  # Calculate the centroid of each cluster
  centroids <- dbs %>%
    filter(Origin == "ImmuneBuilder") %>%
    dplyr::filter(DBS != "0") %>%  # Assuming "0" is noise
    dplyr::group_by(DBS) %>%
    dplyr::summarise(
      Coord1 = mean(Coord1, na.rm = TRUE),
      Coord2 = mean(Coord2, na.rm = TRUE), 
      .groups="drop"
    )
  
  print(centroids)

  # Compute the point closest to the centroid for each cluster
  closest_points <- dbs %>%
    filter(Origin == "ImmuneBuilder") %>%
    dplyr::filter(DBS != "0") %>%  # Assuming "0" is noise
    dplyr::group_by(DBS) %>%
    dplyr::left_join(centroids, by = "DBS", suffix = c("", "_centroid")) %>%
    dplyr::mutate(distance = sqrt((Coord1 - Coord1_centroid)^2 + (Coord2 - Coord2_centroid)^2)) %>%
    dplyr::filter(distance == min(distance)) %>%  # Select the point with minimum distance to the centroid
    dplyr::select(join, DBS) %>%
    dplyr::rename(Centroid = join)
  
  dbs <- dbs %>% left_join(closest_points, by = "DBS") %>% mutate(DBS=as.character(DBS))
  
  pt <- dbs %>% ggplot(aes(x = Coord1, y = Coord2)) +
    geom_point(data = dbs %>% filter(Origin == "ImmuneBuilder"),
               aes(x = Coord1, y = Coord2, color = DBS),
               alpha = 0.25,
               size = 0.5) +
    geom_point(data = dbs %>% filter(join == Centroid),
               aes(x = Coord1, y = Coord2),
               color = "black",
               shape = 4,
               stroke = 1.5,
               size = 1.5) +
    geom_text(data = dbs %>% filter(join == Centroid),
                aes(x = Coord1, y = Coord2, label=DBS),
                color = "black",
                fontface="bold",
                hjust = -0.5, 
                vjust = -0.5,
                size = 6) +
    theme_classic() +
    theme(legend.position = "none",
          aspect.ratio = 5/5,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(size = 20,  color = "black"))  +
    scale_color_manual(name = element_blank(), values = group.colors)

    ggsave(plot = pt, paste0("../output/cdrs/", 
    cdrx, "_", cdr_len, "/",
    cdrx, "_", cdr_len, "_", "dbs.svg"), 
         device = "svg", 
         height = 5, width = 4)
  return()
}

plot_minidbs <-  function(df, dist, chain, cdrx, cdr_len){
  
  df <- df %>% filter(Origin == "ImmuneBuilder") %>%
    mutate(DBS=as.character(DBS))

  pt <- df %>%
    ggplot(aes(x = Coord1, y = Coord2, color = DBS)) +
    geom_point(data = df %>% filter(DBS == "0"),
               aes(x = Coord1, y = Coord2, color = DBS),
               alpha = 0.5, size = 0.5) +
    geom_point(data = df %>% filter(!DBS == "0"),
               aes(x = Coord1, y = Coord2, color = DBS),
               alpha = 1.0, size = 0.5) +
    theme_classic() +
    theme(legend.position = "none",
          aspect.ratio = 5/5,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          text = element_blank())  +
    scale_color_manual(name = element_blank(), values = group.colors)

    ggsave(plot = pt, paste0("../output/cdrs/", 
    cdrx, "_", cdr_len, "/",
    cdrx, "_", cdr_len, "_", "minidbs.svg"), 
         device = "svg", 
         height = 3, width = 3)
  return()
}

plot_tables <- function(df, dist, chain, cdrx, cdr_len){
  
  dbs <- df
  # Calculate the centroid of each cluster
  centroids <- dbs %>%
    dplyr::filter(DBS != "0") %>%  # Assuming "0" is noise
    dplyr::group_by(DBS) %>%
    dplyr::summarise(
      Coord1 = mean(Coord1, na.rm = TRUE),
      Coord2 = mean(Coord2, na.rm = TRUE), 
      .groups="drop"
    )
  
  # Compute the point closest to the centroid for each cluster
  closest_points <- dbs %>%
    filter(Origin == "ImmuneBuilder") %>%
    dplyr::filter(DBS != "0") %>%  # Assuming "0" is noise
    dplyr::group_by(DBS) %>%
    dplyr::left_join(centroids, by = "DBS", suffix = c("", "_centroid")) %>%
    dplyr::mutate(distance = sqrt((Coord1 - Coord1_centroid)^2 + (Coord2 - Coord2_centroid)^2)) %>%
    dplyr::filter(distance == min(distance)) %>%  # Select the point with minimum distance to the centroid
    dplyr::select(join, DBS) %>%
    dplyr::rename(Centroid = join)
  
  dbs <- dbs %>% left_join(closest_points, by = "DBS") %>% mutate(DBS=as.character(DBS))

  plots <- get_numbers(dbs, dist)

  pt <- grid.arrange(grobs = plots, layout_matrix = rbind(c(1,1,2,2), c(1,1,2,2)), top = title)

  ggsave(plot = pt, paste0("../output/cdrs/", 
    cdrx, "_", cdr_len, "/",
    cdrx, "_", cdr_len, "_", "tables.svg"), 
         device = "svg", 
         height = 7, width = 6)

}

get_numbers <- function(df, dist){

   centroids <- df %>%
    filter(join == Centroid) %>%
    select(join, DBS) %>%
    arrange(DBS)

  dt <- as.matrix(read_csv(dist, show_col_types=F))

  meanVal <- mean(dt[upper.tri(dt)])
  medVal <- median(dt[upper.tri(dt)])
  sdVal <- sd(dt[upper.tri(dt)])
  
  if (nrow(centroids) > 1) {
    rownames(dt) <- colnames(dt)
    dt <- dt[, colnames(dt) %in% centroids$join]
    dt <- dt[rownames(dt) %in% centroids$join, ]
    
    read_values <- dt %>% as_tibble(rownames = "join") %>% 
      left_join(centroids, by = "join") %>%
      pivot_longer(2:(nrow(centroids)+1), names_to = "target", values_to = "RMSD") %>% 
      rename("Centroid1" = DBS) %>%
      left_join(centroids, by = c("target"="join")) %>%
      select(-join, -target) %>% 
      rename("Centroid2" = DBS) %>%
      select(Centroid1, Centroid2, RMSD) %>%
      arrange(Centroid1, Centroid2) %>% 
      filter(Centroid1 < Centroid2)
    
  }else {
    read_values <- tibble(Centroid1 = c(), Centroid2 = c(), RMSD = c())
  }
  
  summary_vals <- tibble(Centroid1 = c("Mean", "Median", "SD"), Centroid2 = c("", "", ""), RMSD = c(meanVal, medVal, sdVal))
  
  # Bind the new tibble to the bottom of read_values
  read_values <- bind_rows(read_values, summary_vals) %>% mutate(RMSD = round(RMSD, 2))
  # Convert the read_values table to a plot
  table_plot <- ggpubr::ggtexttable(df %>% group_by(DBS, Origin) %>% summarise(N = n(), .groups="drop") %>% arrange(Origin), rows = NULL)
  plots <- list(table_plot)

  table_plot2 <- ggpubr::ggtexttable(read_values, rows = NULL)
  plots <- c(plots, list(table_plot2))

  return(plots)

}

plot_logo <- function(df, cdrx, cdr_len){
  
  dbs <- df
  print(nrow(df))
  seqs_col_name <- paste0("cdr", str_extract(cdrx, pattern = ".$"), "_aa_",
                          if_else(str_extract(cdrx, pattern = "^.") == "H", "heavy", "light"))
  
  seqs <- cdr_seqs %>%
    select(ID, all_of(seqs_col_name)) %>%
    filter(ID %in% dbs$PDB)  %>%
    mutate(AA = !!sym(seqs_col_name)) %>%
    select(-all_of(seqs_col_name))
  
  dbs <- dbs %>% left_join(seqs, by = c("PDB" = "ID"))
  print(nrow(dbs))

  # logo plots
  lp_df <- dbs %>% select(DBS, AA) %>%
    filter(!is.na(AA), !DBS==0) %>%
    arrange(DBS) %>%
    mutate(DBS = factor(DBS, levels = unique(DBS)))
  
  lp <- unstack(lp_df, AA ~ DBS)
  
  lpt <- ggseqlogo(lp, ncol = 1) + theme_light() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          legend.text = element_text(size = 20,  color = "black"),
          text = element_text(size = 20,  color = "black")) +
    theme(strip.background =element_rect(colour = 'black')) +
    theme(strip.text = element_text(colour = 'black', size =20))
  
  g <- ggplot_gtable(ggplot_build(lpt))
  stripr <- rev(which(grepl('strip-t', g$layout$name)))
  fills <- c("blue","red","green","yellow","magenta","purple")
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  
  ggsave(plot = g, paste0("../output/cdrs/", 
  cdrx, "_", cdr_len, "/", cdrx, "_", cdr_len, "_", "logo.svg"), 
  device = "svg", height = 3.33*length(unique(lp_df$DBS)), width = 5)

  return()

}
