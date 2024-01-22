library(tidyverse)
library(dbscan)
library(data.table)
library(grid)
library(gridExtra)
library(png)

plot_loops <- function(csv){
  
  fl <- basename(csv)
  
  CDRX <- paste0(str_extract(csv,
                             pattern = "mds_.") %>%
                   str_remove(pattern = "mds_"),
                 str_extract(csv,
                             pattern = "[1-9]_length") %>%
                   str_remove(pattern = "_length"))
  
  chain <- str_extract(CDRX, pattern = "^.")
  
  cdr_id <- str_extract(CDRX, pattern = ".$") %>%
    paste0("cdr", . , "_len", collapse = "")
  
  ## FOR THIS TO WORK IT MUST BE UNDERSCORE SAMPLE AT THE END...................
  cdr_length <- str_extract(csv, pattern = "length_[0-9]{1,2}") %>%
    str_remove("length_") %>%
    str_remove("_") %>%
    as.integer
  
  # Plot
  immb_loops <- dir(paste0("../output/aligned/pics/",
                           "chain", chain, "_", cdr_id, "gth", cdr_length,
                           "/ImmuneBuilder/"), full.names = T)
  
  sabdab_loops <- dir(paste0("../output/aligned/pics/",
                             "chain", chain, "_", cdr_id, "gth", cdr_length,
                             "/SAbDab-Dubrack/"), full.names = T)
  
  draw_plots <- function(x){
    rl <- lapply(x, png::readPNG)
    gl <- lapply(rl, grid::rasterGrob)
    gtit <- deparse(substitute(x)) %>% toupper
    gdp <- do.call(gridExtra::grid.arrange, c(grobs = gl, ncol=3, top = gtit))


    ggsave(plot = gdp,
           filename = paste0("../output/aligned/loop_plots/",
                             deparse(substitute(x)),
                             ".svg"),
           device = "svg", height = 22, width = 6, dpi = 300)
    dev.off()
  }
  
  draw_plots(immb_loops)
  draw_plots(sabdab_loops)
  
}
