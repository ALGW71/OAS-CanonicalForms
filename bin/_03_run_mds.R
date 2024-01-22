library(dplyr)
library(stringr)
library(data.table)
library(lmds)
library(tidyverse)

files <- dir("../output/dist_mats/", recursive = T, pattern = ".csv")
files <- files[ !grepl("full", files) ]

for (csv in files){
  
  # it reads sample name - so if you append the DTW output with sample name.
  # it should be fine
  ##############################################################################
  fl <- basename(csv)
  CDRX <- paste0(str_extract(csv, pattern = "chain_.") %>%
                   str_remove(pattern = "chain_"),
                 str_extract(csv, pattern = "cdr_.") %>%
                   str_remove(pattern = "cdr_"))
  ##############################################################################
  
  
  ##############################################################################
  print(Sys.time())
  print("XXX - Reading CSV")
  input <- fread(paste0("../output/dist_mats/", csv), nThread=30)
  print(Sys.time())
  print("XXX - CSV file read.")
  gc()
  ##############################################################################
  
  
  ##############################################################################
  print(Sys.time())
  print("XXX - Building Upper tri.")
  file_names <- input$pdb1 %>% unique
  
  # Define chunk size
  chunk_size <- 1000
  num_chunks <- ceiling(length(file_names) / chunk_size)
  
  print(paste0("XXX - Number of chunks: ", num_chunks))
  
  # Iterate over chunks
  for (chunk_idx in 1:num_chunks) {
    cat("XXX - Processing chunk", chunk_idx, "of", num_chunks, "\n")
    
    # Determine start and end indices for the current chunk
    start_idx <- (chunk_idx - 1) * chunk_size + 1
    end_idx <- min(chunk_idx * chunk_size, length(file_names))
    
    # Initialize an empty list to store data.tables for each iteration
    output_list <- vector("list", end_idx - start_idx + 1)
    
    # Use a for loop to iterate over file_names in the current chunk
    for (x in start_idx:end_idx) {
      # Filter the rows with the current file_name and then dcast
      df <- dcast(input[pdb1 == file_names[x]], pdb1 ~ pdb2, value.var = "dist")
      
      # Store the resulting data.table in the list
      output_list[[x - start_idx + 1]] <- df
    }
    
    # Combine the list of data.tables into a single data.table
    output_chunk <- rbindlist(output_list, fill = TRUE)
    
    # Save the output_chunk to a file
    fwrite(output_chunk, paste0("../output/tmp/output_chunk_", chunk_idx, ".csv"), nThread=30)
    
    # Clear variables and perform garbage collection
    remove(output_list, output_chunk, df)
    gc()
  }
  
  remove(input)
  gc()
  print(Sys.time())
  print("XXX - Files read")
  ##############################################################################
  
  dir <- "../output/tmp/"
  if (!dir.exists(dir)) dir.create(dir)
  
  ##############################################################################
  # combine the csv files into mega matrix
  output <- fread(paste0("../output/tmp/output_chunk_", 1, ".csv"), nThread=30)
  output_list <- list(output)
  
  # Iterate over chunks
  if (num_chunks > 1){
    for (chunk_idx in 2:num_chunks) {
      cat("Binding chunk", chunk_idx, "of", num_chunks, "\n")
      ch <- fread(paste0("../output/tmp/output_chunk_", chunk_idx, ".csv"), nThread=30)
      output_list[[chunk_idx]] <- ch
      remove(ch)
      file.remove(paste0("../output/tmp/output_chunk_", chunk_idx, ".csv"))
      gc()
    }
  }
  
  file.remove(paste0("../output/tmp/output_chunk_1.csv"))
  
  output <- rbindlist(output_list, fill = TRUE)
  remove(output_list)
  gc()
  ##############################################################################
  
  # HERE - check they are identical, then work out index to delete
  print(identical(output$pdb1, colnames(output[,c(-1)])))
  #print(output$pdb1[(nrow(output)-10):nrow(output)])
  #print(colnames(output[,c(-1)])[(ncol(output[,c(-1)])-10):ncol(output[,c(-1)])])
  
  ##############################################################################
  print(identical(output$pdb1, colnames(output[,c(-1)])))
  mega_mat <- as.matrix(output[,-1])
  rownames(mega_mat) <- str_remove_all(string = output$pdb1, pattern = ".pdb")
  colnames(mega_mat) <- str_remove_all(string = colnames(mega_mat), pattern = ".pdb")
  dim(mega_mat)
  remove(output)
  gc()
  print(Sys.time())
  print("XXX - Upper tri created.")
  ##############################################################################
  
  
  ##############################################################################
  tm <- t(mega_mat)
  tm[upper.tri(tm)] <- mega_mat[upper.tri(mega_mat)]
  mega_mat <- tm
  remove(tm)
  gc()
  
  missing_logical <- mega_mat == "Missing"
  table(missing_logical)
  # get the rowsums and col sums
  rows <- rowSums(missing_logical)
  cols <- colSums(missing_logical)
  
  # get the median value
  md_rows <- median(rows)
  md_cols <- median(cols)
  
  rows_to_drop <- names(rows[rows > md_rows])
  cols_to_drop <- names(cols[cols > md_cols])
  
  print("Dims before drop missing")
  print(dim(mega_mat))
  mega_mat <- mega_mat[!(rownames(mega_mat) %in% rows_to_drop),
                       !(colnames(mega_mat) %in% cols_to_drop)]
  
  row_names <- rownames(mega_mat)
  col_names <- colnames(mega_mat)
  mega_mat <- matrix(as.numeric(mega_mat), nrow = nrow(mega_mat))
  rownames(mega_mat) <- row_names
  colnames(mega_mat) <- col_names
  
  print("Dims after drop missing")
  print(dim(mega_mat))
  
  print(Sys.time())
  fwrite(mega_mat, paste0("../output/dist_mats/fullmatrix_", fl), nThread=30)
  print("XXX - Full matrix created.")
  ##############################################################################
  
  # At this point you need to get the EPS value. Plot these and make conclusions on the clusters.
  # You want to decide on the 
  
  ##############################################################################
  print("XXX - Starting landmark MDS")
  ldmks <- nrow(mega_mat)
  
  # select all ----------------------------------------------------------------
  landmarks <- mega_mat[1:ldmks, ]
  print(dim(landmarks))
  gc()
  fit <- cmdscale_landmarks(landmarks, rescale = FALSE, ndim = 2)
  print(dim(fit))
  gc()
  vis_mds <- data.frame("PDB" = colnames(mega_mat),"Coord1" = fit[,1],"Coord2" = fit[,2])
  # save this MDS data...
  
  dir <- "../output/mds_dists/"
  if (!dir.exists(dir)) dir.create(dir)
  
  # save according to the sample number
  fwrite(vis_mds, paste0("../output/mds_dists/vis_mds_", fl),
         nThread=30)
  print(Sys.time())
  print("XXX - Finished landmark MDS")
  gc()
  ##############################################################################
  
}
