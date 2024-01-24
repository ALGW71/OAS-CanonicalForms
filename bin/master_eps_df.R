library(dplyr)

master_df <- tibble(

  csv = c(
    "../output/knn_dbs/eps_knn_H_1_length_8.csv", "../output/knn_dbs/eps_knn_H_1_length_9.csv", "../output/knn_dbs/eps_knn_H_1_length_10.csv",
    "../output/knn_dbs/eps_knn_H_2_length_7.csv", "../output/knn_dbs/eps_knn_H_2_length_8.csv", "../output/knn_dbs/eps_knn_H_2_length_10.csv",

    "../output/knn_dbs/eps_knn_L_1_length_6.csv", "../output/knn_dbs/eps_knn_L_1_length_7.csv", "../output/knn_dbs/eps_knn_L_1_length_8.csv",
    "../output/knn_dbs/eps_knn_L_1_length_9.csv", "../output/knn_dbs/eps_knn_L_1_length_11.csv", "../output/knn_dbs/eps_knn_L_1_length_12.csv",
    "../output/knn_dbs/eps_knn_L_2_length_3.csv",
    "../output/knn_dbs/eps_knn_L_3_length_8.csv", "../output/knn_dbs/eps_knn_L_3_length_9.csv", 
    "../output/knn_dbs/eps_knn_L_3_length_10.csv", "../output/knn_dbs/eps_knn_L_3_length_11.csv",
    
    "../output/knn_dbs/eps_knn_H_3_length_11.csv", "../output/knn_dbs/eps_knn_H_3_length_12.csv", "../output/knn_dbs/eps_knn_H_3_length_13.csv",
    "../output/knn_dbs/eps_knn_H_3_length_14.csv", "../output/knn_dbs/eps_knn_H_3_length_15.csv", "../output/knn_dbs/eps_knn_H_3_length_16.csv",
    "../output/knn_dbs/eps_knn_H_3_length_17.csv", "../output/knn_dbs/eps_knn_H_3_length_18.csv", "../output/knn_dbs/eps_knn_H_3_length_19.csv"),

  dist = c(
    "../output/dist_mats/fullmatrix_H_1_length_8.csv", "../output/dist_mats/fullmatrix_H_1_length_9.csv", "../output/dist_mats/fullmatrix_H_1_length_10.csv",
    "../output/dist_mats/fullmatrix_H_2_length_7.csv", "../output/dist_mats/fullmatrix_H_2_length_8.csv", "../output/dist_mats/fullmatrix_H_2_length_10.csv",
    
    "../output/dist_mats/fullmatrix_L_1_length_6.csv", "../output/dist_mats/fullmatrix_L_1_length_7.csv", "../output/dist_mats/fullmatrix_L_1_length_8.csv",
    "../output/dist_mats/fullmatrix_L_1_length_9.csv", "../output/dist_mats/fullmatrix_L_1_length_11.csv", "../output/dist_mats/fullmatrix_L_1_length_12.csv",
    "../output/dist_mats/fullmatrix_L_2_length_3.csv",
    "../output/dist_mats/fullmatrix_L_3_length_8.csv", "../output/dist_mats/fullmatrix_L_3_length_9.csv", 
    "../output/dist_mats/fullmatrix_L_3_length_10.csv", "../output/dist_mats/fullmatrix_L_3_length_11.csv",
    
    "../output/dist_mats/fullmatrix_H_3_length_11.csv", "../output/dist_mats/fullmatrix_H_3_length_12.csv", "../output/dist_mats/fullmatrix_H_3_length_13.csv",
    "../output/dist_mats/fullmatrix_H_3_length_14.csv", "../output/dist_mats/fullmatrix_H_3_length_15.csv", "../output/dist_mats/fullmatrix_H_3_length_16.csv",
    "../output/dist_mats/fullmatrix_H_3_length_17.csv", "../output/dist_mats/fullmatrix_H_3_length_18.csv", "../output/dist_mats/fullmatrix_H_3_length_19.csv"), 

  dunbrack_length = c(
    13, 14, 15, 9, 10, 12,
    11, 12, 13, 14, 16, 17,
    8,
    8, 9, 10, 11,
    11,12,13,14,15,16,17,18,19),

  chain = c("H","H","H","H","H","H",
            "L","L","L","L","L","L",
            "L",
            "L","L","L","L",
            "H","H","H", "H","H","H", "H","H","H"),

  cdrx = c("H1","H1","H1","H2","H2","H2",
           "L1","L1","L1","L1","L1","L1",
           "L2",
           "L3","L3","L3","L3",
           "H3","H3","H3", "H3","H3","H3", "H3","H3","H3"),
          
  cdr_len = c(8,9,10,
              7,8,10,
              6,7,8,
              9,11,12,
              3,
              8,9,10,11,
              11, 12, 13, 14, 15, 16, 17, 18, 19))

# If only running a specific file
# files <- "eps_knn_H_2_length_8.csv"
# master_df <- master_df %>% filter(grepl(pattern = files, x =master_df$csv))

