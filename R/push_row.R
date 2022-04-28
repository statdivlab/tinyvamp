push_row <- function(row_vec,
                     row_num,
                     target_matrix){
  target_matrix[row_num,] <- row_vec
  return(target_matrix)
}
