create_sparse_matrix <- function(array){
## Purpose: Create a sparse matrix from a 2D array in R, by saving the indices of the non-zero values
## Input: 2D array 
## Output: Sparse matrix
## NOTE: This function is used to create sparse footprints for STILT specifically,
## therefore, index_x is the column index and index_y is the row index, since this is
## how STILT stores the footprints. 

indices <- which(array != 0, arr.ind = T)
infl <- array[indices]
index_x = indices[, 'col']
index_y = indices[, 'row']    

return(list("infl" = infl, "index_x" = index_x, "index_y" = index_y))
}