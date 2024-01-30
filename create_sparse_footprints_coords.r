create_sparse_matrix_coords <- function(array, lons, lats){
## Purpose: Create a sparse matrix from a 2D array in R, by saving the longitude and latitude coordinates of each of the non-zero values
## Input: 2D array 
## Output: Sparse matrix
## NOTE: This function is used to create sparse footprints for STILT specifically,
## therefore, the longitudes are reprsented by the columns and the latitues by the rows, since this is
## how STILT stores the footprints. 

indices <- which(array != 0, arr.ind = T)
infl <- array[indices]
lon = lons[indices[, 'col']]
lat = lats[indices[, 'row']]

return(list("infl" = infl, "lon" = lon, "lat" = lat))
}