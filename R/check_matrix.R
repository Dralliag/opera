#' Function to check and modify the input class and type
#'
#' @param mat \code{data.frame, data.table, tibble}. Object to be cast to matrix.
#' @param name \code{character}. Name of the object to be cast. 
#'
#' @return a 3d array if a 3d array is provided, else a matrix.
#'
check_matrix <- function(mat, name) {
  # check class
  if (! is.null(mat) && ! "array" %in% class(mat)) {
    mat <- tryCatch({
      if (is.vector(mat)) {
        t(as.matrix(mat))
      } else {
        as.matrix(mat)
      }
    }, error = function(e) {
      stop("Error when casting ", name, " to matrix : \n", 
           e$message)
    })
  }
  # check type 
  if (! is.null(mat) && ! typeof(mat) %in% c("integer", "double")) { 
    stop("Type of ", mat, " must be either integer or double")
    
  } else if (! is.null(mat) && typeof(mat) == "integer") {
    storage.mode(mat) <- "numeric"
  }
  
  return(mat)
}