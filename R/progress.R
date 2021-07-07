init_progress <- function(T) {
  # init progress
  cat("[---------|---------|----------|---------|---------|---------|---------|---------|---------|---------]\r");
  cat("[");
  
  steps = table(ceiling((1:100)*T/100))
  
  return(steps)
}

update_progress <- function(t, steps) {
  if (t %in% names(steps))
    cat(paste0(rep("=", steps[[as.character(t)]]), collapse = ""));
}

end_progress <- function() {
  cat("]\n");
}