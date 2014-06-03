mae <-
function (x,y,awake=NULL, na.rm=TRUE) {
  if (is.null(awake)) {
    mean(abs(x-y), na.rm=na.rm)
  }
  else {
    sum(abs(x-y) * awake, na.rm=na.rm) / sum(awake * (!is.na(x-y)), na.rm=na.rm)
  }
}
