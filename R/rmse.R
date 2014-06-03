rmse <-
function (x,y, awake=NULL, na.rm=TRUE) 
{
  if (is.null(awake)) {
    sqrt(mean((x-y)^2, na.rm=na.rm))
  } else {
    sqrt(sum((x-y)^2 * awake, na.rm=na.rm) / sum(awake * (!is.na(x-y))))
  }
}
