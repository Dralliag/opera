mape <-
function(x,y, na.rm=TRUE) {mean(abs(x-y)/y, na.rm=na.rm)}
