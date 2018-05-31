# Reverse caret package preprocess function (only for center and scaling)

unPreProc <- function(preProc, data){
  stopifnot(class(preProc) == "preProcess")
  stopifnot("data.frame" %in% class(data))
  for(i in names(preProc$mean)){
    tmp <- data[, i] * preProc$std[[i]] + preProc$mean[[i]]
    data[, i] <- tmp
  }
  return(data)  
}
