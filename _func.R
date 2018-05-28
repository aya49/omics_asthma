## helper functions
## aya43@sfu.ca
## created 20180517
## last modified 20180517



## input: Sys.time() value
## output: time elapsed since input
time_output = function(start, tz="GMT") {
  start = as.POSIXct(start)
  end = Sys.time()
  time_elapsed = difftime(end, start, units="secs")
  cat(ft(start,tz=tz), "-", ft(end,tz=tz), ">", ft(time_elapsed,tz=tz))
}

## input: Sys.time() value
## output: formatted time as string; used in time_output function
ft = function(time, tz="GMT") {
  return( format(.POSIXct(time,tz=tz), "%H:%M:%S") )
}

## input: x=vector of indices; n=cores or number of parts you want to split x into
## ouput: list of n indices vector
loop_ind_f <- function(x,n) {
  if (n==1) return(list(x))
  return(split(x, ceiling(seq_along(x)/ceiling(length(x)/n))))
}

## input: matrix
## output: returns u1 (col index with only 1 unique element), ua (col index where every row is a unique element), prints colnames and its unique elements if there is less than n unique elements
col_probe = function(m,n=15) {
  u1 = NULL
  ua = NULL
  nm = nrow(m)
  for (col in 1:ncol(m)) {
    coln = colnames(m)[col]
    if (is.data.table(m)) {
      a = unique(as.matrix(m[,..col]))
    } else {
      a = as.vector(unique(as.matrix(m[,col])))
    }
    la = length(a)
    if (la == 1) u1 = append(u1,col)
    if (la == nm) ua = append(ua,col)
    
    cat("\n", coln, ": ", la)
    if (la<n) cat("; ",a)
  }
  return(list(u1=u1,ua=ua))
}






## input: file path and a file extension 
## output: List of all file names in given path with the given extension
fileNames = function(pathList, ext="fcs") {
    temp.str = sapply(str_split(pathList, "/"), function(x) x[length(x)])
    pathList = sub(paste0(".", ext,"$"), "", temp.str, ignore.case=T)
  return(pathList)
}
