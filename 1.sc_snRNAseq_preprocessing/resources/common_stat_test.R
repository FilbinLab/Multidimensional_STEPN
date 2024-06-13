fisher_test <- function(a, b, bg){
  a = unlist(a); b = unlist(b)
  x = length(intersect(a,b))
  m = length(setdiff(a,b))
  n = length(setdiff(b,a))
  k = length(setdiff(bg, union(a,b)))
  mat = matrix(c(x,m,n,k),2,2)
  print(mat)
  return(fisher.test(mat, alternative="greater")$p.value)
}

jaccard_index <- function(a, b){
    a = unlist(a); b=unlist(b)
    return(length(intersect(a,b))/length(union(a,b)))
}