# not sure what if any validity checking needs to be done
# for this object
setClass(
  Class = "popStructure",
  representation = representation(
    grm = "matrix",
    pcs = "matrix",
    eigenvalues = "numeric"))

# this needs work
setMethod(
  f = plot,
  signature = "popStructure",
  definition = function(x, y, ...)
    {
      plot(x@pcs[,1], x@pcs[,2], xlab = "Principle Component 1",
           ylab = "Principle Component 2")
  })

# this needs work
setMethod(
  f = "print",
  signature = "popStructure",
  definition = function(x, ...)
    {
      cat("Results of PCA-seq Population Sturcture Inference \n")
      cat("Principle Components \n")
      print(x@pcs)
      cat("Eigenvalues \n")
      print(x@eigenvalues)
  })

# this needs work
setMethod(
  f = "show",
  signature = "popStructure",
  definition = function(object)
    {
      cat("Results of PCA-seq Population Sturcture Inference \n")
      cat("Principle Components \n")
      print(x@pcs)
      cat("Eigenvalues \n")
      print(x@eigenvalues)
  })