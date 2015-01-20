# Testing popStruct code

# basic test
test <- matrix(c(1,0,1,0,1,1,2,1,1,
                 2,0,0,1,0,1,2,0,1,
                 1,1,1,0,1,2,1,0,1), ncol = 3)
popStruct(test, method = "eigen")
popStruct(test, method = "pcaseq")


# test genotype counts > 2
# Without data checking
# pcaseq returns numbers, but this isn't interpretable(?)
# eigenstrat returns Infs
test <- matrix(c(1,0,1,0,1,1,2,1,1,
                 2,0,0,3,0,1,2,0,1,
                 1,1,1,0,5,2,1,0,1), ncol = 3)
popStruct(test, method = "eigen")
popStruct(test, method = "pcaseq")


# test genotype counts < 0
test <- matrix(c(1,0,1,0,1,1,2,1,1,
                 2,0,0,2,0,1,2,0,1,
                 1,1,1,0,1, -2,1,0,1), ncol = 3)
popStruct(test, method = "eigen")
popStruct(test, method = "pcaseq")


# test genotype counts non-integer
test <- matrix(c(1,0,1,0,1,1,2,1,1,
                 2,0,0,2,0,1.7,2,0,1,
                 1,1,1,0,1, 2,1,0,1), ncol = 3)
popStruct(test, method = "eigen")
popStruct(test, method = "pcaseq")

#test missing genotype data
test <- matrix(c(1,0,1,0,1,1,2,1,1,
                 2,0,0,2,0,1,2,0,1,
                 1,1,1,NA,1,2,1,0,1), ncol = 3)
popStruct(test, method = "eigen")
popStruct(test, method = "pcaseq")

#test letters instead of numbers
test <- matrix(c("A",0,1,0,1,1,2,1,1,
                 2,0,0,2,0,1,2,0,1,
                 1,1,1,0,1,2,1,0,1), ncol = 3)
popStruct(test, method = "eigen")
popStruct(test, method = "pcaseq")

# only one subject
test <- matrix(c(1,1,1,0,1,2,1,0), ncol = 1)
popStruct(test, method = "eigen")
popStruct(test, method = "pcaseq")


# only one locus
test <- matrix(c(1,1,1,1,2,1,0,1), ncol = 8)
popStruct(test, method = "eigen")
popStruct(test, method = "pcaseq")


# non-matrix data

# vector
# returns an error from rowMeans
test <- c(0,1,0,1,1,2,1,1)
popStruct(test, method = "eigen")
popStruct(test, method = "pcaseq")

# data.frame
# this returns an error because
# of the test for 0,1,2
test <- matrix(c(2,0,1,0,1,1,2,1,1,
                 2,0,0,2,0,1,2,0,1,
                 1,1,1,1,2,1,0,1,0), ncol = 3)
test <- data.frame(test)
popStruct(test, method = "eigen")
popStruct(test, method = "pcaseq")