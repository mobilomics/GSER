args <- commandArgs(trailingOnly=TRUE)
filename <- args[1]
samplesize <- as.integer(args[2])
library(ShortRead)
sampler <- FastqSampler(filename,samplesize)
fq <- yield(sampler)
samplemean <- mean(width(fq))

cat(filename,samplemean,sep="\t")
cat("\n")
