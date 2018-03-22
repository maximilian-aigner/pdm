set.seed(123)

# choose 'causal' gene loci
k = 10
causal.genes <- sample(14431249:49588215, k, replace = FALSE)
effect.sizes <- runif(k, min = 1.25, max = 1.5)

