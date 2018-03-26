source('knockoffs_gaussian.R')
source('threshold.R')
library(glmnet)
library(MASS)
#library(corrplot)

n <- 3000 # observations
p <- 100   # variables available
s <- 20   # non-nulls, excluding interactions
k <- 5 # interactions among non-nulls (k <= s)

# generate data matrix
sim.mu <- rep(0, p)
# sim.sig <- matrix(.4, nrow = p, ncol = p) + diag(p) * .6
A = matrix(runif(p^2)*2-1, ncol = p)
sim.sig <- t(A) %*% A
y.var <- 1.5
#sim.sig <- diag(p)
X <- mvrnorm(n, mu = sim.mu, Sigma = sim.sig)
colnames(X) <- sapply(1:p, function(i) paste("X", i, sep="_"))
X <- scale(X)

# generate Ys
selected.vars <- sample(1:p, s, replace = FALSE)
# true.beta <- matrix(1*rep(1, s))
# y <- mvrnorm(1, mu = X[, selected.vars] %*% true.beta, Sigma = diag(n)*y.var)

interaction.vars <- sample(selected.vars, size = 2*k, replace = FALSE)
# first k multiply the last k to form k interaction variables
interactions <- X[, interaction.vars[1:k]] * X[, interaction.vars[(k+1):(2*k)]]
colnames(interactions) <- paste(colnames(X)[interaction.vars[1:k]],
                                colnames(X)[interaction.vars[(k+1):(2*k)]], sep = ".")
all.selected <- c(selected.vars, (p+1):(p+k))
interaction.X <- cbind(X, interactions) # n x (s + k) matrix
true.beta <- matrix(c(3.5 * rep(1, s), (3.5)^2 * rep(1, k)))
true.beta <- mvrnorm(1, mu = true.beta, Sigma = diag(s+k)*1)
y <- mvrnorm(1, mu = interaction.X[, all.selected] %*% true.beta, Sigma = diag(n)*y.var)

# fit lm to augmented data
aug.X <- knockoffs.g(X, shrinkage = TRUE)

aug.reg.fit <- cv.glmnet(aug.X, y, alpha = 0, intercept = FALSE)
Z <- coef(aug.reg.fit, lambda = aug.reg.fit$lambda.min)
W <- abs(Z[2:(p+1)]) - abs(Z[(p+2) : (2*p+1)])
names(W) <- colnames(X)
# W <- abs(Z[2:(p+k+1)]) - abs(Z[(p+k+2):(2*p + 2*k + 1)])
# names(W) <- colnames(interaction.X)

# pdf(file="var_corrs.pdf")
# corrplot(cor(aug.X), method = "circle")
# dev.off()

# compute threshold
threshold <- knockoff.threshold(W, 0.1)
falsely_discovered <- ((abs(W) > threshold) & (1:length(W) %in% selected.vars))
empirical_fdr <- sum(falsely_discovered) / length(W)

par(mfrow = c(1, 2))
#pdf(file="var_check.pdf")
plot(1:length(W), W)
#points(all.selected, W[all.selected], col = "red", pch = 3)
points(selected.vars, W[selected.vars], col = "red", pch = 3)
abline(h=threshold)
#dev.off()

# Compare with direct fit
direct.fit <- cv.glmnet(X, y, alpha = 0, intercept = FALSE)
direct.beta <- as(coef(direct.fit, lambda = direct.fit$lambda.1se), "numeric")[2:(p+1)]
names(direct.beta) <- colnames(X)
plot(direct.beta)
points(selected.vars, direct.beta[selected.vars], col = "red", pch = 3)



