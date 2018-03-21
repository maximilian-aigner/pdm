source('knockoffs_gaussian.R')
library(glmnet)
library(MASS)
library(corrplot)

n <- 3000 # observations
p <- 100   # variables available
s <- 15   # non-nulls (simulation)

# generate data matrix
sim.mu <- rep(0, p)
sim.sig <- matrix(.8, nrow = p, ncol = p) + diag(p) * .2
#sim.sig <- diag(p)
X <- mvrnorm(n, mu = sim.mu, Sigma = sim.sig)
X <- scale(X)

# generate Ys
selected.vars <- sample(1:p, s, replace = FALSE)
true.beta <- matrix(1*rep(1, s))
y <- mvrnorm(1, mu = X[, selected.vars] %*% true.beta, Sigma = diag(n))

# fit lm to augmented data
aug.X <- knockoffs.g(X)

aug.reg.fit <- cv.glmnet(aug.X, y, alpha = 0, intercept = FALSE)
Z <- coef(aug.reg.fit, lambda = aug.reg.fit$best.lambda)
W <- abs(Z[2:(p+1)]) - abs(Z[(p+2) : (2*p+1)])

pdf(file="var_check.pdf")
plot(1:length(W), W)
points(selected.vars, W[selected.vars], col = "red", pch = 3)
dev.off()

# pdf(file="var_corrs.pdf")
# corrplot(cor(aug.X), method = "circle")
# dev.off()

# compute threshold
threshold <- knockoff.threshold(W, 0.1)
falsely_discovered <- ((abs(W) > threshold) & (1:length(W) %in% selected.vars))
empirical_fdr <- sum(falsely_discovered) / length(W)

