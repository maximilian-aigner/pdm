source('knockoffs_gaussian.R')
library(glmnet)
library(MASS)
library(corrplot)

n <- 500 # observations
p <- 100   # variables available
s <- floor(0.1*p)    # non-nulls (simulation)

# generate data matrix
sim.mu <- rep(0, p)
# sim.sig <- matrix(.2, nrow = p, ncol = p) + diag(p) * .8
sim.sig <- diag(p)
X <- mvrnorm(n, mu = sim.mu, Sigma = sim.sig)

# generate Ys
selected.vars <- sample(1:p, s, replace = FALSE)
true.beta <- matrix(5*rep(1, s))
y <- X[, selected.vars] %*% true.beta

# fit lm to augmented data
aug.X <- knockoffs.g(X)

aug.reg.fit <- cv.glmnet(aug.X, y, alpha = 0)
Z <- abs(coef(aug.reg.fit, lambda = aug.reg.fit$best.lambda))
W <- Z[1:p] - Z[(p+1) : (2*p)]

pdf(file="var_check.pdf")
plot(1:length(W), W)
points(selected.vars, W[selected.vars], col = "red", pch = 3)
dev.off()

pdf(file="var_corrs.pdf")
corrplot(cor(aug.X), method = "circle")
dev.off()
# aug.fit <- lm(y ~ ., as.data.frame(aug.X))
