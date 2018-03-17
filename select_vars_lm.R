source('knockoffs_gaussian.R')
library(glmnet)

n <- 5000 # observations
p <- 20   # variables available
s <- 4    # non-nulls (simulation)

# generate data matrix
X <- matrix(rnorm(n*p, sd = 3), nrow = n, ncol = p)

# generate Ys
selected.vars <- sample(1:p, s, replace = FALSE)
# true.beta <- matrix(abs(3*rt(s, 1)))
true.beta <- matrix(5*rep(1, s))
y <- X[, selected.vars] %*% true.beta
# fit lm to augmented data
aug.X <- knockoffs.g(X)

aug.reg.fit <- cv.glmnet(aug.X, y, alpha = 0)
Z <- abs(coef(aug.reg.fit, lambda = aug.reg.fit$best.lambda))
W <- Z[1:p] - Z[(p+1) : (2*p)]
print(W)

pdf(file="var_check.pdf")
plot(1:length(W), W)
points(selected.vars, W[selected.vars], col = "red", pch = 3)
dev.off()


library(corrplot)
pdf(file="var_corrs.pdf")
corrplot(cor(aug.X), method = "circle")
dev.off()
# aug.fit <- lm(y ~ ., as.data.frame(aug.X))
