dim.b <- 10

colors <- brewer.pal(n = 5, name = "Set1")

base_b <- rep(1, dim.b)
mults_exp <- seq(-10, 1, length.out = 200)
mults_pos <- 2^mults_exp
mults <- c(-rev(mults_pos), mults_pos)
#mults <- seq(-2, 2, length.out = 100)

gamma = 2

losses <- list(
  #none <- list(f = function(x, l) { return(0); }, name = "No penalisation"),
  #ridge <- list(f = function(x, l) { return(l * x^2); }, color=colors[1], name = "Ridge penalty", lty=2),
  lasso <- list(f = function(x, l) { return(l * abs(x)); }, color=colors[1], name = "Lasso",lty=3),
  mcp <- list(f = function(x, l) { 
      if(abs(x) < gamma*l)
        return(l * abs(x) - x^2 / (2*gamma))
      else
        return(1/2 * gamma * l^2);
  }, name = "MCP", color = colors[2], lty=4),
  scad <- list(
    f = function(x, l) {
      if (abs(x)<l) 
        return(l * abs(x))
      else if (l < abs(x) && abs(x) < gamma*l)
        return((2*gamma*l*abs(x) - x^2 - l^2)/(2*(gamma-1)))
      else
        return(l^2 * (gamma+1)/2)
    }, name = "SCAD", color=colors[3], lty=5)
)

library(extrafont)
loadfonts()
pdf('figures/losses.pdf')
par(cex.main = 1.5, cex.lab = 1.5)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(mults, rep(0, length(mults)), type = 'l', lty = 1, ylim = c(0, 2), xlab = expression(beta), ylab = expression(loss(beta)),
     main = "Comparison of penalties")
lapply(losses, function(loss) {
  y <- sapply(mults, function(t) loss$f(t, 1))
  lines(mults, y, col = loss$color)
})
legend(-2, .5, legend = c("Lasso", "MCP", "SCAD"), col = sapply(losses, function(l) l$color), lty = 1)
dev.off()

