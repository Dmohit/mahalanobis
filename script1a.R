library(MASS)
library(dplyr)
library(ggplot2)
library(ConfidenceEllipse)
library(scales)

Sigma <- matrix(c(2,0.5,0.5,2), nrow = 2)
mu <- c(1, 4)
set.seed(3400)

pt_at_dist <- function(origin, slope, d)
	return(c(origin[1] + d / sqrt(1 + slope ^ 2), origin[2] + d * slope / sqrt(1 + slope ^ 2)))

get_intercept <- function(x1, y1, slope)
{
	# Point slope formula (y - y1) = slope(x - x1)
    y_intercept = slope * (- x1) + y1
    return(y_intercept)
}

X <- mvrnorm(200, mu, Sigma)
df = data.frame(x = X[,1], y = X[,2])
head(df)

ex_point <- c(-2.8, 5)

Sx = cov(X)
mm <- colMeans(X)

rows <- c(3, 22, 30, 35, 71, 77, 153, 114, 127)

ellipse_95 <- confidence_ellipse(df, x, y, conf_level = 0.95)

pch_vec <- rep(16, nrow(X))
pch_vec[rows] <- 3

cex_vec <- rep(0.7, nrow(X))
cex_vec[rows] <- 1.2

col_vec <- rep('blue', nrow(X))
col_vec[rows] <- alpha('steelblue', 0.7)

plot(X, pch = pch_vec, cex = cex_vec, col = col_vec,
	xlim = c(-5,6), ylim = c(-0.7, 9), xlab = 'X', ylab = 'Y')
lines(ellipse_95, lty = 'dashed', col = 'steelblue')
points(x = ex_point[1], y = ex_point[2], pch = 4, cex = 1.05)
grid(10, 10)
legend('bottomleft', lty = 'dashed', col = 'steelblue', '95% Confidence Ellipse', bty = 'n')
title(main = 'Scatterplot')

subset <- data.frame(X[rows,], row.names = rows)
subset <- rbind(subset, ex_point)
row.names(subset)[10] <- 'ex_pt'

#down, left, up right :: 1,2,3,4
text(subset, row.names(subset), pos = c(1, 1, 2, 2, 2, 2, 2, 3, 2, 2), col = 'steelblue')

points(x = mm[1], y = mm[2], pch = 17, cex = 2.4, col = alpha('red', 0.5))

eigv <- eigen(Sx)$vectors
eigc <- eigen(Sx)$values

slope <- eigv[2,] / eigv[1,]
#abline(get_intercept(x1 = mm[1], y1 = mm[2], slope = slope[1]), slope[1], lty = 'dashed', col = 'steelblue')
#abline(get_intercept(x1 = mm[1], y1 = mm[2], slope = slope[2]), slope[2], lty = 'dashed', col = 'steelblue')

arrowh1 <- pt_at_dist(mm, slope[1], eigc[1])
arrowh2 <- pt_at_dist(mm, slope[2], eigc[2])

arrows(mm[1], mm[2], arrowh1[1], arrowh1[2], col = 'steelblue', angle = 25, length = eigc[1] * 0.1)
arrows(mm[1], mm[2], arrowh2[1], arrowh2[2], col = 'steelblue', angle = 25, length = eigc[1] * 0.1)

mh <- mahalanobis(subset, mm, Sx)

euc_dist <- function(x1, x2) sum((x1 - x2) ^ 2)

eu <- apply(subset, MARGIN = 1, FUN = function(x) euc_dist(x, mm))

key = order(mh)

dev.new()
par(mfrow = c(2, 1))
barplot(mh[key], names.arg = row.names(subset)[key], horiz = FALSE, col = 'steelblue', ylim = c(0, 16),
			ylab = 'Mahalanobis Distance', xlab = 'PointID', main = 'Squared Mahalanobis Distance')
barplot(eu[key], names.arg = row.names(subset)[key], horiz = FALSE, col = 'steelblue', ylim = c(0, 16),
			ylab = 'Euclidean Distance', xlab = 'PointID', main = 'Squared Euclidean Distance')

pt1 <- pt_at_dist(mm, slope[1], 6)
pt2 <- pt_at_dist(mm, slope[2], 6)

eu2 <- c(euc_dist(mm, pt1), euc_dist(mm, pt2))
mh2 <- c(mahalanobis(pt1, mm, Sx), mahalanobis(pt2, mm, Sx))

key2 = order(mh2)

dev.new()
par(mfrow = c(2, 1))
barplot(mh2[key2], names.arg = c(1,2)[key2], horiz = FALSE, col = 'steelblue', ylim = c(0, 37),
			ylab = 'Mahalanobis Distance', xlab = 'PointID', main = 'Squared Mahalanobis Distance')
barplot(eu2[key2], names.arg = c(1,2)[key2], horiz = FALSE, col = 'steelblue', ylim = c(0, 37),
			ylab = 'Euclidean Distance', xlab = 'PointID', main = 'Squared Euclidean Distance')

dev.new()
par(mfrow = c(2, 1))
barplot(mh[key], names.arg = row.names(subset)[key], horiz = FALSE, col = 'steelblue', ylim = c(0, 16),
			ylab = 'Mahalanobis Distance', xlab = 'PointID', main = 'Squared Mahalanobis Distance')
barplot(eu[key], names.arg = row.names(subset)[key], horiz = FALSE, col = 'steelblue', ylim = c(0, 16),
			ylab = 'Euclidean Distance', xlab = 'PointID', main = 'Squared Euclidean Distance')

dev.new()
plot(mm[1], mm[2], pch = 17, cex = 2.4, col = alpha('red', 0.5),
		xlim = c(-5,6.9), ylim = c(-5.6, 9), xlab = 'X', ylab = 'Y')
lines(ellipse_95, lty = 'dashed', col = 'steelblue')
arrows(mm[1], mm[2], arrowh1[1], arrowh1[2], col = 'steelblue')
arrows(mm[1], mm[2], arrowh2[1], arrowh2[2], col = 'steelblue')
points(pt1[1], pt1[2], pch = 4, cex = 1.05)
points(pt2[1], pt2[2], pch = 4, cex = 1.05)
