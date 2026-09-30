# Utilidades compartidas. Nada aquí consume la semilla del experimento.

fmt <- function(x, d = 3) {
  y <- as.numeric(x)
  out <- trimws(formatC(y, digits = d, format = "f"))
  out[!is.na(y) & y == 0] <- "0"
  out
}

fmt_cor <- function(x) sub("^(-?)0", "\\1", fmt(x, 3))

fmt_sci <- function(x, d = 1) {
  trimws(formatC(as.numeric(x), digits = d, format = "e"))
}

y_lista <- function(x) {
  x <- as.character(x)
  if (!length(x)) return("(ninguna)")
  if (length(x) == 1L) return(x)
  paste(paste(x[-length(x)], collapse = ", "), "y", x[length(x)])
}

rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))

vif_uno <- function(X, j = 1L) {
  y <- X[, j]
  A <- cbind(1, X[, -j, drop = FALSE])
  ajuste <- qr.fitted(qr(A), y)
  r2 <- 1 - sum((y - ajuste)^2) / sum((y - mean(y))^2)
  1 / (1 - r2)
}

ols_coef <- function(X, y) {
  A <- cbind("(Intercept)" = 1, X)
  b <- qr.coef(qr(A), y)
  names(b) <- colnames(A)
  b
}

# Seudoinversa de Moore-Penrose de (1 | X). Es el MCO de norma mínima
# cuando la matriz no tiene rango completo.
pinv_coef <- function(X, y, tol_rel = 1e-8) {
  A <- cbind("(Intercept)" = 1, X)
  des <- svd(A)
  umbral <- tol_rel * max(des$d)
  pos <- des$d > umbral
  coef <- des$v[, pos, drop = FALSE] %*%
    ((t(des$u[, pos, drop = FALSE]) %*% y) / des$d[pos])
  list(
    coef = stats::setNames(drop(coef), colnames(A)),
    rango = sum(pos),
    umbral = umbral
  )
}

predecir_lineal <- function(X, coef) {
  as.numeric(cbind(1, X) %*% coef)
}

rss_de <- function(y, X, coef) sum((y - predecir_lineal(X, coef))^2)

pendientes <- function(b) b[names(b) != "(Intercept)"]

norma_l2 <- function(b) sqrt(sum(pendientes(b)^2))
norma_l1 <- function(b) sum(abs(pendientes(b)))

vivos <- function(b) {
  s <- pendientes(b)
  names(s)[s != 0]
}

coef_glm <- function(ajuste, s = NULL) {
  b <- if (is.null(s)) coef(ajuste) else coef(ajuste, s = s)
  b <- drop(as.matrix(b))
  stats::setNames(as.numeric(b), rownames(as.matrix(coef(ajuste))))
}

autovalor_min <- function(M) {
  min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)
}

# Ridge de libro: X centrada, y centrada, intercepto fuera de la pena.
# beta = (X'X + λ I)^{-1} X'y. λ = 0 recupera MCO si X'X es invertible.
ridge_cerrado <- function(X, y, lambda) {
  xc <- scale(X, center = TRUE, scale = FALSE)
  yc <- y - mean(y)
  xtx <- crossprod(xc)
  b <- solve(xtx + lambda * diag(ncol(X)), crossprod(xc, yc))
  b <- drop(b)
  names(b) <- colnames(X)
  b0 <- mean(y) - sum(colMeans(X) * b)
  c("(Intercept)" = unname(b0), b)
}
