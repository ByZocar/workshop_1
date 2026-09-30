# Diseño sintético.
#
# La semilla es 20260908. La desviación idiosincrática 0.16 no es un
# resultado: es la decisión que deja la correlación poblacional del
# bloque cerca de .975. Cada escenario vuelve a fijar la semilla, así
# que se puede correr solo.

SEMILLA <- 20260908L
SIGMA <- 1.2
SD_IDIO <- 0.16

beta_verdadero <- function(p) {
  b <- rep(0, p)
  names(b) <- sprintf("x%02d", seq_len(p))
  b[c("x01", "x02", "x05", "x08")] <- c(3, -2, 4, 1.5)
  b
}

# x01, x02 y x03 comparten la latente z.
# x0j = z + sd * e_j, con e_j normales estándar.
# El resto de columnas es ruido N(0, 1), independiente.
# y = X beta + sigma * eps.
fabricar <- function(n, p, n_train,
                     sd_idio = SD_IDIO,
                     sd_trampa = sd_idio,
                     sigma = SIGMA,
                     semilla = SEMILLA) {
  if (n_train >= n) stop("n_train tiene que ser menor que n")
  if (p < 8L) stop("p tiene que alcanzar a x08")
  set.seed(semilla)
  z <- rnorm(n)
  x01 <- z + sd_idio * rnorm(n)
  x02 <- z + sd_idio * rnorm(n)
  x03 <- z + sd_trampa * rnorm(n)
  resto <- matrix(rnorm(n * (p - 3L)), n, p - 3L)
  X <- cbind(x01, x02, x03, resto)
  colnames(X) <- sprintf("x%02d", seq_len(p))
  beta <- beta_verdadero(p)
  y <- as.numeric(X %*% beta + sigma * rnorm(n))
  entrena <- sample.int(n, n_train)
  pliegues <- sample(rep(seq_len(5L), length.out = n_train))
  list(
    X = X,
    y = y,
    beta = beta,
    entrena = entrena,
    prueba = setdiff(seq_len(n), entrena),
    pliegues = pliegues,
    n = n,
    p = p,
    n_train = n_train,
    sd_idio = sd_idio,
    sd_trampa = sd_trampa,
    sigma = sigma,
    semilla = semilla
  )
}

partir <- function(datos) {
  list(
    Xtr = datos$X[datos$entrena, , drop = FALSE],
    ytr = datos$y[datos$entrena],
    Xte = datos$X[datos$prueba, , drop = FALSE],
    yte = datos$y[datos$prueba],
    pliegues = datos$pliegues,
    beta = datos$beta
  )
}
