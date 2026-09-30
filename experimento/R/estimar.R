# Ajustes que usa la bitácora: MCO, Ridge, LASSO y la fórmula cerrada.

ajustar_cv <- function(X, y, alpha, pliegues, standardize = TRUE) {
  cv.glmnet(
    X, y,
    alpha = alpha,
    foldid = pliegues,
    standardize = standardize,
    intercept = TRUE,
    family = "gaussian"
  )
}

resumen_ajuste <- function(nombre, coef, ytr, Xtr, yte, Xte, lambda = NA_real_) {
  list(
    nombre = nombre,
    coef = coef,
    lambda = unname(lambda),
    rmse = rmse(yte, predecir_lineal(Xte, coef)),
    rss = rss_de(ytr, Xtr, coef),
    l1 = norma_l1(coef),
    l2 = norma_l2(coef),
    vivos = vivos(coef),
    intercepto = unname(coef[["(Intercept)"]])
  )
}

ajustar_escenario <- function(datos) {
  p <- partir(datos)
  b_ols <- ols_coef(p$Xtr, p$ytr)
  cv_ridge <- ajustar_cv(p$Xtr, p$ytr, alpha = 0, p$pliegues)
  cv_lasso <- ajustar_cv(p$Xtr, p$ytr, alpha = 1, p$pliegues)
  b_ridge <- coef_glm(cv_ridge, s = "lambda.min")
  b_lasso <- coef_glm(cv_lasso, s = "lambda.min")
  b_1se <- coef_glm(cv_lasso, s = "lambda.1se")
  ajuste_1 <- glmnet(
    p$Xtr, p$ytr,
    alpha = 1, lambda = 1,
    standardize = TRUE, intercept = TRUE, family = "gaussian"
  )
  b_lam1 <- coef_glm(ajuste_1)

  modelos <- list(
    mco = resumen_ajuste("MCO", b_ols, p$ytr, p$Xtr, p$yte, p$Xte),
    ridge = resumen_ajuste(
      "Ridge", b_ridge, p$ytr, p$Xtr, p$yte, p$Xte, cv_ridge$lambda.min
    ),
    lasso = resumen_ajuste(
      "LASSO", b_lasso, p$ytr, p$Xtr, p$yte, p$Xte, cv_lasso$lambda.min
    ),
    lasso_1se = resumen_ajuste(
      "LASSO_1se", b_1se, p$ytr, p$Xtr, p$yte, p$Xte, cv_lasso$lambda.1se
    ),
    lasso_lambda1 = resumen_ajuste(
      "LASSO_lambda1", b_lam1, p$ytr, p$Xtr, p$yte, p$Xte, 1
    )
  )

  xc <- scale(p$Xtr, center = TRUE, scale = FALSE)
  xtx <- crossprod(xc)
  list(
    datos = datos,
    partes = p,
    modelos = modelos,
    cv_ridge = cv_ridge,
    cv_lasso = cv_lasso,
    cor12 = cor(p$Xtr[, "x01"], p$Xtr[, "x02"]),
    cor13 = cor(p$Xtr[, "x01"], p$Xtr[, "x03"]),
    vif = vif_uno(p$Xtr, 1L),
    kappa = kappa(xtx, exact = TRUE),
    autovalor_min = autovalor_min(xtx),
    rango_aumentada = qr(cbind(1, p$Xtr))$rank
  )
}

# Mismas columnas, otras unidades: x05 * 1000, x08 / 100.
# Cada versión elige su propio λmín por validación cruzada de cinco pliegues.
ajustar_escalas <- function(escenario) {
  p <- escenario$partes
  X <- p$Xtr
  X[, "x05"] <- X[, "x05"] * 1000
  X[, "x08"] <- X[, "x08"] / 100
  cv_si <- ajustar_cv(X, p$ytr, alpha = 1, p$pliegues, standardize = TRUE)
  cv_no <- ajustar_cv(X, p$ytr, alpha = 1, p$pliegues, standardize = FALSE)
  list(
    X = X,
    estandarizado = list(
      cv = cv_si,
      coef = coef_glm(cv_si, s = "lambda.min"),
      lambda = cv_si$lambda.min
    ),
    crudo = list(
      cv = cv_no,
      coef = coef_glm(cv_no, s = "lambda.min"),
      lambda = cv_no$lambda.min
    )
  )
}

ajustar_pn <- function(datos) {
  p <- partir(datos)
  rango <- qr(cbind(1, p$Xtr))$rank
  pinv <- pinv_coef(p$Xtr, p$ytr)
  df <- data.frame(p$Xtr)
  df$y <- p$ytr
  avisos <- character()
  ajuste_lm <- withCallingHandlers(
    lm(y ~ ., data = df),
    warning = function(w) {
      avisos <<- c(avisos, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  texto_resumen <- paste(capture.output(summary(ajuste_lm)), collapse = "\n")
  aviso_prediccion <- character()
  pred_lm <- withCallingHandlers(
    as.numeric(predict(ajuste_lm, newdata = data.frame(p$Xte), rankdeficient = "simple")),
    warning = function(w) {
      aviso_prediccion <<- c(aviso_prediccion, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  b_lm <- coef(ajuste_lm)
  b_lm[is.na(b_lm)] <- 0
  cv_ridge <- ajustar_cv(p$Xtr, p$ytr, alpha = 0, p$pliegues)
  cv_lasso <- ajustar_cv(p$Xtr, p$ytr, alpha = 1, p$pliegues)
  b_ridge <- coef_glm(cv_ridge, s = "lambda.min")
  b_lasso <- coef_glm(cv_lasso, s = "lambda.min")
  xc <- scale(p$Xtr, center = TRUE, scale = FALSE)
  xtx <- crossprod(xc)
  list(
    datos = datos,
    partes = p,
    rango = rango,
    n_parametros = ncol(p$Xtr) + 1L,
    n_indefinidos = sum(is.na(coef(ajuste_lm))),
    texto_resumen = texto_resumen,
    avisos_lm = unique(avisos),
    aviso_prediccion = unique(aviso_prediccion),
    rss_lm = sum(residuals(ajuste_lm)^2),
    rmse_lm = rmse(p$yte, pred_lm),
    pinv = pinv,
    rmse_pinv = rmse(p$yte, predecir_lineal(p$Xte, pinv$coef)),
    rss_pinv = rss_de(p$ytr, p$Xtr, pinv$coef),
    cv_ridge = cv_ridge,
    cv_lasso = cv_lasso,
    ridge = resumen_ajuste(
      "Ridge", b_ridge, p$ytr, p$Xtr, p$yte, p$Xte, cv_ridge$lambda.min
    ),
    lasso = resumen_ajuste(
      "LASSO", b_lasso, p$ytr, p$Xtr, p$yte, p$Xte, cv_lasso$lambda.min
    ),
    autovalor_min = autovalor_min(xtx),
    autovalor_min_penalizado = autovalor_min(xtx + diag(ncol(p$Xtr)))
  )
}

# Misma semilla y mismas extracciones; solo cambia la desviación del bloque.
sensibilidad_sd <- function(sds = c(0.16, 0.20, 0.40)) {
  filas <- lapply(sds, function(sd) {
    esc <- ajustar_escenario(fabricar(80L, 20L, 55L, sd_idio = sd, sd_trampa = sd))
    b <- esc$modelos$lasso$coef
    data.frame(
      sd = sd,
      cor = esc$cor12,
      vif = esc$vif,
      rmse_mco = esc$modelos$mco$rmse,
      rmse_lasso = esc$modelos$lasso$rmse,
      b_x02 = unname(b[["x02"]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, filas)
}

# La fórmula del capítulo, en el entrenamiento de n > p.
ridge_de_libro <- function(escenario, lambdas = c(0, 1, 10, 100)) {
  p <- escenario$partes
  ols <- ols_coef(p$Xtr, p$ytr)
  filas <- lapply(lambdas, function(lambda) {
    b <- ridge_cerrado(p$Xtr, p$ytr, lambda)
    data.frame(
      lambda = lambda,
      x01 = unname(b[["x01"]]),
      x02 = unname(b[["x02"]]),
      x03 = unname(b[["x03"]]),
      x05 = unname(b[["x05"]]),
      x08 = unname(b[["x08"]]),
      max_abs_dif_mco = if (lambda == 0) max(abs(b - ols)) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, filas)
}
