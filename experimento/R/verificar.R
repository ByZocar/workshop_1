# El informe afirma estas cosas en prosa. Si una falla, no se reescribe
# la bitácora: hay que corregir el diseño o la frase.

exigir <- function(condicion, mensaje) {
  if (!isTRUE(condicion)) stop("Falla el relato: ", mensaje, call. = FALSE)
  cat("  ok  ", mensaje, "\n", sep = "")
}

verificar_relato <- function(res) {
  cat("Verificación del relato\n")
  np <- res$np
  m <- np$modelos
  exigir(as.character(getRversion()) == "4.6.1", "R es 4.6.1")
  exigir(packageVersion("glmnet") == "5.0", "glmnet es 5.0")
  exigir(np$datos$semilla == 20260908L, "la semilla es 20260908")
  exigir(np$datos$n == 80L && np$datos$p == 20L, "n = 80 y p = 20")
  exigir(np$datos$n_train == 55L && length(np$datos$prueba) == 25L, "55 de entrenamiento y 25 de prueba")
  exigir(abs(np$datos$sigma - 1.2) < 1e-12, "el ruido de y tiene desviación 1.2")
  exigir(length(unique(np$datos$pliegues)) == 5L, "la validación cruzada usa cinco pliegues")

  exigir(np$cor12 > 0.95, "x01 y x02 están muy correlacionadas")
  exigir(np$vif > 20, "el VIF de x01 es alto")
  exigir(np$kappa > 30, "X'X centrada está mal condicionada")
  exigir(np$rango_aumentada == 21L, "con n > p la matriz aumentada tiene rango completo")

  b <- np$datos$beta
  exigir(all(b[c("x01", "x02", "x05", "x08")] == c(3, -2, 4, 1.5)), "las cuatro señales son las del diseño")
  exigir(b[["x03"]] == 0, "x03 es la trampa, con coeficiente verdadero 0")

  exigir(m$mco$rmse < m$lasso$rmse && m$lasso$rmse < m$ridge$rmse,
         "en n > p, MCO gana en RMSE y LASSO queda entre MCO y Ridge")
  exigir(m$ridge$rss > m$mco$rss && m$lasso$rss > m$mco$rss,
         "Ridge y LASSO aceptan un RSS de entrenamiento peor que MCO")
  exigir(m$ridge$l2 < m$mco$l2, "Ridge encoge la norma L2")
  exigir(m$lasso$l1 < m$mco$l1, "LASSO encoge la norma L1")

  exigir(all(sign(m$mco$coef[c("x01", "x02", "x05", "x08")]) == c(1, -1, 1, 1)),
         "MCO conserva el signo de las cuatro señales")

  exigir(m$lasso$coef[["x02"]] == 0, "LASSO en λmín apaga x02")
  exigir(m$lasso$coef[["x03"]] == 0, "LASSO en λmín apaga la trampa x03")
  exigir(m$lasso$coef[["x01"]] != 0 && m$lasso$coef[["x05"]] != 0 && m$lasso$coef[["x08"]] != 0,
         "LASSO en λmín conserva x01, x05 y x08")
  exigir(length(m$lasso$vivos) < 20L, "LASSO pone ceros exactos")
  exigir(length(m$ridge$vivos) == 20L && all(pendientes(m$ridge$coef) != 0),
         "Ridge no pone ningún cero exacto")
  exigir(m$ridge$coef[["x02"]] != 0 && m$ridge$coef[["x03"]] != 0,
         "Ridge no apaga ni a x02 ni a la trampa")
  exigir(abs(m$ridge$coef[["x01"]]) < abs(m$mco$coef[["x01"]]),
         "Ridge encoge x01 respecto de MCO")

  exigir(setequal(m$lasso_lambda1$vivos, c("x01", "x05", "x08")),
         "λ = 1 deja solo x01, x05 y x08")
  exigir(m$lasso_lambda1$coef[["x02"]] == 0, "con λ = 1, x02 ya está muerta")
  exigir(abs(m$lasso_lambda1$coef[["x01"]]) < abs(m$lasso$coef[["x01"]]),
         "λ = 1 encoge x01 más que λmín")
  exigir(abs(log(m$lasso$lambda) - log(1)) > log(2),
         "λmín de LASSO no está cerca de 1")

  exigir(setequal(m$lasso_1se$vivos, c("x01", "x05", "x08")),
         "λ1EE deja solo x01, x05 y x08")
  exigir(length(m$lasso_1se$vivos) < length(m$lasso$vivos),
         "λ1EE es un modelo más corto que λmín")
  exigir(m$lasso_1se$lambda > m$lasso$lambda, "λ1EE es mayor que λmín")

  falsos <- intersect(m$lasso$vivos, names(b)[b == 0])
  exigir(length(falsos) >= 1L, "λmín de LASSO deja al menos un falso positivo")
  exigir(max(abs(m$lasso$coef[falsos])) < min(abs(m$lasso$coef[c("x05", "x08")])),
         "los falsos positivos son más chicos que x05 y x08")

  esc <- res$escalas
  exigir(setequal(vivos(esc$crudo$coef), "x05"),
         "sin estandarizar solo sobrevive x05")
  exigir(esc$crudo$coef[["x08"]] == 0, "sin estandarizar x08 queda en cero")
  exigir(abs(esc$estandarizado$coef[["x08"]]) > 100,
         "con estandarizar x08 vuelve con un coeficiente enorme")
  exigir(
    abs(esc$estandarizado$coef[["x08"]] - m$lasso$coef[["x08"]] * 100) < 1e-6 &&
      abs(esc$estandarizado$coef[["x05"]] - m$lasso$coef[["x05"]] / 1000) < 1e-8,
    "estandarizar devuelve el mismo LASSO en las unidades nuevas"
  )

  pn <- res$pn
  exigir(pn$datos$n_train == 32L && pn$datos$p == 40L, "el segundo escenario tiene n = 32 y p = 40")
  exigir(pn$rango == 32L, "el rango de (1 | X) es 32")
  exigir(pn$n_parametros == 41L, "hay 41 parámetros, p + 1")
  exigir(pn$n_indefinidos == 9L, "lm deja 9 coeficientes sin definir")
  exigir(pn$rss_lm < 1e-8, "el lm de rango deficiente interpola el entrenamiento")
  exigir(length(pn$aviso_prediccion) >= 1L, "predict avisa que el ajuste es de rango deficiente")
  exigir(abs(pn$autovalor_min) < 1e-8, "X'X centrada es numéricamente singular")
  exigir(pn$autovalor_min_penalizado > 0.5, "sumar λI con λ = 1 vuelve invertible la matriz")
  exigir(pn$lasso$rmse < pn$rmse_pinv && pn$lasso$rmse < pn$ridge$rmse && pn$lasso$rmse < pn$rmse_lm,
         "cuando p > n, LASSO tiene menor RMSE que la seudoinversa, Ridge y lm")
  exigir(pn$ridge$rmse > pn$rmse_pinv, "Ridge no mejora a la seudoinversa")
  exigir(pn$lasso$coef[["x02"]] == 0 && pn$lasso$coef[["x01"]] != 0,
         "en p > n, LASSO otra vez se queda con x01 y apaga x02")

  libro <- res$ridge_libro
  fila0 <- libro[libro$lambda == 0, ]
  exigir(fila0$max_abs_dif_mco < 1e-8, "la fórmula cerrada con λ = 0 coincide con MCO")
  b_grande <- ridge_cerrado(np$partes$Xtr, np$partes$ytr, 1e8)
  exigir(max(abs(pendientes(b_grande))) < 1e-4, "con λ enorme la fórmula cerrada manda β a cero")

  invisible(TRUE)
}
