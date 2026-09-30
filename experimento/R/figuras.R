# Figuras 1 a 9. La figura 1 es un esquema: no usa los datos del laboratorio.

tema_informe <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey25"),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank()
    )
}

col_foco <- c(
  x01 = "#1B4F72",
  x02 = "#922B21",
  x03 = "#B9770E",
  x05 = "#1D8348",
  x08 = "#6C3483"
)

col_metodo <- c(
  Verdadero = "#1C2833",
  MCO = "#5D6D7E",
  Ridge = "#2471A3",
  LASSO = "#CA6F1E"
)

guardar <- function(plot, ruta, ancho = 8, alto = 5.2) {
  ggplot2::ggsave(ruta, plot, width = ancho, height = alto, dpi = 160, bg = "white")
}

# Descenso por coordenadas para el esquema 2D.
# f(b) = 0.5 (b - mu)' A (b - mu) + λ (|b1| + |b2|).
lasso_2d <- function(A, mu, lambda, iter = 300L) {
  b <- c(0, 0)
  for (k in seq_len(iter)) {
    for (j in 1:2) {
      otro <- 3L - j
      r <- A[j, otro] * b[otro] - sum(A[j, ] * mu)
      a <- A[j, j]
      candidato <- -r / a
      b[j] <- sign(candidato) * max(abs(candidato) - lambda / a, 0)
    }
  }
  b
}

ridge_2d <- function(A, mu, lambda) {
  drop(solve(A + lambda * diag(2), A %*% mu))
}

rss_2d <- function(b, A, mu) drop(t(b - mu) %*% A %*% (b - mu))

comprobar_lasso_2d <- function() {
  # Con A = I la solución es el umbral suave de mu.
  b <- lasso_2d(diag(2), c(2, 0.3), 0.5)
  if (max(abs(b - c(1.5, 0))) > 1e-8) {
    stop("El esquema LASSO 2D no recuperó el umbral suave")
  }
}

elegir_contacto <- function() {
  mu <- c(2.35, 0.85)
  A <- matrix(c(1.25, 0.32, 0.32, 0.95), 2, 2)
  lams <- seq(0.02, 4, length.out = 250)
  ridge <- lapply(lams, function(lam) {
    b <- ridge_2d(A, mu, lam)
    c(b, rss = rss_2d(b, A, mu), radio = sqrt(sum(b^2)))
  })
  lasso <- lapply(lams, function(lam) {
    b <- lasso_2d(A, mu, lam)
    c(b, rss = rss_2d(b, A, mu), l1 = sum(abs(b)))
  })
  R <- do.call(rbind, ridge)
  L <- do.call(rbind, lasso)
  ok_r <- R[, 1] > 0.25 & R[, 2] > 0.25
  ok_l <- (abs(L[, 1]) < 1e-6 & L[, 2] > 0.35) | (abs(L[, 2]) < 1e-6 & L[, 1] > 0.35)
  if (!any(ok_r) || !any(ok_l)) stop("No hubo contacto geométrico usable")
  pares <- which(outer(ok_r, ok_l, FUN = `&`), arr.ind = TRUE)
  dif <- abs(R[pares[, 1], "rss"] - L[pares[, 2], "rss"])
  mejor <- pares[which.min(dif), ]
  if (min(dif) > 0.04) stop("No coincidió el nivel de RSS del círculo y del diamante")
  list(
    mu = mu, A = A,
    b_ridge = R[mejor[1], 1:2],
    b_lasso = L[mejor[2], 1:2],
    rss = mean(c(R[mejor[1], "rss"], L[mejor[2], "rss"]))
  )
}

figura_geometria <- function(ruta) {
  comprobar_lasso_2d()
  contacto <- elegir_contacto()
  mu <- contacto$mu
  A <- contacto$A
  des <- eigen(A, symmetric = TRUE)
  radio_elipse <- sqrt(contacto$rss / des$values)
  base <- des$vectors %*% diag(radio_elipse, 2)
  ang <- seq(0, 2 * pi, length.out = 400)
  elipse <- data.frame(
    b1 = mu[1] + base[1, 1] * cos(ang) + base[1, 2] * sin(ang),
    b2 = mu[2] + base[2, 1] * cos(ang) + base[2, 2] * sin(ang)
  )
  radio <- sqrt(sum(contacto$b_ridge^2))
  t_l1 <- sum(abs(contacto$b_lasso))
  circulo <- data.frame(b1 = radio * cos(ang), b2 = radio * sin(ang), restriccion = "Círculo L2")
  diamante <- data.frame(
    b1 = c(t_l1, 0, -t_l1, 0, t_l1),
    b2 = c(0, t_l1, 0, -t_l1, 0),
    restriccion = "Diamante L1"
  )
  puntos <- data.frame(
    b1 = c(mu[1], contacto$b_ridge[1], contacto$b_lasso[1]),
    b2 = c(mu[2], contacto$b_ridge[2], contacto$b_lasso[2]),
    punto = c("Mínimo sin restricción", "Contacto Ridge", "Contacto LASSO")
  )
  margen <- 0.35
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
    ggplot2::geom_path(
      data = elipse,
      ggplot2::aes(b1, b2),
      colour = "grey35",
      linetype = "dashed",
      linewidth = 0.7
    ) +
    ggplot2::geom_path(
      data = circulo,
      ggplot2::aes(b1, b2, colour = restriccion),
      linewidth = 1
    ) +
    ggplot2::geom_path(
      data = diamante,
      ggplot2::aes(b1, b2, colour = restriccion),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = puntos,
      ggplot2::aes(b1, b2, shape = punto),
      size = 2.6
    ) +
    ggplot2::scale_colour_manual(values = c("Círculo L2" = "#2471A3", "Diamante L1" = "#CA6F1E")) +
    ggplot2::coord_equal(
      xlim = range(elipse$b1, circulo$b1, diamante$b1) + c(-margen, margen),
      ylim = range(elipse$b2, circulo$b2, diamante$b2) + c(-margen, margen)
    ) +
    ggplot2::labs(
      title = "Misma elipse de RSS, dos restricciones",
      subtitle = "Esquema. No usa datos del laboratorio. El círculo se toca fuera de los ejes; el diamante, en una esquina.",
      x = expression(beta[1]),
      y = expression(beta[2])
    ) +
    tema_informe()
  guardar(p, ruta, ancho = 8, alto = 5.6)
}

figura_bloque <- function(escenario, ruta) {
  nodos <- data.frame(
    x = c(0, -1.35, 0, 1.35),
    y = c(1.15, 0, 0, 0),
    etiqueta = c(
      "z\nlatente no observada",
      "x01\nseñal,  \u03b2 = 3",
      "x02\nseñal,  \u03b2 = \u22122",
      "x03\ntrampa,  \u03b2 = 0"
    )
  )
  flechas <- data.frame(
    x = 0, y = 0.92,
    xend = c(-1.35, 0, 1.35),
    yend = 0.28
  )
  subtitulo <- sprintf(
    "Entrenamiento: r(x01, x02) = %s;   r(x01, x03) = %s;   VIF(x01) = %s",
    fmt_cor(escenario$cor12), fmt_cor(escenario$cor13), fmt(escenario$vif, 2)
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = flechas,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.18, "cm"), type = "closed"),
      linewidth = 0.6,
      colour = "#1C2833"
    ) +
    ggplot2::geom_label(
      data = nodos,
      ggplot2::aes(x, y, label = etiqueta),
      linewidth = 0.4,
      size = 3.6,
      fill = "white",
      label.padding = ggplot2::unit(0.35, "lines")
    ) +
    ggplot2::coord_cartesian(xlim = c(-2.1, 2.1), ylim = c(-0.45, 1.7)) +
    ggplot2::labs(
      title = "Bloque colineal fabricado",
      subtitle = subtitulo
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey20"),
      plot.margin = ggplot2::margin(12, 12, 12, 12)
    )
  guardar(p, ruta, ancho = 8, alto = 4.6)
}

tabla_coeficientes_foco <- function(escenario) {
  vars <- c("x01", "x02", "x03", "x05", "x08")
  data.frame(
    variable = vars,
    Verdadero = escenario$datos$beta[vars],
    MCO = escenario$modelos$mco$coef[vars],
    Ridge = escenario$modelos$ridge$coef[vars],
    LASSO = escenario$modelos$lasso$coef[vars],
    row.names = NULL
  )
}

figura_coeficientes <- function(escenario, ruta) {
  tab <- tabla_coeficientes_foco(escenario)
  largo <- data.frame(
    variable = rep(tab$variable, 4),
    metodo = rep(c("Verdadero", "MCO", "Ridge", "LASSO"), each = nrow(tab)),
    valor = c(tab$Verdadero, tab$MCO, tab$Ridge, tab$LASSO)
  )
  largo$metodo <- factor(largo$metodo, levels = c("Verdadero", "MCO", "Ridge", "LASSO"))
  p <- ggplot2::ggplot(largo, ggplot2::aes(variable, valor, fill = metodo)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78), width = 0.72) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = col_metodo) +
    ggplot2::labs(
      title = "Coeficientes en las variables de interés",
      subtitle = "Ridge y LASSO en \u03bbmín de validación cruzada de cinco pliegues",
      x = NULL,
      y = "Coeficiente"
    ) +
    tema_informe()
  guardar(p, ruta)
}

camino_largo <- function(ajuste_cv) {
  beta <- as.matrix(coef(ajuste_cv$glmnet.fit))
  beta <- beta[rownames(beta) != "(Intercept)", , drop = FALSE]
  data.frame(
    variable = rep(rownames(beta), ncol(beta)),
    log_lambda = rep(log(ajuste_cv$lambda), each = nrow(beta)),
    lambda = rep(ajuste_cv$lambda, each = nrow(beta)),
    valor = as.numeric(beta),
    stringsAsFactors = FALSE
  )
}

figura_camino <- function(ajuste_cv, lambda_min, titulo, subtitulo, ruta) {
  d <- camino_largo(ajuste_cv)
  d$foco <- d$variable %in% names(col_foco)
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = d[!d$foco, ],
      ggplot2::aes(log_lambda, valor, group = variable),
      colour = "grey80",
      linewidth = 0.4
    ) +
    ggplot2::geom_line(
      data = d[d$foco, ],
      ggplot2::aes(log_lambda, valor, colour = variable, group = variable),
      linewidth = 0.9
    ) +
    ggplot2::geom_vline(xintercept = log(lambda_min), linetype = "dashed", linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = col_foco) +
    ggplot2::labs(
      title = titulo,
      subtitle = subtitulo,
      x = expression(log(lambda)),
      y = "Coeficiente"
    ) +
    tema_informe()
  guardar(p, ruta)
}

figura_validacion <- function(escenario, ruta) {
  curva <- function(cv, modelo) {
    data.frame(
      log_lambda = log(cv$lambda),
      cvm = cv$cvm,
      cvsd = cv$cvsd,
      modelo = modelo
    )
  }
  d <- rbind(
    curva(escenario$cv_ridge, "Ridge"),
    curva(escenario$cv_lasso, "LASSO")
  )
  guias <- data.frame(
    log_lambda = log(c(
      escenario$cv_lasso$lambda.min,
      escenario$cv_ridge$lambda.min,
      1
    )),
    etiqueta = c("\u03bbm\u00edn LASSO", "\u03bbm\u00edn Ridge", "\u03bb = 1"),
    colour = c("#CA6F1E", "#2471A3", "#1C2833")
  )
  p <- ggplot2::ggplot(d, ggplot2::aes(log_lambda, cvm, colour = modelo, fill = modelo)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = cvm - cvsd, ymax = cvm + cvsd),
      alpha = 0.15, colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_vline(
      data = guias,
      ggplot2::aes(xintercept = log_lambda),
      colour = guias$colour,
      linetype = "dashed",
      linewidth = 0.45,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = guias,
      ggplot2::aes(x = log_lambda, y = 11.4, label = etiqueta),
      colour = guias$colour,
      inherit.aes = FALSE,
      angle = 90,
      vjust = -0.4,
      hjust = 1,
      size = 3.1
    ) +
    ggplot2::scale_colour_manual(values = c(Ridge = "#2471A3", LASSO = "#CA6F1E")) +
    ggplot2::scale_fill_manual(values = c(Ridge = "#2471A3", LASSO = "#CA6F1E")) +
    ggplot2::coord_cartesian(ylim = c(0, 12)) +
    ggplot2::labs(
      title = "Error de validación cruzada",
      subtitle = "La banda es un error estándar. \u03bb = 1 queda a la derecha del valle.",
      x = expression(log(lambda)),
      y = "ECM de validación"
    ) +
    tema_informe()
  guardar(p, ruta, alto = 5.6)
}

figura_rmse <- function(nombres, valores, titulo, subtitulo, ruta) {
  d <- data.frame(
    modelo = factor(nombres, levels = nombres),
    rmse = valores
  )
  p <- ggplot2::ggplot(d, ggplot2::aes(modelo, rmse, fill = modelo)) +
    ggplot2::geom_col(width = 0.72, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = fmt(rmse, 3)), vjust = -0.4, size = 3.6) +
    ggplot2::scale_fill_manual(values = c(
      "MCO" = "#5D6D7E",
      "Seudoinversa" = "#5D6D7E",
      "lm" = "#85929E",
      "Ridge" = "#2471A3",
      "LASSO" = "#CA6F1E"
    )) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(title = titulo, subtitle = subtitulo, x = NULL, y = "RMSE de prueba") +
    tema_informe()
  guardar(p, ruta, alto = 4.8)
}

figura_escalas <- function(escalas, ruta) {
  armar <- function(bloque, etiqueta) {
    b <- pendientes(bloque$coef)
    data.frame(
      variable = c("x05", "x08"),
      valor = as.numeric(b[c("x05", "x08")]),
      ajuste = etiqueta,
      stringsAsFactors = FALSE
    )
  }
  d <- rbind(
    armar(escalas$crudo, "Sin estandarizar"),
    armar(escalas$estandarizado, "Con estandarizar")
  )
  d$ajuste <- factor(d$ajuste, levels = c("Sin estandarizar", "Con estandarizar"))
  d$etiqueta <- ifelse(d$valor == 0, "0 exacto", fmt(abs(d$valor), 4))
  techos <- tapply(abs(d$valor), d$ajuste, max)
  d$y_etiqueta <- abs(d$valor)
  for (aj in names(techos)) {
    sel <- d$ajuste == aj & abs(d$valor) < 0.05 * techos[[aj]]
    d$y_etiqueta[sel] <- 0.14 * techos[[aj]]
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(variable, abs(valor), fill = variable)) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(y = y_etiqueta, label = etiqueta),
      vjust = -0.2,
      size = 3.5
    ) +
    ggplot2::facet_wrap(~ajuste, scales = "free_y") +
    ggplot2::scale_fill_manual(values = col_foco) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(
      title = "LASSO después de cambiar las unidades",
      subtitle = "x05 se multiplicó por 1000 y x08 se dividió por 100. Cada panel usa su \u03bbmín. El eje es |\u03b2| en la escala medida.",
      x = NULL,
      y = expression(group("|", beta, "|"))
    ) +
    tema_informe()
  guardar(p, ruta)
}

hacer_figuras <- function(res, dir_fig) {
  dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
  figura_geometria(file.path(dir_fig, "figura_01_geometria.png"))
  figura_bloque(res$np, file.path(dir_fig, "figura_02_bloque.png"))
  figura_coeficientes(res$np, file.path(dir_fig, "figura_03_coeficientes.png"))
  figura_camino(
    res$np$cv_ridge, res$np$cv_ridge$lambda.min,
    "Camino de coeficientes de Ridge",
    "El bloque se encoge junto. En λmín ninguna señal está en cero exacto. La línea vertical es λmín.",
    file.path(dir_fig, "figura_04_camino_ridge.png")
  )
  figura_camino(
    res$np$cv_lasso, res$np$cv_lasso$lambda.min,
    "Camino de coeficientes de LASSO",
    "x02, señal gemela de x01, llega a cero exacto. La línea vertical es \u03bbmín.",
    file.path(dir_fig, "figura_05_camino_lasso.png")
  )
  figura_validacion(res$np, file.path(dir_fig, "figura_06_validacion.png"))
  figura_rmse(
    c("MCO", "Ridge", "LASSO"),
    c(res$np$modelos$mco$rmse, res$np$modelos$ridge$rmse, res$np$modelos$lasso$rmse),
    "RMSE de prueba cuando n de entrenamiento supera a p",
    "n = 80, p = 20, entrenamiento = 55, prueba = 25",
    file.path(dir_fig, "figura_07_rmse_np.png")
  )
  figura_escalas(res$escalas, file.path(dir_fig, "figura_08_escalas.png"))
  figura_rmse(
    c("Seudoinversa", "lm", "Ridge", "LASSO"),
    c(res$pn$rmse_pinv, res$pn$rmse_lm, res$pn$ridge$rmse, res$pn$lasso$rmse),
    "RMSE de prueba cuando p supera a n",
    "Entrenamiento = 32, p = 40, rango de (1 | X) = 32",
    file.path(dir_fig, "figura_09_rmse_pn.png")
  )
}
