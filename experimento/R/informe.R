# Rellena la bitácora con los números de esta corrida y la pasa a HTML.

md_tabla <- function(df, digitos = NULL) {
  columnas <- names(df)
  encabezado <- paste0("| ", paste(columnas, collapse = " | "), " |")
  separador <- paste0("| ", paste(rep("---", length(columnas)), collapse = " | "), " |")
  filas <- apply(df, 1, function(fila) {
    paste0("| ", paste(fila, collapse = " | "), " |")
  })
  paste(c(encabezado, separador, filas), collapse = "\n")
}

redondear_columnas <- function(df, digitos) {
  out <- df
  for (nm in names(digitos)) {
    out[[nm]] <- fmt(out[[nm]], digitos[[nm]])
  }
  out
}

parrafo_sensibilidad <- function(tab) {
  base <- tab[tab$sd == SD_IDIO, ]
  otras <- tab[tab$sd != SD_IDIO, ]
  trozos <- apply(otras, 1, function(fila) {
    sd <- fmt(as.numeric(fila[["sd"]]), 2)
    cor <- fmt_cor(as.numeric(fila[["cor"]]))
    b02 <- as.numeric(fila[["b_x02"]])
    destino <- if (b02 == 0) {
      "LASSO siguió dejando x02 en cero"
    } else {
      sprintf("LASSO ya no apagó x02: el coeficiente fue %s", fmt(b02, 3))
    }
    rmse_l <- as.numeric(fila[["rmse_lasso"]])
    rmse_m <- as.numeric(fila[["rmse_mco"]])
    duelo <- if (rmse_l < rmse_m) {
      sprintf("en test LASSO quedó por debajo de MCO (%s frente a %s)", fmt(rmse_l, 3), fmt(rmse_m, 3))
    } else {
      sprintf("MCO siguió por debajo en test (%s frente a %s de LASSO)", fmt(rmse_m, 3), fmt(rmse_l, 3))
    }
    sprintf("Con desviación idiosincrática %s la correlación fue %s, %s, y %s.", sd, cor, destino, duelo)
  })
  paste(
    sprintf(
      "El diseño principal usa desviación idiosincrática %s (correlación %s, x02 de LASSO en %s).",
      fmt(base$sd, 2), fmt_cor(base$cor), fmt(base$b_x02, 3)
    ),
    paste(trozos, collapse = " "),
    "La misma semilla y las mismas extracciones: solo cambia cuánto ruido propio tiene el bloque. El diamante no mata a la gemela por decreto."
  )
}

parrafo_lm <- function(pn) {
  sprintf(
    "lm deja el RSS de entrenamiento en %s y %d coeficientes sin definir por singularidades. En prueba, ese ajuste marca RMSE %s. La seudoinversa, la solución de norma mínima entre las que minimizan el RSS, marca %s.",
    fmt(pn$rss_lm, 3), pn$n_indefinidos, fmt(pn$rmse_lm, 3), fmt(pn$rmse_pinv, 3)
  )
}

falsos_texto <- function(modelo, beta) {
  falsos <- intersect(modelo$vivos, names(beta)[beta == 0])
  if (!length(falsos)) return("(ninguno)")
  partes <- paste0(falsos, " = ", fmt(modelo$coef[falsos], 3))
  y_lista(partes)
}

textos_del_informe <- function(res) {
  m <- res$np$modelos
  b <- res$np$datos$beta
  falsos <- intersect(m$lasso$vivos, names(b)[b == 0])
  resto <- setdiff(names(b), c("x01", "x02", "x03", "x05", "x08"))
  ceros_resto <- sum(m$lasso$coef[resto] == 0)
  mensaje <- res$pn$aviso_prediccion[1]
  texto <- list(
    r_version = as.character(getRversion()),
    glmnet_version = as.character(packageVersion("glmnet")),
    ggplot2_version = as.character(packageVersion("ggplot2")),
    semilla = "20260908",
    sd_idio = fmt(SD_IDIO, 2),
    sigma = fmt(SIGMA, 1),
    cor_pob = fmt_cor(1 / (1 + SD_IDIO^2)),
    n = "80", p = "20", n_train = "55", n_test = "25",
    cor12 = fmt_cor(res$np$cor12),
    cor13 = fmt_cor(res$np$cor13),
    vif = fmt(res$np$vif, 2),
    kappa = fmt(res$np$kappa, 1),
    rmse_mco = fmt(m$mco$rmse, 3),
    rmse_ridge = fmt(m$ridge$rmse, 3),
    rmse_lasso = fmt(m$lasso$rmse, 3),
    lam_ridge = fmt(m$ridge$lambda, 3),
    lam_lasso = fmt(m$lasso$lambda, 3),
    lam_1se = fmt(m$lasso_1se$lambda, 3),
    nnz_lasso = as.character(length(m$lasso$vivos)),
    nnz_1se = as.character(length(m$lasso_1se$vivos)),
    nnz_lam1 = as.character(length(m$lasso_lambda1$vivos)),
    n_ceros_lasso = as.character(20L - length(m$lasso$vivos)),
    n_ceros_resto = as.character(ceros_resto),
    soporte_lasso_txt = y_lista(m$lasso$vivos),
    soporte_1se_txt = y_lista(m$lasso_1se$vivos),
    soporte_lam1_txt = y_lista(m$lasso_lambda1$vivos),
    falsos_txt = falsos_texto(m$lasso, b),
    falso_max = fmt(max(abs(m$lasso$coef[falsos])), 3),
    intercepto_lasso = fmt(m$lasso$intercepto, 3),
    rss_mco = fmt(m$mco$rss, 2),
    rss_ridge = fmt(m$ridge$rss, 2),
    rss_lasso = fmt(m$lasso$rss, 2),
    l1_mco = fmt(m$mco$l1, 2),
    l1_lasso = fmt(m$lasso$l1, 2),
    l2_mco = fmt(m$mco$l2, 2),
    l2_ridge = fmt(m$ridge$l2, 2),
    max_resto_mco = fmt(max(abs(m$mco$coef[resto])), 3),
    max_resto_ridge = fmt(max(abs(m$ridge$coef[resto])), 3),
    max_resto_lasso = fmt(max(abs(m$lasso$coef[resto])), 3),
    b_lam1_x01 = fmt(m$lasso_lambda1$coef[["x01"]], 3),
    b_lam1_x05 = fmt(m$lasso_lambda1$coef[["x05"]], 3),
    b_lam1_x08 = fmt(m$lasso_lambda1$coef[["x08"]], 3),
    esc_lam_std = fmt(res$escalas$estandarizado$lambda, 3),
    esc_lam_raw = fmt(res$escalas$crudo$lambda, 3),
    esc_x05_std = fmt(res$escalas$estandarizado$coef[["x05"]], 5),
    esc_x08_std = fmt(res$escalas$estandarizado$coef[["x08"]], 3),
    esc_x05_raw = fmt(res$escalas$crudo$coef[["x05"]], 5),
    esc_x08_raw = fmt(res$escalas$crudo$coef[["x08"]], 3),
    esc_nnz_std = as.character(length(vivos(res$escalas$estandarizado$coef))),
    esc_nnz_raw = as.character(length(vivos(res$escalas$crudo$coef))),
    pn_n = as.character(res$pn$datos$n),
    pn_p = "40",
    pn_n_train = "32",
    pn_n_test = as.character(length(res$pn$datos$prueba)),
    pn_rango = as.character(res$pn$rango),
    pn_n_param = as.character(res$pn$n_parametros),
    pn_indefinidos = as.character(res$pn$n_indefinidos),
    pn_mensaje_lm = mensaje,
    pn_rmse_mco = fmt(res$pn$rmse_pinv, 3),
    pn_rmse_lm = fmt(res$pn$rmse_lm, 3),
    pn_rmse_ridge = fmt(res$pn$ridge$rmse, 3),
    pn_rmse_lasso = fmt(res$pn$lasso$rmse, 3),
    pn_nnz_lasso = as.character(length(res$pn$lasso$vivos)),
    pn_lam_ridge = fmt(res$pn$ridge$lambda, 3),
    pn_lam_lasso = fmt(res$pn$lasso$lambda, 3),
    pn_soporte_txt = y_lista(res$pn$lasso$vivos),
    pn_b_x01 = fmt(res$pn$lasso$coef[["x01"]], 3),
    pn_b_x02 = fmt(res$pn$lasso$coef[["x02"]], 3),
    pn_b_x05 = fmt(res$pn$lasso$coef[["x05"]], 3),
    pn_b_x08 = fmt(res$pn$lasso$coef[["x08"]], 3),
    pn_min_eigen = fmt_sci(res$pn$autovalor_min, 1),
    pn_min_eigen_pen = fmt(res$pn$autovalor_min_penalizado, 3),
    parrafo_sensibilidad = parrafo_sensibilidad(res$sensibilidad),
    parrafo_lm = parrafo_lm(res$pn),
    tabla_diseno = md_tabla(data.frame(
      Pieza = c(
        "Tamaño, n > p",
        "Bloque colineal",
        "Señales",
        "Trampa",
        "Ruido de y",
        "Validación",
        "Tamaño, p > n"
      ),
      `Decisión` = c(
        "n = 80, p = 20, entrenamiento = 55, prueba = 25",
        "x01, x02 y x03 = z + 0.16 e; correlación poblacional ≈ .975",
        "x01 = 3, x02 = −2, x05 = 4, x08 = 1.5",
        "x03 colineal y coeficiente verdadero 0",
        "desviación 1.2",
        "cinco pliegues, los mismos para Ridge y LASSO",
        "n = 57, p = 40, entrenamiento = 32, prueba = 25"
      ),
      check.names = FALSE
    )),
    tabla_coef = local({
      tab <- redondear_columnas(
        tabla_coeficientes_foco(res$np),
        c(MCO = 3, Ridge = 3, LASSO = 3)
      )
      tab$Verdadero <- c("3", "-2", "0", "4", "1.5")
      names(tab)[1] <- "Variable"
      md_tabla(tab)
    }),
    tabla_rmse = md_tabla(data.frame(
      Modelo = c("MCO", "Ridge", "LASSO"),
      RMSE = fmt(c(m$mco$rmse, m$ridge$rmse, m$lasso$rmse), 3),
      lambda = c("—", fmt(m$ridge$lambda, 3), fmt(m$lasso$lambda, 3)),
      `No ceros` = c("20", "20", as.character(length(m$lasso$vivos))),
      check.names = FALSE
    )),
    tabla_escalas = md_tabla(data.frame(
      Ajuste = c("Sin estandarizar", "Con estandarizar"),
      `lambda min` = fmt(c(res$escalas$crudo$lambda, res$escalas$estandarizado$lambda), 3),
      x05 = fmt(c(res$escalas$crudo$coef[["x05"]], res$escalas$estandarizado$coef[["x05"]]), 5),
      x08 = fmt(c(res$escalas$crudo$coef[["x08"]], res$escalas$estandarizado$coef[["x08"]]), 3),
      `No ceros` = c(
        as.character(length(vivos(res$escalas$crudo$coef))),
        as.character(length(vivos(res$escalas$estandarizado$coef)))
      ),
      check.names = FALSE
    )),
    tabla_pn = md_tabla(data.frame(
      Modelo = c("Seudoinversa", "lm", "Ridge", "LASSO"),
      RMSE = fmt(c(res$pn$rmse_pinv, res$pn$rmse_lm, res$pn$ridge$rmse, res$pn$lasso$rmse), 3),
      lambda = c("—", "—", fmt(res$pn$ridge$lambda, 3), fmt(res$pn$lasso$lambda, 3)),
      `No ceros` = c(
        "—",
        "—",
        as.character(length(res$pn$ridge$vivos)),
        as.character(length(res$pn$lasso$vivos))
      ),
      check.names = FALSE
    )),
    tabla_sensibilidad = md_tabla(redondear_columnas(
      data.frame(
        sd = res$sensibilidad$sd,
        cor = res$sensibilidad$cor,
        VIF = res$sensibilidad$vif,
        `RMSE MCO` = res$sensibilidad$rmse_mco,
        `RMSE LASSO` = res$sensibilidad$rmse_lasso,
        `LASSO x02` = res$sensibilidad$b_x02,
        check.names = FALSE
      ),
      c(sd = 2, cor = 3, VIF = 2, `RMSE MCO` = 3, `RMSE LASSO` = 3, `LASSO x02` = 3)
    )),
    tabla_ridge_cerrado = md_tabla(redondear_columnas(
      res$ridge_libro[, c("lambda", "x01", "x02", "x03", "x05", "x08")],
      c(lambda = 0, x01 = 3, x02 = 3, x03 = 3, x05 = 3, x08 = 3)
    ))
  )
  foco <- c("x01", "x02", "x03", "x05", "x08")
  for (met in c("mco", "ridge", "lasso")) {
    for (v in foco) {
      texto[[paste0("b_", met, "_", v)]] <- fmt(m[[met]]$coef[[v]], 3)
    }
  }
  texto
}

escribir_informe <- function(res, raiz) {
  plantilla <- file.path(raiz, "informe", "bitacora.md.tmpl")
  salida_md <- file.path(raiz, "informe", "bitacora.md")
  lineas <- readLines(plantilla, encoding = "UTF-8", warn = FALSE)
  txt <- paste(lineas, collapse = "\n")
  vals <- textos_del_informe(res)
  for (nm in names(vals)) {
    txt <- gsub(paste0("{{", nm, "}}"), vals[[nm]], txt, fixed = TRUE)
  }
  faltan <- unique(regmatches(txt, gregexpr("\\{\\{[A-Za-z0-9_]+\\}\\}", txt))[[1]])
  if (length(faltan)) {
    stop("Quedaron placeholders sin rellenar: ", paste(faltan, collapse = ", "))
  }
  writeLines(txt, salida_md, useBytes = TRUE)
  anterior <- getwd()
  setwd(file.path(raiz, "informe"))
  status <- system2(
    "pandoc",
    c(
      "bitacora.md",
      "-o", "bitacora.html",
      "--standalone",
      "--embed-resources",
      "--css", "estilo.css"
    )
  )
  setwd(anterior)
  if (!identical(status, 0L)) stop("pandoc no pudo escribir el HTML")
  exportar_pdf(raiz)
  invisible(salida_md)
}

exportar_pdf <- function(raiz) {
  html <- normalizePath(file.path(raiz, "informe", "bitacora.html"))
  pdf <- file.path(raiz, "informe", "Informe_Ridge_LASSO.pdf")
  chrome <- Sys.which("google-chrome")
  if (!nzchar(chrome)) chrome <- Sys.which("google-chrome-stable")
  if (!nzchar(chrome)) stop("No está google-chrome; no pude escribir el PDF")
  perfil <- tempfile("chrome-pdf-")
  dir.create(perfil)
  # Chrome a veces escribe el PDF y no termina. timeout lo cierra.
  status <- system2(
    "timeout",
    c(
      "--signal=KILL", "25",
      chrome,
      "--headless", "--disable-gpu", "--no-sandbox",
      "--disable-dev-shm-usage", "--disable-extensions",
      "--no-pdf-header-footer",
      paste0("--user-data-dir=", perfil),
      paste0("--print-to-pdf=", pdf),
      paste0("file://", html)
    )
  )
  if (!file.exists(pdf) || file.info(pdf)$size < 50000) {
    stop("Chrome no pudo escribir el PDF")
  }
  if (!status %in% c(0L, 124L, 137L)) {
    warning("Chrome terminó con código ", status, "; el PDF sí quedó escrito")
  }
  invisible(pdf)
}
