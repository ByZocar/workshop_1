# Corre el experimento completo y regenera figuras, tablas e informe.
# Desde la raíz del repositorio:
#   Rscript experimento/run_all.R

ruta_de_este_archivo <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  marca <- grep("^--file=", args, value = TRUE)
  if (length(marca)) {
    return(normalizePath(sub("^--file=", "", marca[[1]])))
  }
  for (i in rev(seq_len(sys.nframe()))) {
    candidato <- sys.frame(i)$ofile
    if (!is.null(candidato)) return(normalizePath(candidato))
  }
  normalizePath("experimento/run_all.R")
}

raiz <- normalizePath(file.path(dirname(ruta_de_este_archivo()), ".."))
setwd(raiz)

suppressPackageStartupMessages({
  library(glmnet)
  library(ggplot2)
})
options(warn = 1)

sys.source(file.path(raiz, "experimento", "R", "util.R"), envir = globalenv())
sys.source(file.path(raiz, "experimento", "R", "diseno.R"), envir = globalenv())
sys.source(file.path(raiz, "experimento", "R", "estimar.R"), envir = globalenv())
sys.source(file.path(raiz, "experimento", "R", "figuras.R"), envir = globalenv())
sys.source(file.path(raiz, "experimento", "R", "verificar.R"), envir = globalenv())
sys.source(file.path(raiz, "experimento", "R", "informe.R"), envir = globalenv())

cat("Escenario n > p\n")
np <- ajustar_escenario(fabricar(80L, 20L, 55L))
cat("Escalas\n")
escalas <- ajustar_escalas(np)
cat("Escenario p > n\n")
pn <- ajustar_pn(fabricar(57L, 40L, 32L))
cat("Sensibilidad y fórmula cerrada\n")
sensibilidad <- sensibilidad_sd()
ridge_libro <- ridge_de_libro(np)

res <- list(
  np = np,
  escalas = escalas,
  pn = pn,
  sensibilidad = sensibilidad,
  ridge_libro = ridge_libro
)

dir_salida <- file.path(raiz, "experimento", "salida")
dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(dir_salida, "resultados.rds"))

cat("Figuras\n")
hacer_figuras(res, file.path(dir_salida, "figuras"))

cat("Tablas\n")
dir_tab <- file.path(dir_salida, "tablas")
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)
write.csv(tabla_coeficientes_foco(np), file.path(dir_tab, "coeficientes_foco.csv"), row.names = FALSE)
write.csv(sensibilidad, file.path(dir_tab, "sensibilidad_sd.csv"), row.names = FALSE)
write.csv(ridge_libro, file.path(dir_tab, "ridge_cerrado.csv"), row.names = FALSE)

verificar_relato(res)
escribir_informe(res, raiz)
cat("Listo. Informe en informe/bitacora.md\n")
