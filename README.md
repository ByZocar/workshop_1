# Ridge y LASSO: un experimento reproducible en R

Bitácora de Estadística y Probabilidad 2. El repositorio reemplaza el taller de ETL que estaba aquí. Ahora contiene el experimento sintético de regularización y el informe cuyos números salen de ese código.

El diseño conserva el ejemplo del trabajo: cuatro señales (`x01 = 3`, `x02 = -2`, `x05 = 4`, `x08 = 1.5`), una trampa colineal (`x03` con coeficiente verdadero 0), ruido de desviación 1.2, semilla `20260908`, glmnet y validación cruzada de cinco pliegues. Un segundo escenario usa más columnas que filas de entrenamiento. Las cifras del informe son las de esta corrida, no una tabla copiada a mano.

## Qué sostiene el código

- Con n de entrenamiento = 55 y p = 20, MCO gana el RMSE de prueba. Ridge encoge el bloque colineal y no pone ceros exactos. LASSO apaga `x02`, que era señal, y conserva `x01`.
- λ = 1 no es el valle de la validación cruzada. λmín y la regla de un error estándar dejan modelos de distinto tamaño.
- Si `x05` se mide en unidades mil veces más grandes y `x08` en unidades cien veces más chicas, LASSO sin estandarizar se queda solo con `x05`. Con la estandarización de glmnet reaparece `x08`.
- Con n de entrenamiento = 32 y p = 40, la matriz `(1 | X)` tiene rango 32. `lm` deja coeficientes sin definir. LASSO predice mejor que la seudoinversa y que Ridge.

`experimento/R/verificar.R` exige esas afirmaciones. Si una corrida deja de cumplirlas, el script se detiene antes de reescribir el informe.

## Cómo correrlo

Requisitos: R 4.6.1, glmnet 5.0 y ggplot2. Pandoc hace falta solo para el HTML.

```bash
Rscript experimento/run_all.R
```

El comando regenera:

- `experimento/salida/figuras/` (figuras 1 a 9)
- `experimento/salida/tablas/`
- `experimento/salida/resultados.rds`
- `informe/bitacora.md` y `informe/bitacora.html`

El informe para entregar es [`informe/Informe_Ridge_LASSO.pdf`](informe/Informe_Ridge_LASSO.pdf). La misma bitácora está en [`informe/bitacora.md`](informe/bitacora.md).

## Estructura

| Ruta | Contenido |
| --- | --- |
| `experimento/run_all.R` | Orquesta la corrida |
| `experimento/R/diseno.R` | Proceso generador y semilla |
| `experimento/R/estimar.R` | MCO, glmnet, escalas, p > n, fórmula cerrada |
| `experimento/R/figuras.R` | Figuras 1 a 9 |
| `experimento/R/verificar.R` | Contraste del relato |
| `experimento/R/informe.R` | Rellena la bitácora |
| `informe/bitacora.md.tmpl` | Texto del informe, con los huecos numéricos |
