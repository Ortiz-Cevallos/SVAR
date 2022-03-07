---
title: "VECTORES AUTOREGRESIVOS (VAR)"
author: "[Luis Ortiz-Cevallos](https://ortiz-cevallos.github.io/MYSELF/)"
date: "2022-03-07"
site: bookdown::bookdown_site
output: bookdown::gitbook
code_download: true
code_folding: hide
documentclass: book
bibliography: [book.bib, packages.bib]
biblio-style: apalike
link-citations: yes
---





# INTRODUCCIÓN

En esta nota para la obtención de un Vector Autorregresivo Estructural (SVAR) se utilizará el paquete en R *SVARS* elaborado por @SVAR21. A continuación se muestra la carga de este paquete junto a otros que serán de utilidad cuando se analicen los resultados. 


```r
#install.packages("svars")
library("svars")
library("ggplot2")
```

Enseguida se cargará una base de datos en frecuencia trimestral que contiene la brecha producto, la inflación de trimestre a trimestre del deflactor del PIB y una tasa de interés nominal del fondo federal.

```r
data("USA")
usa<-as.zoo(USA)
autoplot(usa, facets = T) + theme_bw() + ylab('')
```

![](01-conceptos_files/figure-epub3/unnamed-chunk-2-1.png)<!-- -->

A continuación se replicará los resultados obtenidos por @Herwartz2016, específicamente se estimaran los shocks estructurales a través de la metodología de cambios en volatilidades. Ello se logra a través de la función *id.cv()* dentro del cual hay que indicarle la fecha a partir de la cual se dio ese cambio estructural.

@Herwartz2016 argumentan que el punto de quiebra se dio en el tercer trimestre de 1979 donde aconteció un cambio en la política de la Reserva Federal que redujo la volatilidad de las variables macroeconómicas (ver @Stock2003)

Lo anterior es posible visualizarlo.
Un primer paso es visualizar cada una de las series.

```r
autoplot(usa, facets = T) + theme_bw() + ylab('Evolución de series de USA')
```

![](01-conceptos_files/figure-epub3/unnamed-chunk-3-1.png)<!-- -->

Seguidamente se estimara la forma reducida de un VAR, aplicaremos una especificación que incluya intercepto y hasta seis rezagos (p=6)

Un primer paso es visualizar cada una de las series.


```r
plain.var <- vars::VAR(USA, p = 6, type = 'const')
```



# EJERCICIO

## BLA Bla

La estructura del modelo determina el orden de prelación de cada una de las variables y las restricciones de corto plazo.
El vector está dado por X=[y, $\pi$, i, NPL, L, e] y se establece un SVAR de la forma $AX=X_{-1}+u$, en el que A contiene las restricciones y con el cual se identifica los shocks.
Los datos utilizados tienen la siguiente características:}

\begin{eqnarray}
AX&=\begin{bmatrix}
1 & 0 & 1 & 0 & 0 & 0  \\
1 & 1 & 0 & 0 & 0 & 0 \\
0 & 1 & 1 & 1 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 \\
1 & 1 & 1 & 1 & 1 & 0 \\
1 & 1 & 1 & 1 & 1 & 1 \end{bmatrix}
\begin{bmatrix}
y \\
\pi \\
i\\
NPL\\
L\\
e\end{bmatrix}
\end{eqnarray}



