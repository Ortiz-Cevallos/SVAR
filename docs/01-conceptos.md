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

## Paquete en R

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
```

A continuación se replicará los resultados obtenidos por @Herwartz2016, específicamente se estimaran los shocks estructurales a través de la metodología de cambios en volatilidades. Ello se logra a través de la función *id.cv()* dentro del cual hay que indicarle la fecha a partir de la cual se dio ese cambio estructural.

@Herwartz2016 argumentan que el punto de quiebra se dio en el tercer trimestre de 1979 donde aconteció un cambio en la política de la Reserva Federal que redujo la volatilidad de las variables macroeconómicas (ver @Stock2003)

Lo anterior es posible visualizarlo.
Un primer paso es visualizar cada una de las series.

```r
autoplot(usa, facets = T) + theme_bw() + ylab('Evolución de series de USA')
```

![](01-conceptos_files/figure-epub3/unnamed-chunk-3-1.png)<!-- -->

Seguidamente se estimara la forma reducida de un VAR, aplicaremos una especificación que incluya intercepto y hasta seis rezagos (p=6)


```r
plain.var <- vars::VAR(USA, p = 6, type = 'const')
```
Con base en el objeto VAR podemos estimar su forma estructural con la función *id.cv()* la cual introduce el *cambio en la varianza* que se observa en las series, sin embargo la aplicación de esta función requiere específicar el argumento *SB* en formato *ts*:

```r
usa.cv <- id.cv(plain.var, SB = c(1979, 3))
summary(usa.cv)
```

```
## 
## Identification Results
## ---------------------- 
## 
## Method: Changes in Volatility
## Sample size: 169
## Log-Likelihood: -564.2994
## AIC: 1268.599
## Structural Break: At Observation Number 59 during 1979 Q3
## Number of GLS estimations: 22
## Number of Restrictions: 0
## 
## Estimated unconditional Heteroscedasticity Matrix (Lambda):
##         [,1]     [,2]     [,3]
## x  0.3925906 0.000000 0.000000
## pi 0.0000000 0.191641 0.000000
## i  0.0000000 0.000000 1.244348
## 
## Standard Errors of Lambda:
##         [,1]       [,2]      [,3]
## x  0.0926582 0.00000000 0.0000000
## pi 0.0000000 0.04527264 0.0000000
## i  0.0000000 0.00000000 0.2935572
## 
## Estimated B Matrix (unique decomposition of the covariance matrix): 
##           [,1]       [,2]      [,3]
## x   0.61193300 -0.5931964 0.2241237
## pi  0.75559400  1.2987520 0.1131134
## i  -0.02899916  0.1572953 0.7084709
## 
## Standard Errors of B:
##         [,1]      [,2]       [,3]
## x  0.1330924 0.1955350 0.07101215
## pi 0.2498465 0.2600375 0.09960246
## i  0.1559672 0.1213446 0.07004431
## 
## Identification Wald Test of equal Eigenvalues:
## [1] 1.2443485 0.3925906 0.1916410
##                              Test statistic dof p-value  
## lambda_ 1 =lambda_2                  5.3828   2 0.06779 .
## lambda_ 1 =lambda_2=lambda_3        15.0586   5 0.01011 *
## lambda_ 2 =lambda_3                  2.1465   2 0.34189  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
Si el VAR puede ser escrito con matrices y como un sistema de media móviles (como se verá más adelante) de la siguiente forma:

\begin{align}
x_{t}&=B(L)\epsilon_{t},\;\;B(0)=I,\;\;E(\epsilon_{t}\epsilon_{t}^{'})=\Sigma 
\end{align}   
Los resultados muestran la matriz de covarianza estimada $\hat{B}$, la matriz de cambio en la convarianza estiamda $\hat{\lambda}$.

Es de notar que el orden de las columnas en  $\hat{B}$ es arbitrario. Lo importante es que ese ordenamiento se haga teniendo un sentido económico.

A pesar de lo anterior, @Herwartz2016 hacen un ordenamiento en las de acuerdo con el patrón de signo único que indica la dirección
de los choques. El siguiente código ordena las columnas de la misma manera. 

```r
usa.cv$B <- usa.cv$B[, c(3, 2, 1)]
usa.cv$B[,3] <- usa.cv$B[, 3] * (-1)
usa.cv$B_SE <- usa.cv$B_SE[, c(3, 2, 1)]
usa.cv$Lambda <- diag(diag(usa.cv$Lambda)[c(3, 2, 1)])
usa.cv$Lambda_SE <- diag(diag(usa.cv$Lambda_SE)[c(3, 2, 1)])
round(usa.cv$B, 3)
```

```
##     [,1]   [,2]   [,3]
## x  0.224 -0.593 -0.612
## pi 0.113  1.299 -0.756
## i  0.708  0.157  0.029
```

```r
round(usa.cv$Lambda, 3)
```

```
##       [,1]  [,2]  [,3]
## [1,] 1.244 0.000 0.000
## [2,] 0.000 0.192 0.000
## [3,] 0.000 0.000 0.393
```


## Función Impulso Respuesta

La función impulso respuesta es la senda que sigue una serie cuando enfrenta un shock unitario.

Es interesantes por:

1. Caracteriza el comportamiento del modelo.
2. Permite la discusión de quién causa a quién.


Para una proceso AR(1) $y_{t}=\phi y_{t-1}+\epsilon_{t}$ el cual lo expresamos en MA($\infty$) como $y_{t}=\sum_{j=0}^{\infty}\phi^{j}  \epsilon_{t-j}$, su impulso respuesta es:


$\begin{array}{cccccccc}
\epsilon_{t} &\colon &0&0&1&0&0&0 \\
x_{t}        &\colon &0&0&\phi&\phi^{2}&\phi^{3}&\cdots
\end{array}$

Similarmente para un proceso MA($\infty$): $y_{t}=\sum_{j=0}^{\infty}\theta \epsilon_{t-j}$, su impulso respuesta es:

$\begin{array}{cccccccc}
	\epsilon_{t} &\colon &0&0&1&0&0&0 \\
	x_{t}        &\colon &0&0&\theta&\theta_{2}&\theta_{3}&\cdots
\end{array}$

¿Por qué es importante definir un proceso como un MA$(\infty)$?

1. La representación MA infinito de todo proceso es *la función impulso respuesta*
2. La función impulso respuesta es equivalente a $E_{t}(x_{t+h})-E_{t-1}(x_{t+h})$

 
Representación de un proceso ARMA como un vector AR(1)

Supongamos que tenemos un proceso ARMA(2,1):


$y_t = \phi_1y_{t-1} + \phi_2y_{t-2} + \epsilon_t + \theta_1\epsilon_{t-1}$



Podemos reescribirlo como:
\begin{equation}
\left( \begin{array}{c}
y_{t} \\
y_{t-1} \\
\epsilon_{t} 
\end{array}
\right)=\left( \begin{array}{ccc}
\phi_1 & \phi_2 & \theta_1\\
1 & 0 & 0\\
0 & 0 & 0
\end{array}
\right) \left( \begin{array}{c}
y_{t-1} \\
y_{t-2} \\
\epsilon_{t-1} 
\end{array}
\right)+\left( \begin{array}{c}
1\\
0 \\
1
\end{array}
\right) \epsilon_{t} 
\nonumber
\end{equation} 


Lo cual podemos escribir en forma AR(1) como:

$x_t = Ax_{t−1} + Cw_t$



## Función Impulso Respuesta en una representación VAR

Similarmente en VAR: $x_t = Ax_{t−1} + Cw_t$, su impulso respuesta es:

$C, \;\; AC, \;\; A^{2}C, \;\; \cdots, \;\;A^{k}C, \;\;\cdots$

## Ortogonalización

La función impulso respuesta de un VAR es ligeramente ambigua.

Sabemos que podemos representar cualquier serie de tiempo como una arbitraria combinación lineal de un conjunto de funciones de impulso respuestas.

Entonces la *Ortogonalización* se refiere al proceso de seleccionar uno de las posibles funciones de impulsos respuestas que sea más interesante para el análisis económico. 

Comenzamos con un VAR expresado en notación vectorial, siendo el vector $x_{t}$ función de sus rezagos:

\begin{align}
A(L)x_{t}&=\epsilon_{t},\;\;A(0)=I,\;\;E(\epsilon_{t}\epsilon_{t}^{'})=\Sigma 
\end{align}    

En notación MA se escribe:

\begin{align}
x_{t}&=B(L)\epsilon_{t},\;\;B(0)=I,\;\;E(\epsilon_{t}\epsilon_{t}^{'})=\Sigma 
\end{align}   



Donde $B(L)=A(L)^{-1}$, siendo $B(L)$ la respuesta de $x_{t}$ a una unidad de impulso de cada elemento de $\epsilon_{t}$    

Ahora supongase que queremos calcular la respuesta de $x_{t}$ a un nuevo shock, el cual es una combinación lineal de viejos schocks, por ejemplo usted quiere calcular el efecto de $x_{t}$ de un schocks tal que: $\eta_{1t}=\epsilon_{yt}$ y  $\eta_{2t}=5\epsilon_{yt}+\epsilon_{zt}$ ó:
\begin{equation}
\eta_{t}=Q\epsilon_{t},\;\;Q=\left( \begin{array}{cc}
    1 & 0\\
    5 & 1
\end{array}
\right) \nonumber
\end{equation} 

Con ello podemos reescribir el VAR en notación MA:
\begin{align}
x_{t}&=C(L)\eta_{t},\;\;C(L)=B(L)Q^{-1}
\end{align}

Donde C(L) representa la respuesta de $x_{t}$ a los schocks $\eta_{t}$, ó visto de una forma alternativa C(L) sería una combinación lineal de los impulsos respuestas originales B(L).

Pero ¿Qué combinación lineal deberíamos observar?

Es de notar que los datos no nos pueden ayudar ya que tanto la representación son equivalentes y producen las mismas series. Así que somos nosotros quienes decidimos cual combinación lineal es más interesante, en base a un conjunto de supuestos llamados *supuestos de ortogonalización*.

## Supuestos de Ortogonalización

El primero y quizás el más popular supuesto es que los shocks son ortogonales  (no están correlacionados).

Si dos shocks $\epsilon_{yt}$ y $\epsilon_{zt}$ están correlacionados entonces no tiene sentido preguntarse, ¿qué pasa si $\epsilon_{yt}$ se impulsa en una unidad?, pues obviamente también se estaría moviendo $\epsilon_{zt}$  al mismo tiempo.

Un ejemplo concreto, nosotros queremos saber la función impulso respuesta en términos causales, es decir el *efecto* de las remesa en el PIB. Pero si el shock de las remesas está correlacionado al shock del PIB, no podemos saber si lo que estamos viendo como respuesta del PIB a un shock de remesas es debido a un shock en el PIB del país de origen de las remesas.
    
Adicionalmente es conveniente reescalar el shock para que este se expresa en unidades de la varianza de las series.\\
El supuesto de shocks ortogonales busca encontrar una matriz Q, tal que $E(\eta_{t}\eta_{t}^{'})=I$. Por lo cual:

\begin{align}
Q^{-1}Q^{-1'}&=\Sigma \\
E(\eta_{t}\eta_{t}^{'})&=E(Q\epsilon_{t}\epsilon_{t}^{'}Q^{'})=Q\Sigma Q^{'}=I\nonumber
\end{align}    



## Ortogonalización de Sims

Una vía de construir la matriz Q es por la descomposición de Choleski.\\
Desafortunadamente, existe muchas matrices Q que pueden hacer cumplir la identidad (22), entonces la pregunta es: ¿cuál de esas matrices Q hay que escoger?: La respuesta dependerá de las propiedades que impongamos del proceso MA ó C(L), en concreto se utilizara la teoría económica para especificar C(0) y C(1).

Sims sugiere que se especifica las propiedades de C(0).

Enfatizando el hecho de que sólo en el caso en que $\Sigma$ sea diagonal, cualquier diagonalización de la matriz Q debería tener elementos en la posición no diagonal y por tanto $C(0)\neq I$, implicando que algún shocks debería tener efecto contemporáneo en más de una variable.

En especifico Sims sugiere escoger una matriz triangular inferior C(0) tal que:

\begin{equation}
	\left( \begin{array}{c}
	y_{t} \\
	z_{t} 
	\end{array}\right)=\left( \begin{array}{cc}
	C_{0yy} & 0\\
	C_{0zy} & C_{0zz}
	\end{array}
	\right) \left( \begin{array}{c}
	\eta_{1t}\\
	\eta_{2t}
	\end{array}
	\right) +C_{1}\eta_{t-1}+\cdots
\end{equation} 

La matriz triangular inferior C(0) implica que $y_{t}$  aparece en la ecuación de $z_{t}$, pero $z_{t}$ no en la ecuación de $y_{t}$. Para observar lo anterior podemos representar el sistema dado por (23) como un proceso autorregresivo y ortogonalizado: $D(L)x_{t}=\eta_{t}$, donde $D(L)=C(L)^{-1}$; recordando que la inversa de una matriz triangular inferior es también una matriz triangular inferior, tenemos:

\begin{equation}
      \left( \begin{array}{cc}
      D_{0yy} & 0\\
      D_{0zy} & D_{0zz}
      \end{array}
      \right)
      \left( \begin{array}{c}
      y_{t} \\
      z_{t} 
      \end{array}\right)+D_{1}x_{t-1}+\cdots=\eta_{t}\nonumber
\end{equation}   
      
Ó
       
\begin{align}
D_{0yy}y_{t}&=\;\;\;\;\;\;\;\;\;\;\;\;-D_{1yy}y_{t-1}-D_{1yz}z_{t-1}+\eta_{1t}\nonumber \\
D_{0zz}z_{t}&=D_{0zy}y_{t}-D_{1zy}y_{t-1}-D_{1zz}z_{t-1}+\eta_{2t}\nonumber
\end{align}
       
 
 Es de notar que la sugerencia de Sims es equivalente a estimar un sistema $(y_{t}, z_{t})$ por OLS, en que la ecuación de $z_{t}$ tenga como regresor  $y_{t}$, pero no viceversa y escalando cada ecuación para que la varianza del error sea 1.

\begin{align}
y_{t}&=\;\;\;\;\;\;\;\;\;\;\;\;\;\;\alpha_{1yy}y_{t-1}+\cdots+\alpha_{1yz}z_{t-1}+\cdots+\eta_{1t}\nonumber \\
z_{t}&=\alpha_{0zy}y_{t}+\alpha_{1zy}y_{t-1}+\cdots+\alpha_{1zz}z_{t-1}+\cdots +\eta_{2t}\nonumber
\end{align}

En resumen se puede unicamente especificar la matriz Q como una combinación lineal de los shocks originales $(\epsilon_{t})$ y elaborar los impulsos respuesta de acuerdo a:


1. El error $(\epsilon_{t})$ son ortogonales.
2. Que la respuesta instantánea de una variable a un shock es cero, siendo este al caso de estimar un VAR por OLS con efecto contemporaneo de y sobre z y no viceversa.

Es de notar que si se específica C(0) se garantiza que al aplicar la descomposición de Choleski se produce una única matriz triangular inferior Q:
\begin{align}
C(0)&=B(0)Q^{-1}=Q^{-1} \nonumber\\
\end{align}

Entonces lo importante es decidir el orden de prelación *De lo exógeno a lo endógeno* de las variables en el VAR; auxiliándose de la teoría económica.
    
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



