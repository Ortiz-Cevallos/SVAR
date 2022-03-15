---
title: "VECTORES AUTOREGRESIVOS (VAR)"
author: "[Luis Ortiz-Cevallos](https://ortiz-cevallos.github.io/MYSELF/)"
date: "2022-03-15"
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

## El modelo de vectores autoregresivos (VAR)

Un VAR consiste en un conjunto de k variables endógenas $Y_{t}=\{y_{1t},\cdots, y_{kt}\}$.

En el caso de un proceso ARMA(2,1):
$Y_t = \Phi_1Y_{t-1} + \Phi_2Y_{t-2} + u_t + \theta_1u_{t-1}$

Podemos reescribirlo como:
\begin{equation}
\left( \begin{array}{c}
Y_{t} \\
Y_{t-1} \\
u_{t} 
\end{array}
\right)=\left( \begin{array}{ccc}
\Phi_1 & \Phi_2 & \theta_1\\
1 & 0 & 0\\
0 & 0 & 0
\end{array}
\right) \left( \begin{array}{c}
Y_{t-1} \\
Y_{t-2} \\
u_{t-1} 
\end{array}
\right)+\left( \begin{array}{c}
1\\
0 \\
1
\end{array}
\right) \left( \begin{array}{c}
u_{t} \\
0 \\
0
\end{array}
\right) 
\nonumber
\end{equation} 


Lo anterior no es otra cosa que un proceso AR(1), matricialmente:

$x_t = Ax_{t−1} + Cw_t$

Este proceso tiene los siguientes supuestos:

1. $E(w_{t})=0$

2. $E(w_{t}w_{t}^{\tau})=\Sigma_{w}$ una matriz de covarianza invariante en el tiempo la cual  es definida positiva.


## Función Impulso Respuesta

La función impulso respuesta es la senda que sigue una serie cuando enfrenta un shock unitario.

Es interesante por:

1. Caracteriza el comportamiento del modelo.
2. Permite la discusión de quién causa a quién.


Para una proceso AR(1) $y_{t}=\phi y_{t-1}+\epsilon_{t}$ el cual lo expresamos en MA($\infty$) como $y_{t}=\sum_{j=0}^{\infty}\phi^{j}  \epsilon_{t-j}$, su impulso respuesta es:


$\begin{array}{cccccccc}
\epsilon_{t} &\colon &0&0&1&0&0&0 \\
x_{t}        &\colon &0&0&\phi&\phi^{2}&\phi^{3}&\cdots
\end{array}$

¿Por qué es importante definir un proceso como un MA$(\infty)$?

1. La representación MA infinito de todo proceso es *la función impulso respuesta*
2. La función impulso respuesta es equivalente a $E_{t}(x_{t+h})-E_{t-1}(x_{t+h})$ que no es otra cosa que el error de pronóstico hacia h pasos adelante.

Noten que los pronósticos para horizontes $h\geq1$ de un proceso VAR(p) puede ser generado de manera recursiva.
$x_{t+h/t} = A_{1}x_{t+h-1/t} + \cdots+A_{p}x_{t+h-p/t}$


## Vectores autoregresivos estructurales

Casi todas las variables en economía podrían considerarse 
hasta cierto punto endógenas. Sin embargo, prima un criterio de relevancia. Si bien es cierto que un comerciante al fijar el precio de una manzana está afectando la inflación, parece razonable suponer que la inflación es exógena al precio de las manzanas de este comerciante.

Por tanto si la inflación está dado por $\pi$ y el precio de las manzanas por $m$ tenemos el siguiente VAR

\begin{align}
\alpha_{1,1}^{+}\pi_t+\alpha_{1,2}^{+}m_t &= \beta_{1,1}^{+}\pi_{t-1}+\beta_{1,2}^{+}m_{t-1}+e_{\pi,t}\\
\alpha_{2,1}^{+}\pi_t+\alpha_{2,2}^{+}m_t &= \beta_{2,1}^{+}\pi_{t-1}+\beta_{2,2}^{+}m_{t-1}+e_{m,t}
\end{align}

Sí normalizamos el VAR tenemos lo siguiente:
\begin{align}
\pi_t+\alpha_{1,2}m_t &= \beta_{1,1}\pi_{t-1}+\beta_{1,2}m_{t-1}+e_{\pi,t}\\
\alpha_{2,1}\pi_t+m_t &= \beta_{2,1}\pi_{t-1}+\beta_{2,2}m_{t-1}+e_{m,t}
\end{align}

Si resolvemos el sistema de ecuaciones simultáneas, despejamos las
variables obtenemos un VAR en su forma reducida; comenzamos con la primera ecuación:

\begin{align}
\pi_{t} &=\beta_{1,1}\pi_{t-1}+\beta_{1,2}m_{t-1}+e_{\pi,t}\\
        &- \alpha_{1,2}(\beta_{2,1}\pi_{t-1}+\beta_{2,2}m_{t-1}-\alpha_{2,1}\pi_t+e_{m,t})\\
        \pi_{t} &=\frac{(\beta_{1,1}-\alpha_{1,2}\beta_{2,1})}{(1+\alpha_{1,2}\alpha_{2,1})}\pi_{t-1}+ \frac{(\beta_{1,2}-\alpha_{1,2}\beta_{2,2})}{(1+\alpha_{1,2}\alpha_{2,1})}m_{t-1}+\frac{1}{(1+\alpha_{1,2}\alpha_{2,1})}e_{\pi,t}- \frac{\alpha_{1,2}}{(1+\alpha_{1,2}\alpha_{2,1})}e_{m,t}\\
 \pi_{t}&=\gamma_{\pi}\pi_{t-1}+\gamma_{m}m_{t-1}+\gamma_{e_{\pi}}e_{\pi,t}+\gamma_{e_{m}}e_{m,t}       
\end{align}

De manera análoga para las manzanas:
\begin{align}
m_{t} &=\beta_{2,1}\pi_{t-1}+\beta_{2,2}m_{t-1}+e_{m,t}\\
        &- \alpha_{2,1}(\beta_{1,1}\pi_{t-1}+\beta_{1,2}m_{t-1}-\alpha_{1,2}m_t+e_{\pi,t})\\
        m_{t} &=\frac{(\beta_{2,1}-\alpha_{2,1}\beta_{1,1})}{(1+\alpha_{1,2}\alpha_{2,1})}\pi_{t-1}+ \frac{(\beta_{2,2}-\alpha_{2,1}\beta_{1,2})}{(1+\alpha_{1,2}\alpha_{2,1})}m_{t-1}+\frac{\alpha_{2,1}}{(1+\alpha_{1,2}\alpha_{2,1})}e_{\pi,t}- \frac{1}{(1+\alpha_{1,2}\alpha_{2,1})}e_{m,t}\\
 m_{t}&=\omega_{\pi}\pi_{t-1}+\omega_{m}m_{t-1}+\omega_{e_{\pi}}e_{\pi,t}+\omega_{e_{m}}e_{m,t}       
\end{align}


Entonces ahora tenemos el siguiente sistema:
\begin{align}
\pi_{t}&=\gamma_{\pi}\pi_{t-1}+\gamma_{m}m_{t-1}+\gamma_{e_{\pi}}e_{\pi,t}+\gamma_{e_{m}}e_{m,t}\\   
 m_{t}&=\omega_{\pi}\pi_{t-1}+\omega_{m}m_{t-1}+\omega_{e_{\pi}}e_{\pi,t}+\omega_{e_{m}}e_{m,t}       
\end{align}

El sistema anterior puede escribirse de forma matricial como:
\begin{equation}
	\left( \begin{array}{c}
	\pi_{t} \\
	  m_{t} 
	\end{array}\right)=\left( \begin{array}{cc}
 \gamma_{\pi} & \gamma_{m}\\
 \omega_{\pi} & \omega_{m}
	\end{array}
	\right) \left( \begin{array}{c}
	\pi_{t-1}\\
	  m_{t-1}
	\end{array}
	\right) +
\left( \begin{array}{cc}
 \gamma_{e_{\pi}} & \gamma_{e_{m}}\\
 \omega_{e_{\pi}} & \omega_{e_{m}}
	\end{array}
	\right) \left( \begin{array}{c}
	e_{\pi,t}\\
	e_{m,t} 
	\end{array}
	\right) 
\end{equation} 

Expresándolo de forma matricial tenermos:

\begin{equation}
\label{e1}
x_t = \Gamma x_{t−1} + Zw_t
\end{equation}


La ecuación anterior es la *forma estructural*. Sin embargo, cuando se estima el VAR se hace en su *forma reducida*, la que resulta trás premultiplicar por la matriz A:
\begin{equation}
\label{e2}
Ax_t = A\Gamma x_{t−1} + Bw_t
\end{equation}


Si asumimos que los *errores estructurales* que integran el vector $w_{t}$ son ruido blanco, y que los coeficientes en la matriz $\Gamma$ son *estructurales* los cuales difieran de sus contrapartes en su forma reducida. 

Para ver lo anterior, multiplicamos la ecuación anterior por $A^{-1}$:

\begin{align}
x_t &= A^{-1}A\Gamma x_{t−1} + A^{-1}Bw_t\\
x_t &= \Gamma x_{t−1} + \eta_t
\end{align}

Un modelo SVAR puede ser usado para identificar los shocks y trazar su función impulso respuesta o bien analizar su descomposición de varianza, lo anterior gracias a que podemos imponer ciertas restricciones en las matrices A y/o B. 

Es de recalcar que si el modelo SVAR es un modelo estructural, parte de un modelo VAR(p) en su forma reducida:

\begin{equation}
Ax_t = A\Gamma x_{t−1} + Bw_t
\end{equation}

Donde solo se pueden agregar restricciones para A y B. Además los residuos en su  forma reducida se puede recuperar de un modelo SVAR por $\eta_t=A^{-1}Bw_t=Qw_t$ y su matriz de varianza-covarianza es $\Sigma_\eta=A^{-1}B\Sigma_{w}B^{\tau}(A^{-1})^{\tau}=Q\Sigma_{w}Q^{\tau}=QQ^{\tau}$.

De manera que el impulso respuesta de este modelo ante un shock estrutural en ($w_t$) viene dado por:

$\begin{array}{ccccccccc}
w_{t} &\colon &0&0&I&0&0&0&0 \\
x_{t}        &\colon &0&0&A^{-1}B&B&AB&A^{2}B&\cdots
\end{array}$

Es posible distinguir tres tipos de SVAR según se impongan las restricciones:

1. *Modelo A* La matriz B corresponde a $I_{k}$ lo que implica que el número mínimo de restricciones para su identificación es $k\frac{(k-1)}{2}$

2. *Modelo B* La matriz A corresponde a $I_{k}$ lo que implica que el número mínimo de restricciones para su identificación es $k\frac{(k-1)}{2}$

3. *Modelo AB* Las restricciones pueden ser impuestas en ambas matrices siendo el número mínimo de restricciones para su identificación es $k^{2}+k\frac{(k-1)}{2}$

## Ortogonalización

La función impulso respuesta de un VAR es ligeramente ambigua.

Sabemos que podemos representar cualquier serie de tiempo como una *arbitraria* combinación lineal de un conjunto de funciones de impulso respuestas.

Entonces la *Ortogonalización* se refiere al proceso de seleccionar uno de las posibles funciones de impulsos respuestas que sea más interesante para el análisis económico. 

Si definimos la matriz $A^{*}$ en función del operador de rezagos:
\begin{equation}
A^{*}=\left( \begin{array}{cccc}
L      & NA  & \cdots & NA\\
\beta_{1,2} & L   &        & NA\\
\vdots &     & \ddots & NA\\
\beta_{1,k} & \beta_{2,k}   & \cdots & L\\
\end{array}
\right)
\end{equation}

Un VAR expresado en notación vectorial estaría dado:

\begin{align}
A^{*}(L)x_{t}&=w_{t},\;\;A^{*}(0)=I,\;\;E(w_{t}w_{t}^{'})=\Sigma_w 
\end{align}    

Si definimos la matriz B en función de sus rezagos:
\begin{equation}
B^{*}=\left( \begin{array}{cccc}
L      & NA  & \cdots & NA\\
\Phi_{1,2} & L   &        & NA\\
\vdots &     & \ddots & NA\\
\Phi_{1,k} & \Phi_{2,k}   & \cdots & L\\
\end{array}
\right)
\end{equation}


En notación MA se escribe:

\begin{align}
x_{t}&=B^{*}(L)w_{t},\;\;B^{*}(0)=I,\;\;E(w_{t}w_{t}^{'})=\Sigma_w 
\end{align}   



Donde $B^{*}(L)=A^{*}(L)^{-1}$, siendo $B^{*}(L)$ la respuesta de $x_{t}$ a una unidad de impulso de cada elemento de $e_{t}$    

Ahora supongase que queremos calcular la respuesta de $x_{t}$ a un nuevo shock, el cual es una combinación lineal de viejos schocks, por ejemplo usted quiere calcular el efecto de $x_{t}$ de un schocks tal que: $\eta_{1t}=e_{\pi t}$ y  $\eta_{2t}=5e_{\pi t}+e_{mt}$ ó:
\begin{equation}
\eta_{t}=Qw_t,\;\;Q=\left( \begin{array}{cc}
    1 & 0\\
    5 & 1
\end{array}
\right) \nonumber
\end{equation} 

Con ello podemos reescribir el VAR en notación MA:
\begin{align}
x_{t}&=C(L)\eta_{t},\;\;C(L)=B^{*}(L)Q^{-1}
\end{align}

Donde C(L) representa la respuesta de $x_{t}$ a los schocks $\eta_{t}$, ó visto de una forma alternativa C(L) sería una combinación lineal de los impulsos respuestas originales $B^{*}(L)$.

Pero ¿Qué combinación lineal deberíamos observar?

Es de notar que los datos no nos pueden ayudar ya que tanto la representación $x_{t}=B^{*}(L)w_{t}$ y $x_{t}=C(L)\eta_{t}$ son equivalentes y producen las mismas series. 

Así que somos nosotros quienes decidimos cual combinación lineal es más interesante, en base a un conjunto de supuestos llamados *supuestos de ortogonalización*.

## Supuestos de Ortogonalización

El primero y quizás el más popular supuesto es que los shocks son ortogonales  (no están correlacionados).

Si dos shocks $e_{\pi t}$ y $e_{mt}$ están correlacionados entonces no tiene sentido preguntarse, ¿qué pasa si $e_{\pi t}$ se impulsa en una unidad?, pues obviamente también se estaría moviendo $e_{mt}$  al mismo tiempo.

Un ejemplo concreto, nosotros queremos saber la función impulso respuesta en términos causales, es decir el *efecto* de la inflación en el precio de las manzanas. Pero si el shock de la inflación está correlacionado al shock del precio de las manzanas, no podemos saber si lo que estamos viendo como respuesta del precio de las manzanas  se debe a un shock de inflación ó es debido a un shock en la  oferta de las manzanas.
    
Adicionalmente es conveniente reescalar el shock para que este se exprese en unidades de la varianza de las series.

El supuesto de shocks ortogonales busca encontrar una matriz Q, tal que $E(\eta_{t}\eta_{t}^{'})=I$. Por lo cual:

\begin{align}
Q^{-1}Q^{-1'}&=\Sigma_{\eta} \\
E(\eta_{t}\eta_{t}^{'})&=E(Qw_{t}w_{t}^{'}Q^{'})=Q\Sigma_{w} Q^{'}=I\nonumber
\end{align}    



## Ortogonalización de Sims

Una vía de construir la matriz Q es por la descomposición de Choleski.

Desafortunadamente existen muchas matrices Q que pueden hacer cumplir la identidad $E(\eta_{t}\eta_{t}^{'})=I$, entonces la pregunta es: ¿cuál de esas matrices Q hay que escoger?: La respuesta dependerá de las propiedades que impongamos del proceso MA ó C(L), en concreto se utilizara la teoría económica para especificar C(0) y C(1).

Sims sugiere que se especifica las propiedades de C(0).

Enfatizando el hecho de que sólo en el caso en que $\Sigma_{w}$ sea diagonal, cualquier diagonalización de la matriz Q debería tener elementos en la posición no diagonal y por tanto $C(0)\neq I$, implicando que algún shocks debería tener efecto contemporáneo en más de una variable.

En especifico Sims sugiere escoger una matriz triangular inferior C(0) tal que:

\begin{equation}
	\left( \begin{array}{c}
	\pi_{t} \\
	  m_{t} 
	\end{array}\right)=\left( \begin{array}{cc}
	C_{0\pi \pi} & 0\\
	C_{0m  \pi} & C_{0mm}
	\end{array}
	\right) \left( \begin{array}{c}
	\eta_{1t}\\
	\eta_{2t}
	\end{array}
	\right) +C_{1}\eta_{t-1}+\cdots
\end{equation} 

La matriz triangular inferior C(0) implica que $\pi_{t}$  aparece en la ecuación de $m_{t}$, pero $m_{t}$ no en la ecuación de $\pi_{t}$. Para observar lo anterior podemos representar el sistema dado por la ecuación anterior como un proceso autorregresivo y ortogonalizado: $D(L)x_{t}=\eta_{t}$, donde $D(L)=C(L)^{-1}$; recordando que la inversa de una matriz triangular inferior es también una matriz triangular inferior, tenemos:

\begin{equation}
      \left( \begin{array}{cc}
      D_{0\pi m} & 0\\
      D_{0m \pi} & D_{0mm}
      \end{array}
      \right)
      \left( \begin{array}{c}
      \pi_{t} \\
       m_{t} 
      \end{array}\right)+D_{1}x_{t-1}+\cdots=\eta_{t}\nonumber
\end{equation}   
      
 Es de notar que la sugerencia de Sims es equivalente a estimar un sistema $(\pi_{t}, m_{t})$ por OLS, en que la ecuación de $m_{t}$ tenga como regresor  $\pi_{t}$, pero no viceversa y escalando cada ecuación para que la varianza del error sea 1.

\begin{align}
\pi_{t}&=\;\;\;\;\;\;\;\;\;\;\;\;\;\;\alpha_{1\pi \pi}\pi_{t-1}+\cdots+\alpha_{1\pi m}m_{t-1}+\cdots+\eta_{1t}\nonumber \\
m_{t}&=\alpha_{0m \pi}\pi_{t}+\alpha_{1m \pi}\pi_{t-1}+\cdots+\alpha_{1mm}m_{t-1}+\cdots +\eta_{2t}\nonumber
\end{align}

En resumen se puede unicamente especificar la matriz Q como una combinación lineal de los shocks originales $(w_{t})$ y elaborar los impulsos respuesta de acuerdo a:


1. El error $(w_{t})$ son ortogonales.
2. Que la respuesta instantánea de una variable a un shock es cero, siendo este al caso de estimar un VAR por OLS con efecto contemporaneo de $\pi$ sobre m y no viceversa.

Es de notar que si se específica C(0) se garantiza que al aplicar la descomposición de Choleski se produce una única matriz triangular inferior Q:
\begin{align}
C(0)&=B^{*}(0)Q^{-1}=A^{-1}BQ^{-1}=A^{-1}BQ \nonumber\\
\end{align}

Entonces lo importante es decidir el orden de prelación *De lo exógeno a lo endógeno* de las variables en el VAR; auxiliándose de la teoría económica.


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

A pesar de lo anterior, @Herwartz2016 hacen un ordenamiento de acuerdo con el patrón de signo único que indica la dirección de los choques. El siguiente código ordena las columnas de esa misma manera. 

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

@Herwartz2016 interpretan el efecto de la primera columna de la matriz $\hat{B}$ como un shock de demanda. La segunda y tercera columna representan, respectivamente, un shock de oferta y política monetaria.


```r
round(usa.cv$Lambda, 3)
```

```
##       [,1]  [,2]  [,3]
## [1,] 1.244 0.000 0.000
## [2,] 0.000 0.192 0.000
## [3,] 0.000 0.000 0.393
```
En cambio la matriz $\hat{\lambda}$ representa la varianza de los shocks estructurales en un segundo regimén, esto es posterior al quiebre estructural. Los autores interpretan a partir de la diagonal de $\hat{\lambda}$ que los shocks de oferta y política monetaria tiene baja varianza en comparación a los shocks de demanda.

Los autores comparan estos resultados con los obtenidos por @Sims80 utilizando la descomposición de la matriz de covarianza ($\hat{B}$) como una matriz triangular inferior. La función *id.cv()* permite
probar tales restricciones configurando para ello una matriz triangular inferior como se describe en el código siguiente.

```r
restMat <- matrix(rep(NA, 9), ncol = 3)
restMat[1, c(2, 3)] <- 0
restMat[2, 3] <- 0
restMat
```

```
##      [,1] [,2] [,3]
## [1,]   NA    0    0
## [2,]   NA   NA    0
## [3,]   NA   NA   NA
```

```r
restricted.model <- id.cv(plain.var, SB = c(1979, 3),restriction_matrix = restMat)
summary(restricted.model)
```

```
## 
## Identification Results
## ---------------------- 
## 
## Method: Changes in Volatility
## Sample size: 169
## Log-Likelihood: -568.6664
## AIC: 1277.333
## Structural Break: At Observation Number 59 during 1979 Q3
## Number of GLS estimations: 23
## Number of Restrictions: 3
## 
## Estimated unconditional Heteroscedasticity Matrix (Lambda):
##         [,1]      [,2]      [,3]
## x  0.3501948 0.0000000 0.0000000
## pi 0.0000000 0.2346854 0.0000000
## i  0.0000000 0.0000000 0.9420116
## 
## Standard Errors of Lambda:
##          [,1]       [,2]     [,3]
## x  0.08266738 0.00000000 0.000000
## pi 0.00000000 0.05616318 0.000000
## i  0.00000000 0.00000000 0.227189
## 
## Estimated B Matrix (unique decomposition of the covariance matrix): 
##          [,1]      [,2]      [,3]
## x  0.87988465 0.0000000 0.0000000
## pi 0.08137972 1.5306503 0.0000000
## i  0.31518384 0.2606745 0.7378484
## 
## Standard Errors of B:
##          [,1]       [,2]       [,3]
## x  0.08638851 0.00000000 0.00000000
## pi 0.10334026 0.15169565 0.00000000
## i  0.08527441 0.08620187 0.07354585
## 
## Identification Wald Test of equal Eigenvalues:
## [1] 0.9420116 0.3501948 0.2346854
##                              Test statistic dof p-value
## lambda_ 1 =lambda_2                 4.04784   2  0.1321
## lambda_ 1 =lambda_2=lambda_3        9.15413   5  0.1031
## lambda_ 2 =lambda_3                 0.68408   2  0.7103
## 
## Likelihood Ratio Test: 
##  Test statistic p-value  
##           8.734   0.033 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Si partimos que los shocks structurales fueron identificados con la metodología de cambio en la matriz de covarianza, cualquier restricción impuesta sobre $B$ hace que el modelo esté sobre identificado y ello es posible probarse.

En resumen el "likelihood ratio test" indica que bajo la $H_{0}$
que $B$ la matriz triangular (baja la cual impactan los shocks structurales) es rechazada al $5\%$ de nivel de significancia. @Herwartz2016 argumentan que la identificación de valores iguales a cero para lograr una matriz trinagular inferior es algo contra intuitivo desde el punto de vista económico; por tanto un modelo sin restricciones debería ser preferido.

El siguiente paso es calcular la función impulso respuesta (IFR) con un intervalo de confianza obtenido a partir del método "boostrap". A partir de la IFR es posible investigar los efectos futuros de los shocks estructurales etiquetados como "económicos" sobre las variables
incluidas en el modelo. De hecho de acuerdo a @Herwartz18 aplicar la metodología "boostrap permite evaluar la significancia en los signos que se hacen patentes en la matriz $\hat{B}$.

Para aplicar lo anterior, definimos una lista con los signos que representan nuestras restricciones tanto para los shocks de demanda, oferta y política monetaria.


```r
signrest <- list(demand = c(1, 1, 1), supply = c(-1, 1, 1), monetary_policy = c(-1, -1, 1))
```

Notal que el horizonte de tiempo para el IFR tiene que ser determinado de antemano usando el argumento n.ahead.


```r
cores <- parallel::detectCores() - 1
set.seed(231)
usa.cv.boot <- wild.boot(usa.cv, design = "fixed",distr = "rademacher", nboot = 1000, n.ahead = 15,nc = cores, signrest = signrest)
summary(usa.cv.boot)
```

```
## 
## Bootstrap Results
## ----------------- 
## 
## Method: Wild bootstrap
## Bootstrap iterations: 1000
## Distribution used: rademacher
## Design: fixed
## 
## Point estimates: 
##         [,1]       [,2]        [,3]
## x  0.2241237 -0.5931964 -0.61193300
## pi 0.1131134  1.2987520 -0.75559400
## i  0.7084709  0.1572953  0.02899916
## 
## Bootstrap means: 
##          [,1]       [,2]       [,3]
## x  0.09580740 -0.5030485 -0.6303868
## pi 0.09639712  1.1310905 -0.7720157
## i  0.70100321  0.0364798 -0.1675659
## 
## Bootstrap standard errors: 
##          [,1]      [,2]      [,3]
## x  0.14143383 0.3077349 0.2494872
## pi 0.17500155 0.4556690 0.5925734
## i  0.07511879 0.2273595 0.2178134
## 
## Identified sign patterns: 
## =========================
## Specified sign pattern: 
## 
##    demand supply monetary_policy
## x       1     -1              -1
## pi      1      1              -1
## i       1      1               1
## 
## Unique occurrence of single shocks according to sign pattern: 
## demand : 65.5 % 
## supply : 66 % 
## monetary_policy : 28.4 % 
## 
## Joint occurrence of specified shocks: 12.7 %
```

Enseguida, graficamos para cada variable y shocks su IFR.


```r
plot(usa.cv.boot, lowerq = 0.16, upperq = 0.84)
```

![](01-conceptos_files/figure-epub3/unnamed-chunk-12-1.png)<!-- -->

En resumen solamente el 12.7% de todos los bootstrap estimados están en línea con el patron de signos económicamente sensatos.
El signo del patron de shocks monetario esperado aparece en un 28.4%. Finalmente, el bootstrap muestra que el tercer shock va más en línea con el patrón asociado a un shock de demanda.


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



