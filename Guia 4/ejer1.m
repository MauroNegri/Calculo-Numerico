## Ejercicio 1.
## a) Realice cuatro iteraciones con el metodo de la biseccion para obtener una aproximacion a
## una de las raices de la ecuación f(x) = 3(x + 1)(x − 0.5)(x − 1) en el intervalo [−2, 1.5].
## ¿A cuál de las raices converge el metodo? Luego estime una cota para la precision del resultado obtenido.
##
## Paso 1: Identificamos el signo de f(x) en los extremos del intervalo
## f(−2)=3(−1)(−2.5)(−3)=3(−1)(−2.5)(−3)=−22.5
## f(1.5)=3(2.5)(1)(0.5)=3.75>0
## Como f(−2)<0 y f(1.5)>0, hay al menos una raíz en el intervalo [−2,1.5]
##
## Paso 2: Hacemos 4 iteraciones del método de bisección
## Iteración 1:
## x1 = (-2 + 1.5) / 2 = −0.25, f(−0.25)=3(0.75)(−0.75)(−1.25) ≈ 2.109
## Cambio de signo está entre [−2,−0.25]
##
## Iteración 2:
## x2 = (-2 + (-0.25)) / 2 = -1.125, f(−1.125)=3(−0.125)(−1.625)(−2.125) ≈ −1.294
## Cambio de signo está entre [−1.125,−0.25]
##
## Iteración 3:
## x3 = (-1.125 + (-0.25)) / 2 = −0.6875, f(−0.6875)≈3(0.3125)(−1.1875)(−1.6875) ≈ 1.875
## Cambio de signo está entre [−1.125,−0.6875]
##
## Iteración 4:
## x4 = (-1.125 + (-0.6875)) / 2 = −0.90625, f(−0.90625)≈3(0.09375)(−1.40625)(−1.90625) ≈ 0.754
## Cambio de signo está entre [−1.125,−0.90625]
##
## Resultado después de 4 iteraciones:
## Aproximación de la raíz:
##
## 𝑥4 = −0.90625
## Cota de error:
## Error ≤ (𝑏 - a) / 2 = (1.5 - (-2)) / 2^4 = 3.5 / 16 ≈ 0.21875
##
## ¿A qué raíz converge?
## La función tiene raíces en
## x=−1, x=0.5, y x=1
## Como el resultado de la bisección se está acercando a x=−1, podemos decir que el método converge a la raíz x=−1.
##
## b) Implemente una funcion de Octave function [x,h] = biseccion(f,xmin,xmax,kmax,tol)
## function [p,h] = biseccion(f,a,b,maxit,tol)
##
## c)  Obtenga una cota para el numero de iteraciones que se requieren para
## alcanzar una aproximacion con una exactitud de 10−3 a la solucion de x3 + x − 4 = 0
## que se encuentra en el intervalo [1,4]. Obtenga una aproximacion de la raiz
## con esta exactitud mediante el metodo de la biseccion.
clc; format long ;
fa = @(x) 3*(x+1).*(x-0.5).*(x-1)
fb = @(x) x.^3 + x - 4;
a = 1;
b = 4;
maxit = 12;
tol = 1e-3;

% Cota: (b-a)/2^n < tol  →  (4-1)/2^n < 1e-3
% 2^n > 3000  →  n > log2(3000) ≈ 11.55  →  n_min = 12
n_cota = ceil(log2((b-a)/tol));
fprintf("Cota de iteraciones necesarias: %d\n", n_cota);

[pa, ~, ha, ~] = biseccion(fa, -2, 1.5, 4, 1e-8);
[pb, ~, hb, ~] = biseccion(fb, a, b, maxit, tol);

fprintf("Raíz aproximada de a: %.6f\n", pa);
fprintf("Error final de a: %.6f\n", ha(end));

fprintf("Raíz aproximada de b: %.6f\n", pb);
fprintf("Error final de b: %.6f\n", hb(end));


## d) Escribir una funcion x=rcubica(a) para calcular la raiz cubica de a con
## un error relativo menor a 10^−12 usando biseccion. En la misma debe hacer
## una llamada a la funcion implementada en (b).

##function x = rcubica(a)
##  tol = 1e-12;
##  maxit = 1000;
##  if a == 0; x = 0; return; endif
##  f = @(t) t.^3 - a;
##  if a > 1
##    xmin = 0; xmax = a;
##  elseif a > 0
##    xmin = 0; xmax = 1;
##  else
##    xmin = a; xmax = 0;
##  endif
##  [x, ~, ~, ~] = biseccion(f, xmin, xmax, maxit, tol);  % usa maxit_d y tol_rel
##endfunction

