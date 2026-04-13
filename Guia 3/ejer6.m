A1 = [3 1 1; 1 3 -1; 3 1 -5];
b1 = [5; 3; -1];
A2 = [3 1 1; 3 1 -5; 1 3 -1];
b2 = [5; -1; 3];
x0 = [0; 0; 0];
tol = 1e-6;
maxit = 1000;

# Radios espectrales
[rho_GS1] = radio_espectral(A1, 'gs', []);
[rho_GS2] = radio_espectral(A2, 'gs', []);

disp("Radio espectral Gauss-Seidel A1:");
disp(rho_GS1);
disp("Radio espectral Gauss-Seidel A2:");
disp(rho_GS2);

# Métodos iterativos
[x_GS1, it_GS1, rh_GS1] = gauss_seidel(A1, b1, x0, maxit, tol);
[x_GS2, it_GS2, rh_GS2] = gauss_seidel(A2, b2, x0, maxit, tol);

# Mostrar resultados Gauss-Seidel A1
disp("Método Gauss-Seidel A1:");
disp(["Iteraciones: ", num2str(it_GS1)]);
disp(["Error final: ", num2str(rh_GS1(end))]);

# Mostrar resultados Gauss-Seidel A2
disp("Método Gauss-Seidel A2:");
disp(["Iteraciones: ", num2str(it_GS2)]);
disp(["Error final: ", num2str(rh_GS2(end))]);

# Método de Gauss
[x1] = gauss(A1,b1);
[x2] = gauss(A2,b2);
[x2_p] = gauss_p(A2,b2);

# Mostrar Resultados de Gauss A1
disp("Método de Gauss A1:");
disp(x1);

# Mostrar Resultados de Gauss A2
disp("Método de Gauss A2:");
disp(x2);

# Mostrar Resultados de Gauss A2 con pivoteo
disp("Método de Gauss A2 con pivoteo:");
disp(x2_p);

##Sistema 1 (A1) — Todo funciona
##- Gauss-Seidel: converge en 11 iteraciones con radio espectral 0.2582 < 1.
##  La matriz A1 es diagonalmente dominante
##  (fila 1: 3 > 1+1, fila 2: 3 > 1+1, fila 3: 5 > 3+1), lo que garantiza convergencia.
##- Gauss sin pivoteo: funciona correctamente porque ningún pivote se anula durante la eliminación.

##Sistema 2 (A2) — Mismo sistema, distinto orden de filas
##- A2 es exactamente A1 con las filas 2 y 3 intercambiadas. La solución es la misma [1 1 1],
##  pero el comportamiento numérico cambia completamente.
##- Gauss-Seidel: diverge, radio espectral 18.577 >> 1. Al reordenar las filas la diagonal
##  queda no dominante (fila 2: |1| < |3|+|5|), rompiendo la condición de convergencia.
##  El error llega a NaN porque los valores explotan.
##- Gauss sin pivoteo: devuelve NaN porque durante la eliminación aparece un pivote
##  cero o casi cero, generando divisiones por cero.
##- Gauss con pivoteo: resuelve correctamente porque el pivoteo reordena las filas
##  buscando el mayor elemento, recuperando esencialmente el orden de A1.

##El orden de las ecuaciones importa. Dos sistemas matemáticamente equivalentes
##pueden comportarse de forma completamente distinta numéricamente.
##Para métodos iterativos como Gauss-Seidel es necesario que la matriz sea
##diagonalmente dominante, lo que depende del orden de las filas.
##Para eliminación de Gauss, el pivoteo parcial es necesario cuando el orden
##natural de las filas produce pivotes nulos o muy pequeños.
