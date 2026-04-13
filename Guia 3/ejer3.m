## a) Algoritmo de Jacobi
## b) Algoritmo de Gauss Seidel
## d) Justifique si el metodo de Gauss-Seidel convergera para la matriz del ejercicio 1
## Paso 1: Verificar diagonal dominante A = [3/4 1/6; 1/4 1/4]
## Fila 1 se cumple, pero Fila 2 no, por ser iguales. No es estrictamente dominante

## Paso 2: Ver si Gauss-Seidel converge igual
## Ahora vamos a calcular el radio espectral de la matriz iterativa de Gauss-Seidel. Para eso se usa:
## Tgs = -(D + L)^-1 * U
## D: parte diagonal de A
## L: parte estrictamente inferior
## U: parte estrictamente superior

A = [ 3/4 1/6; 1/4 1/4];

D = diag(diag(A));  # parte diagonal de A
L = tril(A, -1);    # parte estrictamente inferior
U = triu(A, 1);     # parte estrictamente superior

T_GS = -inv(D + L) * U;

rho = max(abs(eig(T_GS)));

disp("Radio espectral de T_GS:");
disp(rho);

## Radio espectral de T_GS: 0.2222
## Como el radio espectral < 1 el método de Gauss-Seidel converge para la matriz del
## ejercicio 1, aunque no sea diagonalmente dominante.

[rhoA] = radio_espectral(A, 'gs', []);
disp("Radio espectral de T_GS:");
disp(rhoA);

