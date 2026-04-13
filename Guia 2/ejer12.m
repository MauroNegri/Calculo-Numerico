# Ax = b, A = LU
# LUx = b, y = Ux
# Ly = b

# EJERCICIO 5
disp("Ejercicio 5: ")
A5 = [ 16.87 0.1650 0.2019 0.3170 0.2340 0.1820 0.1100;
         0 27.70 0.8620 0.0620 0.0730 0.1310 0.1200;
         0 0 22.35 13.05 4.420 6.001 0.3710;
         0 0 0 11.28 0 1.110 0.3710;
         0 0 0 0 9.850 1.1684 2.108;
         0 0 0 0 0.2990 15.98 2.107;
         0 0 0 0 0 0 4.670]; # Matriz de coeficientes de sensibilidad (de la tabla)
b5 = [17.1; 65.1; 186.0; 82.7; 84.2; 63.7; 119.7]; # Vector de alturas de picos (dado en el enunciado)

# Factorización LU
[L5, U5] = doolittle(A5);

# Sustituciones
y5 = sust_adelante(L5, b5);
x5 = sust_atras(U5, y5);

# Resultado
disp("Solución x: "), disp(x5)

# Presión total
P_total = sum(x5);
disp("Presión total calculada: "), disp(P_total)
disp("Presión medida en el experimento: 38.78 μm Hg")

error_rel = abs(P_total - 38.78) / 38.78 * 100;
disp("Error relativo (%): "), disp(error_rel)

disp("-----------------------------------------------------------------")
# EJERCICIO 6
disp("Ejercicio 6: ")

A6 = [80 -50 -30 0; -50 100 -10 -25; -30 -10 65 -20; 0 -25 -20 100];
b6 = [120; 0; 0; 0];

# Factorización LU
[L6, U6] = doolittle(A6);

# Sustituciones
y6 = sust_adelante(L6, b6);
x6 = sust_atras(U6, y6);

disp("Solución x: "), disp(x6)

disp("-----------------------------------------------------------------")
# EJERCICIO 7
disp("Ejercicio 7: ")
N = 10;
A7 = zeros(N,N);
Bin = [-ones(N,1), 2*ones(N,1), -ones(N,1)]; # Matriz donde cada columna es una diagonal
d = [-1, 0, 1]; # Vector con los offsets de cada diagonal
A7 = spdiags(Bin, d, N, N);
A7(1,:) = 0;  A7(1,1) = 1;
A7(N,:) = 0;  A7(N,N) = 1;
A7 = full(A7);

b7 = ones(N,1) / N^2;
b7(1) = 0;
b7(N) = 0;

# Factorización LU
[L7, U7] = doolittle(A7);

# Sustituciones
y7 = sust_adelante(L7, b7);
x7 = sust_atras(U7, y7);

disp("Solución x: "), disp(x7)
