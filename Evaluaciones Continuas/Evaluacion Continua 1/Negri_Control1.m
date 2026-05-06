A = [ 3 1 -1 1 2; 1 3 2 -2 -2; -2 -2 1 -2 -2; 2 -1 -1 3 2; 0 0 1 1 2];
b = [ 1; 4; -1; 2; 2];
# a)
# Metodo de Gauss
[x] = gauss(A,b);
disp(x);

# b)
x0 = zeros(5,1);
tol = 1e-4;
maxit = 500;

# Radios espectrales
[rho_JA] = radio_espectral(A,'ja', []);
[rho_GS] = radio_espectral(A,'gs', []);

disp("Radio espectral Jacobi:");
disp(rho_JA);
disp("Radio espectral Gauss-Seidel:");
disp(rho_GS);

# Métodos iterativos
[X_JA, it_JA, rh_JA, t_JA] = jacobi(A, b, x0, maxit, tol);
[X_GS, it_GS, rh_GS, t_GS] = gauss_seidel(A, b, x0, maxit, tol);

# Resultados Jacobi
disp("Método Jacobi:");
disp(["Iteraciones: ", num2str(it_JA)]);
disp(["Error final: ", num2str(rh_JA(end))]);
disp(["Tiempo: ", num2str(t_JA), " segundos"]);

# Resultados Gauss-Seidel
disp("Método Gauss-Seidel:");
disp(["Iteraciones: ", num2str(it_GS)]);
disp(["Error final: ", num2str(rh_GS(end))]);
disp(["Tiempo: ", num2str(t_GS), " segundos"]);

err = norm(x - X_JA, inf) / norm(x, inf);
printf("Error relativo (norma inf) entre Gauss y Jacobi: %.4g\n", err);
disp(["Error relativo (norma inf) entre Gauss y Jacobi: ", num2str(err, '%.4g')])
