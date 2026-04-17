clc;
A = [0.5 -1 0; -1 3 -1; 0 -1 2];
b = [7; 4; 5];

x0 = b;
maxit = 1000;
tol = 1e-8;

% Radios espectrales
[rho_JA] = radio_espectral(A, 'ja', []);
[rho_GS] = radio_espectral(A,'gs', []); % Calculo el Radio Espectral de GaussSeidel y con eso saco w

disp("Radio espectral Jacobi:");
disp(rho_JA);
disp("Radio espectral Gauss-Seidel:");
disp(rho_GS);

% Parámetro óptimo w
w_opt = 2 / (1 + sqrt(1 - rho_JA^2));
disp("Valor óptimo de w:");
disp(w_opt);

% Métodos iterativos
[X_JA, it_JA, rh_JA, t_JA] = jacobi(A, b, x0, maxit, tol);
[X_GS, it_GS, rh_GS, t_GS] = gauss_seidel(A, b, x0, maxit, tol);
[X_SOR, it_SOR, rh_SOR, t_SOR] = sor(A,b, x0, maxit, tol, w_opt);

% Mostrar resultados Jacobi
disp("Método Jacobi:");
disp(["Iteraciones: ", num2str(it_JA)]);
disp(["Error final: ", num2str(rh_JA(end))]);
disp(["Tiempo: ", num2str(t_JA), " segundos"]);

% Mostrar resultados Gauss-Seidel
disp("Método Gauss-Seidel:");
disp(["Iteraciones: ", num2str(it_GS)]);
disp(["Error final: ", num2str(rh_GS(end))]);
disp(["Tiempo: ", num2str(t_GS), " segundos"]);

% Mostrar resultados SOR
disp("Método SOR:");
disp(["Iteraciones: ", num2str(it_SOR)]); % printf("Iteraciones: %d\n", itSOR); Más preciso y limpio
disp(["Error final: ", num2str(rh_SOR(end))]);
disp(["Tiempo: ", num2str(t_SOR), " segundos"]);
