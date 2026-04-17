n = 40;
%main_diag = 2 * (1:n)';
%B = [ -1*ones(n,1), 2*(1:n)', -1*ones(n,1)];
%d = [-1, 0, 1];

main_diag = 2 * (1:n)';
sub_diag = -1 * ones(n, 1);
super_diag = -1 * ones(n, 1);

B = [sub_diag, main_diag, super_diag];
d = [-1, 0, 1];


A = spdiags(B,d,n,n);
A=full(A);
b = 1.5 * (1:n)' - 6;
x0 = b; % Comenzando las iteraciones con el vector b
tol = 1e-6;
maxit = 100;

% Método de Gauss
##[A_tri1, b_tri1] = gauss(A,b);
[x] = gauss(A,b);

##disp("Matriz A_tri1:");
##disp(A_tri1);
##disp("Matriz b_tri1:");
##disp(b_tri1);


% Resuelvo los sistemas con Sus_atras Ax = b
##x = sust_atras(A_tri1,b_tri1);
##disp("Sistema Resuelto por Gauss");
disp(x);

% x20
disp("x7");
disp(x(7));


% Radios espectrales
[rho_JA] = radio_espectral(A, 'ja', []);
[rho_GS] = radio_espectral(A,'gs', []);


disp("Radio espectral Jacobi:");
disp(rho_JA);
disp("Radio espectral Gauss-Seidel:");
disp(rho_GS);

% Métodos iterativos
[X_GS, it_GS, rh_GS, t_GS] = gauss_seidel(A, b, x0, maxit, tol);
[X_JA, it_JA, rh_JA, t_JA] = jacobi(A, b, x0, maxit, tol);

rja = b - A * X_JA;
residuo_final = norm(b - A * X_JA, inf);
fprintf("Residuo final de Jacobi esperado: %.10e\n", residuo_final);

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

rja = b - A * X_GS;
residuo_final2 = norm(rja, inf);
fprintf("Residuo final de Gauss Seidel esperado: %.10e\n", residuo_final2);
