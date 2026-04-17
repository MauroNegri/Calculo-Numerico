n = 40;
B = [ -1*ones(n,1), 2*ones(n,1), -1*ones(n,1)];
d = [-1, 0, 1];

A = spdiags(B,d,n,n);
b = 1.5 * (1:n)' - 6;
x0 = zeros(n,1);
w = 1.85;
tol = 1e-5;
maxit = 2000;

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

% Valor de x20
disp("Valor de x20");
disp(x(20));

% Métodos iterativos
[X_JA, it_JA, rh_JA, t_JA] = jacobi(A, b, x0, maxit, tol);
[X_GS, it_GS, rh_GS, t_GS] = gauss_seidel(A, b, x0, maxit, tol);
[X_SOR, it_SOR, rh_SOR, t_SOR] = sor(A,b, x0, maxit, tol, w);

% Mostrar resultados Jacobi
disp("Método Jacobi:");
disp(["Iteraciones: ", num2str(it_JA)]);

% Mostrar resultados Gauss-Seidel
disp("Método Gauss-Seidel:");
disp(["Iteraciones: ", num2str(it_GS)]);

% Mostrar resultados SOR
disp("Método SOR:");
disp(["Iteraciones: ", num2str(it_SOR)]);

