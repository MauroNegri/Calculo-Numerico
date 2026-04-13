N = 10;
A = zeros(N,N);

A = diag(-1*ones(N-1,1), 1) + diag(2*ones(N,1),0) + diag(-1*ones(N-1,1), -1);
# Condiciones de borde: x_1 = 0 y x_N = 0
A(1,:) = 0;  A(1,1) = 1;
A(N,:) = 0;  A(N,N) = 1;
disp(A);
# Vector b
b = ones(N,1) / N^2;
b(1) = 0;
b(N) = 0;

# Resolver
x = gauss(A, b)

# Graficar (parte b)
t = linspace(0, 1, N);
plot(t, x);
xlabel('t'); ylabel('x(t)');
title('Solución del sistema');
disp("");

# Con spdiag
Bin = [-ones(N,1), 2*ones(N,1), -ones(N,1)]; # Matriz donde cada columna es una diagonal
d = [-1, 0, 1]; # Vector con los offsets de cada diagonal
A2 = spdiags(Bin, d, N, N);
A2(1,:) = 0;  A2(1,1) = 1;
A2(N,:) = 0;  A2(N,N) = 1;
A2 = full(A2);
# spdiags devuelve una matriz dispersa (sparse)
# full() convierte una matriz sparse en una matriz densa normal
# la funcion gauss fue escrita para operar con matrices densas normales
disp(A2);


