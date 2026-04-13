# Ejercicio 13
N = 100;

# PARTE A
# Matriz del ejercicio 7
A = zeros(N,N);
Bin = [-ones(N,1), 2*ones(N,1), -ones(N,1)]; # Matriz donde cada columna es una diagonal
d = [-1, 0, 1]; # Vector con los offsets de cada diagonal
A = spdiags(Bin, d, N, N);
A(1,:) = 0;  A(1,1) = 1;
A(N,:) = 0;  A(N,N) = 1;
A = full(A);

# Factorización LU (una sola vez, se reutiliza en parte b)
[L, U] = doolittle(A);
disp('=== Parte (a): Factorizacion LU obtenida ===')
disp('L (primeras 5 filas):'); disp(L(1:5, 1:5))
disp('U (primeras 5 filas):'); disp(U(1:5, 1:5))

# PARTE B
disp('=== Parte (b): Maximo de x para cada k ===')
max_x = zeros(10, 1);

for k = 1:10
    # Mismo sistema pero con b = 1/N^k
    b = ones(N,1) / N^k;
    b(1) = 0;
    b(N) = 0;

    # Reutilizamos L y U ya calculadas (ventaja de la factorizacion)
    y = sust_adelante(L, b);
    x = sust_atras(U, y);

    max_x(k) = max(x);
    fprintf('k=%2d  |  1/N^k = %.2e  |  max(x) = %.6e\n', k, 1/N^k, max_x(k));
end

# PARTE C
figure;
semilogy(1:10, max_x, 'o-b', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('k');
ylabel('max(x) - escala logarítmica');
title('Máximo de la solución x en función de k');
grid on;
