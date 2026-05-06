clc;
format short;
% a)
% Funcion
t = 0.02;
z = @(a) 0.04 .* sqrt(a + t) .* (1 - t) - t .* sqrt(3 .* a) - t;

% Datos
a = 19;
b = 20;
maxit = 500;
tol = 1e-6;

# Grafica funcion a)
figure(1);
x = linspace(a,b,100);
plot(x,z(x),'r-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-');
hold on;

[A, ~, ~, ~] = biseccion(z, a, b, maxit, tol);
disp(A);


% b)
z2 = @(t) 0.04 .* sqrt(A + t) .* (1 - t) - t .* sqrt(3 .* A);
a2 = -10;
b2 = 10;

# Grafica funcion b)
figure(2);
x = linspace(a2,b2,100);
plot(x,z2(x),'r-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-');
hold on;

[T, ~, ~, ~] = biseccion(z2, a2, b2, maxit, tol);
disp(T);

