clc;
% Funcion original
%z = @(a) 0.04 .* sqrt(a+0.02).*(1-0.02)-0.02 .*sqrt(3 .*a);

% Funcion con el punto fijo t = 0.02
z1 = @(a) 0.04 .* sqrt(a + 0.02) .* (1 - 0.02) - 0.02 .* sqrt(3 .* a) - 0.02;
aa = 19;
b = 20;
maxit = 500;
tol = 1e-6;

% a) Encontrar el valor positivo del parametro a tal que posea un punto fijo en t = 0.02
[p,~, it, ~] = biseccion(z1, aa, b, maxit, tol);
p % Parametro a

% b) Encuentre la raiz de z(t) con un error de 10^-6
z2 = @(t) 0.04.*sqrt(p+t).*(1-t)-t.*sqrt(3.*p);

% Grafico para encontrar donde corta al eje x
figure(1);
a = linspace(-5,5,100);
plot(a,z2(a),'b-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

[p2,~, it, ~] = biseccion(z2, -5, 5, maxit, tol);
p2 % Raiz de z(t)
