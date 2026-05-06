format long;

p = @(x) x - x.^3 - (4 .* x.^2) +10;
a1 = -4;
b1 = -3;
a2 = 1;
b2 = 2;
maxit = 500;
tol = 1e-7;

#Plot funcion raiz1

figure(1);
x = linspace(a1,b1,100);
plot(x,p(x),'b-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

#Plot funcion raiz2

figure(2);
x = linspace(a2,b2,100);
plot(x,p(x),'b-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

% a) Dos raices, cual es la raiz mas chica y cual la mas grande
# Biseccion
[raizchica,~,~,~] = biseccion(p, a1, b1, maxit, tol);
[raizgrande,~,~,~] = biseccion(p, a2, b2, maxit, tol);
raizchica
raizgrande

#Plot funcion

figure(3);
x = linspace(-3,3,100);
plot(x,p(x),'b-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

% b) El polinomio alcanza un maximo relativo en x = ? y ese valor es ?
dp = @(x) (-3.*x.^2) - (8.*x) +1;

#Plot funcion derivada
figure(4);
x = linspace(-1,1,100);
plot(x,dp(x),'b-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

[raizDerivada,~,~,~] = biseccion(dp, -1, 1, maxit, tol);
valorMaximo = p(raizDerivada);
valorMaximo

% c) El polinomio posee un punto fijo en x = ?
# Punto Fijo
% ppf = @(x) x - x.^3 - (4 .* x.^2) +10 - x; Se restan las x
ppf = @(x) - x.^3 - (4 .* x.^2) +10;

figure(5);
x = linspace(1,2,100);
plot(x,ppf(x),'b-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

[puntofijo,~,~,~] = biseccion(ppf, 1,2,maxit, tol);
puntofijo


