clc;
format long;

a = -5;
b = 0;
a2 = 0;
b2 = 5;
maxit = 500;
tol = 1e-6;

#Datos
p = @(x) (230 .*x.^4) + (18.*x.^3) + (9 .* x.^2) - (222 .*x) -9;

#Plot funcion

figure(1);
x = linspace(a,b,100);
plot(x,p(x),'b-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

#Plot funcion

figure(2);
x = linspace(a2,b2,100);
plot(x,p(x),'b-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

# Biseccion 1
[ptmin,~,~,~,~] = biseccion3(p, a, b, maxit, tol);
ptmin

%[ptmin,~,~,~] = biseccion(p, a, b, maxit, tol);
%ptmin

# Biseccion 2
[ptmax,~,~,~,~] = biseccion3(p, a2, b2, maxit, tol);
ptmax

[ptmax2,~,~,~] = biseccion(p, a2, b2, maxit, tol);
ptmax2
