clc;
format long;

# Datos
g = 9.8;
c = 15;
t = 7;
v = 35;
a = 68;
b = 70;
maxit = 500;
tol = 1e-7; % 5 cifras decimales exactas
% a) Determinar m
f = @(m) (g.*m).*(1-exp(-(c./m).*t)) - (c.*v);

#Plot funcion

figure(1);
x = linspace(a,b,100);
plot(x,f(x),'b-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

# Biseccion
[masa,~,~,~] = biseccion(f, a, b, maxit, tol);
masa

% b) Si la masa es m = 73, el coeficiente de arrastre deberia ser c = ?
m2 = 73;
f2 = @(c) (g.*m2).*(1-exp(-(c./m2).*t)) - (c.*v);
a2 = 15;
b2 = 17;
#Plot funcion

figure(2);
x = linspace(a2,b2,100);
plot(x,f2(x),'b-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

# Biseccion
[c,~,~,~] = biseccion(f2, a2, b2, maxit, tol);
c
