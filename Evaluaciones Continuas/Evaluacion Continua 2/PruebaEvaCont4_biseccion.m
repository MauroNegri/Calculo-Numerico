format long;
# Datos
pt = 3.5;

k = 0.04;
maxit = 500;
tol = 1e-8;
%Cifras decimales exactas despues de la coma entonces va los 8,
% Sino pide cifras significativas entonces es uno menos.
a = -1;
b = 1;

H = @(x) (x ./  (1-x)) .* sqrt((2 .* pt) ./ (2+x)) -k;

#Plot
figure(1);
x = linspace(a,b,100);
plot(x,H(x),'b-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

# Biseccion
[p,~,~,~] = biseccion(H, a, b, maxit, tol);
p

H2 = @(pt) (0.02 ./  (1-0.02)) .* sqrt((2 .* pt) ./ (2+0.02)) -k;

a2 = 0;
b2 = 4;

#Plot
figure(2);
x = linspace(a2,b2,100);
plot(x,H2(x),'b-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

# Biseccion
[p2,~,~,~] = biseccion(H2, a2, b2, maxit, tol);
p2
