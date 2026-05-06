 % Entregable
clear all; clc;

% Inciso 1
% Se define la funcion E(t)-1.5 y su primer derivada
E = @(t) (((t+(1/3)).^3 + (1/3)) .* e.^(-t) - 1.5);
dE = @(t) (e.^(-t) .* ( 3 .* (t+(1/3)).^2 - (t+(1/3)).^3 - (1/3)));

% Grafica de E(t) y E'(t)
 x = [0:0.1:10];
 y_E = E(x);
 y_dE = dE(x);
 z=@(x) 0.*x;
 zz=z(x);
 figure(1);
 plot(x, y_E, x, y_dE,x,zz,'k-','linewidth', 2);
 grid minor
% Maximo de iteraciones y tolerancia
 maxit = 20;
 tol = 1E-12;

 % Calculo de la ra?z en el intervalo [1:2]

 x0_1 = 2;
 [r1, it1, rh1,~]= newton(E, dE, x0_1, maxit, tol);

 % Calculo de la ra?z en el intervalo [3:5]

 x0_2 = 3;
 [r2, it2, rh2,~]= newton(E, dE, x0_2, maxit, tol);

 disp("t -> E = 1.5 | Iteraciones ");
 disp([r1 it1]);
 disp([r2 it2]);

 % Grafica de la convergencia
 figure(2);
 semilogy([1:it1], rh1, [1:it2], rh2);
 grid minor

 % ----------------------------------------------------------------- %

 % Inciso 2
 % Se redefine la funcion de E(t) y se define la segunda derivada
 E = @(t) (((t+(1/3)).^3 + (1/3)) .* e.^(-t));
 d2E = @(t) (e.^(-t) .* ( (t+(1/3)).^3 - 6 .* (t+(1/3)).^2 + 6 .*(t+(1/3)) + (1/3)));

 % Grafica de E(t), E'(t) y E''(t)
 x = [0:0.1:10];
 y_E = E(x);
 y_dE = dE(x);
 y_d2E = d2E(x);

 figure(3);
 plot(x, y_E, x, y_dE, x, y_d2E,x,zz,'k-','linewidth', 2);
 grid minor
 % Calculo del punto cr?tico en el intervalo [2:3]

 x0_3 = 2;
 [r3, it3, rh3,~]= newton(dE, d2E, x0_3, maxit, tol);

 disp("t -> E_max | Iteraciones ");
 disp([r3 it3]);
 disp("Valor de E_max");
 disp(E(r3));

 % Grafica de la convergencia
 figure(4);
 semilogy([1:it3], rh3);
 grid minor
 % ----------------------------------------------------------------- %

 % Inciso 3
 % Se redefine la funcion de E(t) y se define la segunda derivada
 d3E = @(t) (e.^(-t) .* ( -((t+(1/3)).^3) + 9 .* (t+(1/3)).^2 - 18 .*(t+(1/3)) + (17/3)));

 % Grafica de E'(t), E''(t) y E'''(t)
 x = [0:0.1:10];
 y_dE = dE(x);
 y_d2E = d2E(x);
 y_d3E = d3E(x);
 figure(5);
 plot(x, y_dE, x, y_d2E, x, y_d3E,x,zz,'k-','linewidth', 2);
 grid minor
 % Aproximacion al punto cr?tico de E'(t) en el intervalo [0:2]

 xmin = 0;
 xmax = 2;
 maxit_b = 10;
 tol_b = 1E-1;
 [x0_4_aprox, it4_aprox] = biseccion(d2E, xmin, xmax, maxit_b,tol_b);


 [r4, it4, rh4,~]= newton (d2E, d3E, x0_4_aprox, maxit, tol);

 disp("t -> dE_max");
 disp(r4);


