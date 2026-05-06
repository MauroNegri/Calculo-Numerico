format long;
%t = 1.5;
E = @(t) ((t + 1/3).^3 + 1/3).* exp(-t) -1.5;
%dE = @(t) (54.*t.^2+45.*t -1-27.*t.^3-27.*exp(t)) ./ (27 .* exp(t));
dE = @(t) exp(-t) .* (3*(t + 1/3).^2 - ((t + 1/3).^3 + 1/3));

a1 = 0;
b1 = 10;
%t0 = 0.9;
tol = 1e-12;
maxit = 100;

# Grafica funcion a)
figure(1);
x = linspace(a1,b1,100);
plot(x,E(x),'r-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-');
hold on;

% Dos intervalos con cambio de signo
a1 = 0;   b1 = 1.6;    % primera raiz
a2 = 3;   b2 = 4;    % segunda raiz

% Biseccion + Newton para cada una
[x_BIS1, it_BIS1, rh_BIS1, ~] = biseccion(E, a1, b1, maxit, tol);
[x_NEW1, it_NEW1, rh_NEW1, ~] = newton(E, dE, x_BIS1, maxit, tol);

[x_BIS2, it_BIS2, rh_BIS2, ~] = biseccion(E, a2, b2, maxit, tol);
[x_NEW2, it_NEW2, rh_NEW2, ~] = newton(E, dE, x_BIS2, maxit, tol);

fprintf("Raiz 1: t = %.10f,  E(t) = %.6e\n", x_NEW1, E(x_NEW1)+1.5);
fprintf("Raiz 2: t = %.10f,  E(t) = %.6e\n", x_NEW2, E(x_NEW2)+1.5);

##% Mostrar resultados Newton
##disp("Método Newton a):");
##disp(["Raíz: ", num2str(x_NEW, '%.10f')]);
##disp(["Iteraciones: ", num2str(it_NEW)]);
##disp(["Error final: ", num2str(rh_NEW(end))]);
##disp(["Tiempo: ", num2str(t_NEW), " segundos"]);
##fprintf("La energía es 1.5 en t = %.5f\n", x_NEW);

% 2) a dE/dt la evaluo en 0 ahi obtengo el maximo
% ahora aplico newton de vuelta y derivo de nuevo

%ddE = @(t) (27.*t.^3 - 135.*t.^2 + 63.*t + 46) ./ (27 .* exp(t));
##ddE = @(t) exp(-t) .* ( ...
##    - (3*(t + 1/3).^2 - ((t + 1/3).^3 + 1/3)) ...
##    + 6*(t + 1/3) - 3*(t + 1/3).^2 );
ddE = @(t) (e.^(-t) .* ( (t+(1/3)).^3 - 6 .* (t+(1/3)).^2 + 6 .*(t+(1/3)) + (1/3)));

a2 = 0;
b2 = 3;

# Grafica funcion b)
figure(2);
x = linspace(a2,b2,100);
plot(x,dE(x),'r-');
grid on;
title('Funcion');
hold on;
plot(xlim, [0 0], 'k-');
hold on;

% Aproximo con Bisección
%[x_BIS2, it_BIS2, rh_BIS2, t_BIS2] = biseccion(dE, a, b, maxit, tol);
t0 = 2;
% Newton
[x_NEW2, it_NEW2, rh_NEW2, t_NEW2] = newton(dE, ddE, t0, maxit, tol);

disp("-----------------------------------------------------------------");
% Mostrar resultados Newton
disp("Método Newton b):");
disp(["Raíz: ", num2str(x_NEW2, '%.10f')]);
disp(["Iteraciones: ", num2str(it_NEW2)]);
disp(["Error final: ", num2str(rh_NEW2(end))]);
disp(["Tiempo: ", num2str(t_NEW2), " segundos"]);


E_max = E(x_NEW2);
fprintf("Máxima energía = %.10f en t = %.10f\n", E_max, x_NEW2);

% -------------------------------------------------------------------------
% INCISO (c): Encontrar el t donde E'(t) es máxima => E''(t) = 0
% -------------------------------------------------------------------------

% Derivada tercera para usar en Newton
% dddE = derivada de ddE
dddE = @(t) exp(-t) .* ( ...
    -6*(t + 1/3) + 3*(t + 1/3).^2 - (6*(t + 1/3) - 3*(t + 1/3).^2) ...
    - (3*(t + 1/3).^2 - ((t + 1/3).^3 + 1/3)) ...
);

% Gráfica de E''(t) (segunda derivada)
figure(3);
x = linspace(a2, b2, 100);
plot(x, ddE(x), 'b-');
grid on;
title("Derivada segunda E''(t)");
hold on;
plot(xlim, [0 0], 'k--');

% Newton para encontrar raíz de E''(t)
t0_c = 1.0;  % valor inicial razonable para Newton
[x_NEW3, it_NEW3, rh_NEW3, t_NEW3] = newton(ddE, dddE, t0_c, maxit, tol);

% Evaluar E'(t) en el punto encontrado
valor_dE_max = dE(x_NEW3);

% Mostrar resultados
disp("-----------------------------------------------------------------");
disp("Método Newton c):");
disp(["t donde E''(t) = 0: ", num2str(x_NEW3, '%.10f')]);
disp(["Iteraciones: ", num2str(it_NEW3)]);
disp(["Error final: ", num2str(rh_NEW3(end))]);
disp(["Tiempo: ", num2str(t_NEW3), " segundos"]);
fprintf("Máximo valor de E'(t) = %.10f en t = %.10f\n", valor_dE_max, x_NEW3);
