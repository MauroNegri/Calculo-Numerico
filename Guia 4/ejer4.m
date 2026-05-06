clc;
format long;
f = @(x) sin(x) + cos(1+x.^2)-1;
df = @(x) cos(x) - 2*x .* sin(1 + x.^2);  % Derivada para Newton
tol = 1e-10;
maxit = 100;
a1 = 0;
b1 = 2;
a2 = 2.1;
b2 = 3;

% Grafico
figure(1);
x = linspace(a1,b1,1000);
plot(x, f(x), 'b-');
grid on;
xlabel('x'); ylabel('f(x)');
title('Función');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

% Métodos iterativos
[x_BIS, it_BIS, rh_BIS, t_BIS] = biseccion(f, a1, b1, maxit, tol);
[x_BIS2, it_BIS2, rh_BIS2, t_BIS2] = biseccion(f, a2, b2, maxit, tol);
[x_NEW, it_NEW, rh_NEW, t_NEW] = newton(f, df, x_BIS, maxit, tol);
[x_NEW2, it_NEW2, rh_NEW2, t_NEW2] = newton(f, df, x_BIS2, maxit, tol);
[x_SEC, it_SEC, rh_SEC, t_SEC] = secante(f, a1, b1, maxit, tol);
[x_SEC2, it_SEC2, rh_SEC2, t_SEC2] = secante(f, a2, b2, maxit, tol);

% Mostrar resultados Bisección
disp("Método Bisección Punto 1:");
disp(["Raíz 1: ", num2str(x_BIS, '%.10f')]);
disp(["Iteraciones: ", num2str(it_BIS)]);
disp(["Error final: ", num2str(rh_BIS(end))]);
disp(["Tiempo: ", num2str(t_BIS), " segundos"]);
disp("Método Bisección Punto 2:");
disp(["Raíz 2: ", num2str(x_BIS2, '%.10f')]);
disp(["Iteraciones: ", num2str(it_BIS2)]);
disp(["Error final: ", num2str(rh_BIS2(end))]);
disp(["Tiempo: ", num2str(t_BIS2), " segundos"]);

disp("-----------------------------------------------------------------");
% Mostrar resultados Newton
disp("Método Newton Punto 1:");
disp(["Raíz 1: ", num2str(x_NEW, '%.10f')]);
disp(["Iteraciones: ", num2str(it_NEW)]);
disp(["Error final: ", num2str(rh_NEW(end))]);
disp(["Tiempo: ", num2str(t_NEW), " segundos"]);
disp("Método Newton Punto 2:");
disp(["Raíz 2: ", num2str(x_NEW2, '%.10f')]);
disp(["Iteraciones: ", num2str(it_NEW2)]);
disp(["Error final: ", num2str(rh_NEW2(end))]);
disp(["Tiempo: ", num2str(t_NEW2), " segundos"]);

disp("-----------------------------------------------------------------");
% Mostrar resultados Secante
disp("Método Secante Punto 1:");
disp(["Raíz 1: ", num2str(x_SEC, '%.10f')]);
disp(["Iteraciones: ", num2str(it_SEC)]);
disp(["Error final: ", num2str(rh_SEC(end))]);
disp(["Tiempo: ", num2str(t_SEC), " segundos"]);
disp("Método Secante Punto 2:");
disp(["Raíz 2: ", num2str(x_SEC2, '%.10f')]);
disp(["Iteraciones: ", num2str(it_SEC2)]);
disp(["Error final: ", num2str(rh_SEC2(end))]);
disp(["Tiempo: ", num2str(t_SEC2), " segundos"]);


disp("-----------------------------------------------------------------");
% Utilice el metodo de newton partiendo de x0 = 1 y seleccione la respuesta correcta:
% a. Converge a la cuarta raíz positiva.
% b. Converge a la tercera raíz positiva.
% c. Converge a la raíz postiva más cercana a 1.
% d. Converge a la quinta raíz postivia.
% e. Converge a la sexta raíz postivia.
% f. La iteración diverge.

[x_NEW2, it_NEW2, rh_NEW2, t_NEW2] = newton(f, df, 1, maxit, tol);

% Mostrar resultados Newton
disp("Método Newton:");
disp(["Raíz 1: ", num2str(x_NEW2, '%.10f')]);
disp(["Iteraciones: ", num2str(it_NEW2)]);
disp(["Error final: ", num2str(rh_NEW2(end))]);
disp(["Tiempo: ", num2str(t_NEW2), " segundos"]);

% Para ver eso grafico hasta x_NEW2 y cuento cuantas raices hay
% Converge a la cantidad de raices + x_NEW2
figure(2);
x = linspace(0, x_NEW2, 1000);
plot(x, f(x), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('x'); ylabel('f(x)');
title('Raíces de f(x) = sin(x) + cos(1+x^2) - 1');
hold on;
plot(xlim, [0 0], 'k--');

% d. Converge a la quinta raíz postivia.
