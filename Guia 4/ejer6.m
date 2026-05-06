format long;

% Datos
A = 0.401;
B = 42.7e-6;
T = 300;
P = 3.5e+7;
k = 1.3806503e-23;
N = 1000;
tol = 1e-12;
maxit = 100;
a = 0;
b = 1;

% Función
f = @(V) (P + A * (N ./ V).^2) .* (V - N * B) - (k * N * T);
df = @(V) -2 * A * N^2 ./ V.^3 .* (V - N * B) + P + A .* (N ./ V).^2;

% Grafico
figure(1);
x = linspace(a,b,100);
plot(x, f(x), 'b-');
grid on;
xlabel('x'); ylabel('f(x)');
title('Función');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

% Aproximo con Bisección
[x_BIS, it_BIS, rh_BIS, t_BIS] = biseccion(f, a, b, maxit, tol);

% Continuo con Newton
[x_NEW, it_NEW, rh_NEW, t_NEW] = newton(f, df, x_BIS, maxit, tol);

% Mostrar resultados Bisección
disp("Método Bisección:");
disp(["Raíz: ", num2str(x_BIS, '%.10f')]);
disp(["Iteraciones: ", num2str(it_BIS)]);
disp(["Error final: ", num2str(rh_BIS(end))]);
disp(["Tiempo: ", num2str(t_BIS), " segundos"]);

% Mostrar resultados Newton
disp("Método Newton:");
disp(["Raíz: ", num2str(x_NEW, '%.10f')]);
disp(["Iteraciones: ", num2str(it_NEW)]);
disp(["Error final: ", num2str(rh_NEW(end))]);
disp(["Tiempo: ", num2str(t_NEW), " segundos"]);




