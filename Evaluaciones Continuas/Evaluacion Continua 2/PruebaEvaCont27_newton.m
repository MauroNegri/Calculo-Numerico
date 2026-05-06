clc;

f = @(x) x + exp(-10 .*x.^2) .* cos(x);
df = @(x) -20 .* x .* exp(-10 .*x.^2) .* cos(x) + 1 - exp(-10 .*x.^2) .* sin(x);

x0 = 0;
tol = 1e-6;
a = -2;
b = 2;
maxit = 10000;

% Grafico
figure(1);
x = linspace(a,b,1000);
plot(x, f(x), 'b-');
grid on;
xlabel('x'); ylabel('f(x)');
title('Función');
hold on;
plot(xlim, [0 0], 'k-')  % Línea negra horizontal sobre y = 0
hold on;

% Continuo con Newton
[x_NEW, it_NEW, rh_NEW, t_NEW] = newton(f, df, x0, maxit, tol);

% Mostrar resultados Newton
disp("Método Newton:");
disp(["Raíz: ", num2str(x_NEW, '%.10f')]);
disp(["Iteraciones: ", num2str(it_NEW)]);
disp(["Error final: ", num2str(rh_NEW(end))]);
disp(["Tiempo: ", num2str(t_NEW), " segundos"]);

if it_NEW == maxit || rh_NEW(end) > tol
    disp("Newton no converge en el número máximo de iteraciones.");
end
