clear;clc;format long;
% Función a minimizar (distancia al cuadrado)
# El punto en y = x^2 más cercano a (1, 0) minimiza
# d^2(x) = (x-1)^2+ (x^2-0)^2= (x-1)^2
f = @(x) (x-1).^2 + x.^4;

% Derivada y segunda derivada
df = @(x) 2*x + 4*x.^3 - 2;
ddf = @(x) 2 + 12*x.^2;

% Datos
a = 0;
b = 1;
x0 = 0.4;  % Punto inicial
tol = 1e-4;
maxit = 2000;

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

% Gráfico de la primera derivada
figure(2);
plot(x, df(x), 'r-', 'LineWidth', 2);
hold on;
plot(xlim, [0 0], 'k--');  % Línea horizontal en y = 0
grid on;
xlabel('x');
ylabel("f'(x)");
title("Primera derivada f'(x)");
legend("f'(x)", "y = 0");

% Aproximo con Bisección
[x_BIS, it_BIS, rh_BIS, t_BIS] = biseccion(df, a, b, maxit, tol);

% Continuo con Newton
[x_NEW, it_NEW, rh_NEW, t_NEW] = newton(df, ddf, x_BIS, maxit, tol);

% Mostrar resultados Newton
disp("Método Newton:");
disp(["Raíz: ", num2str(x_NEW, '%.10f')]);
disp(["Iteraciones: ", num2str(it_NEW)]);
disp(["Error final: ", num2str(rh_NEW(end))]);
disp(["Tiempo: ", num2str(t_NEW), " segundos"]);
