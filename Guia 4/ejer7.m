format long;

f = @(x) x .* (log(x+3) -17) -1;
df = @(x) log(x+3) + x ./ (x + 3) - 17;

maxit = 50;
tol1 = 1e+2;
tol2 = 1e-12;
a = 0;
b = 1;

while f(a)*f(b) > 0
  b *= 2;
endwhile
disp(b);
## Para bisección se puede calcular exactamente cuántas iteraciones necesita antes de correr:
n_cota = ceil(log2((b-a)/tol1));
fprintf("Cota de iteraciones necesarias: %d\n", n_cota);

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
[x_BIS, it_BIS, rh_BIS, t_BIS] = biseccion(f, a, b, maxit, tol1);

% Continuo con Newton
[x_NEW, it_NEW, rh_NEW, t_NEW] = newton(f, df, x_BIS, maxit, tol2);

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

fprintf("f(x_BIS) = %.6e\n", f(x_BIS));
fprintf("f(x_NEW) = %.6e\n", f(x_NEW));

figure(2);
x = linspace(x_NEW * 0.9999, x_NEW * 1.0001, 1000);
plot(x, f(x), 'b-', 'LineWidth', 1.5);
hold on;
plot(xlim, [0 0], 'k--');
grid on;
xlabel('x'); ylabel('f(x)');
title('Zoom cerca de la raíz');

fprintf("\nVerificacion:\n");
fprintf("Raiz Newton:       %.10f\n", x_NEW);
fprintf("f(raiz):           %.6e  (debe ser ~0)\n", f(x_NEW));
fprintf("Error relativo:    %.6e  (debe ser < %.0e)\n", rh_NEW(end), tol2);
fprintf("Intervalo inicial: [%.0f, %.0f]\n", a, b);
fprintf("b encontrado por duplicacion: %.0f = 2^%.0f\n", b, log2(b));
