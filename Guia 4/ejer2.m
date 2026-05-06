% Definir funciones
clear;clc;
g1 = @(x) 0.5* sqrt(10 - x.^3);
g2 = @(x) sqrt(10 ./ (4+x));

% Definir parámetros
p0 = 1.5;
tol = 1e-3;
maxit = 100; % máximo de iteraciones

% Aplicar método de punto fijo a g1 y g2
[p_g1, it_g1, r_g1, t_g1] = puntofijo(g1,p0, maxit, tol);
[p_g2, it_g2, r_g2, t_g2] = puntofijo(g2,p0, maxit, tol);

% Mostrar resultados g1
disp('Resultados para g1(x):');
disp(['p_g1 = ', num2str(p_g1)]);
disp(['Iteraciones = ', num2str(it_g1)]);
disp(['Tiempo = ', num2str(t_g1), ' segundos']);
disp('Residuos:');
disp(r_g1);
% Mostrar resultados g2
disp('Resultados para g2(x):');
disp(['p_g2 = ', num2str(p_g2)]);
disp(['Iteraciones = ', num2str(it_g2)]);
disp(['Tiempo = ', num2str(t_g2), ' segundos']);
disp('Residuos:');
disp(r_g2);

% Verifico que f(p) = 0 (aprox)
% f(x)
f = @(x) x.^3 + 4*x.^2 - 10;

% Evaluar f en los puntos obtenidos
fp_g1 = f(p_g1);
fp_g2 = f(p_g2);

% Resultados
disp('Chequeo de la raíz para g1:')
disp(['f(p_g1) = ', num2str(fp_g1)])

disp('Chequeo de la raíz para g2:')
disp(['f(p_g2) = ', num2str(fp_g2)])

% --- Gráfica de residuos para g1 y g2 ---

figure;
semilogy(1:it_g1, abs(r_g1), 'r-o', 'DisplayName', 'g1(x)');
hold on;
semilogy(1:it_g2, abs(r_g2), 'b-s', 'DisplayName', 'g2(x)');
xlabel('Iteraciones');
ylabel('Residuo |p_{n+1} - p_n|');
title('Convergencia de método de Punto Fijo');
legend;
grid on;

% --- Comparación de los puntos p_g1 y p_g2 ---
diferencia = abs(p_g1 - p_g2);
disp(['Diferencia entre p_g1 y p_g2: ', num2str(diferencia)]);

% Derivadas de g1 y g2
dg1 = @(x) -0.75 * x.^2 ./ sqrt(10 - x.^3);
dg2 = @(x) -0.5 * sqrt(10) ./ (4 + x).^(3/2);

% Evaluar |g'(x)| en [1, 2]
x_test = linspace(1, 2, 100);
k1 = max(abs(dg1(x_test)));
k2 = max(abs(dg2(x_test)));

fprintf('\n=== Análisis de convergencia (Teorema 2.3) ===\n');
fprintf('Para g1: max|g1''(x)| en [1,2] = %.4f\n', k1);
if k1 < 1
    fprintf('k1 < 1 → g1 CONVERGE\n');
else
    fprintf('k1 >= 1 → g1 NO converge o converge lento\n');
endif

fprintf('\nPara g2: max|g2''(x)| en [1,2] = %.4f\n', k2);
if k2 < 1
    fprintf('k2 < 1 → g2 CONVERGE\n');
else
    fprintf('k2 >= 1 → g2 NO converge o converge lento\n');
endif

% Comparar velocidad de convergencia
% Comparar velocidad de convergencia
fprintf('\nVelocidad de convergencia:\n');
if k1 < k2
    fprintf('g1: k = %.4f → MÁS RÁPIDA\n', k1);
else
    fprintf('g1: k = %.4f → más lenta\n', k1);
endif

if k2 < k1
    fprintf('g2: k = %.4f → MÁS RÁPIDA\n', k2);
else
    fprintf('g2: k = %.4f → más lenta\n', k2);
endif

fprintf('\n=== Corolario 2.4: Cotas de error ===\n');

% Para g1
if k1 < 1
    % Cota |pn - p| <= k^n * max(p0 - a, b - p0)
    cota1_a = k1^it_g1 * max(p0 - 1, 2 - p0);
    % Cota |pn - p| <= (k^n / (1-k)) * |p1 - p0|
    cota1_b = (k1^it_g1 / (1 - k1)) * abs(g1(p0) - p0);

    fprintf('g1:\n');
    fprintf('  Cota (a): |p%d - p| <= %.6e\n', it_g1, cota1_a);
    fprintf('  Cota (b): |p%d - p| <= %.6e\n', it_g1, cota1_b);
    fprintf('  Error real: %.6e\n', abs(fp_g1));
endif

% Para g2
if k2 < 1
    cota2_a = k2^it_g2 * max(p0 - 1, 2 - p0);
    cota2_b = (k2^it_g2 / (1 - k2)) * abs(g2(p0) - p0);

    fprintf('g2:\n');
    fprintf('  Cota (a): |p%d - p| <= %.6e\n', it_g2, cota2_a);
    fprintf('  Cota (b): |p%d - p| <= %.6e\n', it_g2, cota2_b);
    fprintf('  Error real: %.6e\n', abs(fp_g2));
endif
