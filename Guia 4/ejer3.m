% Definir funciones
clear;clc;
g1 = @(x) (3 +x -2*x.^2).^(1/4);
g2 = @(x) ((x + 3 - x.^4) ./ 2).^(1/2);
g3 = @(x) ((x + 3) ./ (x.^2 + 2)).^(1/2);
g4 = @(x) (3*x.^4 +2*x.^2 +3) ./ (4*x.^3 + 4*x -1);

% Definir parámetros
p0 = 1;
tol = 1e-3;
maxit = 50; % máximo de iteraciones

% Aplicar método de punto fijo
[p_g1, it_g1, r_g1, t_g1] = puntofijo(g1,p0,maxit,tol);
[p_g2, it_g2, r_g2, t_g2] = puntofijo(g2,p0,maxit,tol);
[p_g3, it_g3, r_g3, t_g3] = puntofijo(g3,p0,maxit,tol);
[p_g4, it_g4, r_g4, t_g4] = puntofijo(g4,p0,maxit,tol);

% Mostrar resultados g1
disp('Resultados para g1(x):');
disp(['p_g1 = ', num2str(p_g1)]);
disp(['Iteraciones = ', num2str(it_g1)]);
disp(['Tiempo = ', num2str(t_g1), ' segundos']);
disp('Residuos:');
%disp(r_g1);

% Mostrar resultados g2
disp('Resultados para g2(x):');
disp(['p_g2 = ', num2str(p_g2)]);
disp(['Iteraciones = ', num2str(it_g2)]);
disp(['Tiempo = ', num2str(t_g2), ' segundos']);
disp('Residuos:');
%disp(r_g2);

% Mostrar resultados g3
disp('Resultados para g3(x):');
disp(['p_g3 = ', num2str(p_g3)]);
disp(['Iteraciones = ', num2str(it_g3)]);
disp(['Tiempo = ', num2str(t_g3), ' segundos']);
disp('Residuos:');
%disp(r_g3);

% Mostrar resultados g4
disp('Resultados para g4(x):');
disp(['p_g4 = ', num2str(p_g4)]);
disp(['Iteraciones = ', num2str(it_g4)]);
disp(['Tiempo = ', num2str(t_g4), ' segundos']);
disp('Residuos:');
%disp(r_g4);

% Verifico que f(p) = 0 (aprox)
% f(x)
f = @(x) x.^4 + 2*x.^2 - x - 3;

% Evaluar f en los puntos obtenidos
fp_g1 = f(p_g1);
fp_g2 = f(p_g2);
fp_g3 = f(p_g3);
fp_g4 = f(p_g4);

% Resultados
disp('Chequeo de la raíz para g1:')
disp(['f(p_g1) = ', num2str(fp_g1)])

disp('Chequeo de la raíz para g2:')
disp(['f(p_g2) = ', num2str(fp_g2)])

disp('Chequeo de la raíz para g3:')
disp(['f(p_g3) = ', num2str(fp_g3)])

disp('Chequeo de la raíz para g4:')
disp(['f(p_g4) = ', num2str(fp_g4)])

% Derivadas
dg1 = @(x) (1/4).*(1 - 4*x).*(3 + x - 2*x.^2).^(-3/4);
dg2 = @(x) (1/4).*(1 - 4*x.^3)./sqrt((x + 3 - x.^4)/2);
dg3 = @(x) (1/2).*((x.^2+2).*(1)-(x+3).*(2*x))./((x.^2+2).^2 * sqrt((x+3)/(x.^2+2)));
dg4 = @(x) ((12*x.^3 + 4*x).*(4*x.^3 + 4*x -1) - (3*x.^4 + 2*x.^2 + 3).*(12*x.^2 + 4)) ./ (4*x.^3 + 4*x -1).^2;

x = linspace(1,1.5,500);

figure;
plot(x, abs(dg1(x)), 'r', 'DisplayName', '|g1''(x)|'); hold on;
plot(x, abs(dg2(x)), 'b', 'DisplayName', '|g2''(x)|');
plot(x, abs(dg3(x)), 'g', 'DisplayName', '|g3''(x)|');
plot(x, abs(dg4(x)), 'm', 'DisplayName', '|g4''(x)|');
plot(1,'k--', 'DisplayName', 'y=1');
xlabel('x');
ylabel('|g''(x)|');
legend();
title('|g''(x)| en [1, 3/2]');
grid on;

% Evaluar |g'(x)| en [1, 1.5]
x_test = linspace(1, 1.5, 100);
k1 = max(abs(dg1(x_test)));
k2 = max(abs(dg2(x_test)));
k3 = max(abs(dg3(x_test)));
k4 = max(abs(dg4(x_test)));

fprintf('\n=== Teorema 2.3: Convergencia ===\n');

fprintf('g1: max|g1''(x)| = %.4f → ', k1);
if k1 < 1
    fprintf('CONVERGE\n');
else
    fprintf('NO converge\n');
endif

fprintf('g2: max|g2''(x)| = %.4f → ', k2);
if k2 < 1
    fprintf('CONVERGE\n');
else
    fprintf('NO converge\n');
endif

fprintf('g3: max|g3''(x)| = %.4f → ', k3);
if k3 < 1
    fprintf('CONVERGE\n');
else
    fprintf('NO converge\n');
endif

fprintf('g4: max|g4''(x)| = %.4f → ', k4);
if k4 < 1
    fprintf('CONVERGE\n');
else
    fprintf('NO converge\n');
endif

% Encontrar la más rápida
k_vals = [k1, k2, k3, k4];
[k_min, idx_min] = min(k_vals);
fprintf('\n→ g%d converge MÁS RÁPIDO (k = %.4f)\n', idx_min, k_min);
