function analizar_funcion(f, df, a, b)
  % Analiza características de f en [a,b] y recomienda método
  % f: función anónima, ej: @(x) x.^3 + x - 4
  % df: derivada de f, ej: @(x) 3*x.^2 + 1
  % a, b: intervalo

  fprintf('\n=== ANÁLISIS DE LA FUNCIÓN ===\n\n');

  % 1. Graficar
  x = linspace(a, b, 1000);
  y = f(x);

  figure;
  plot(x, y, 'b-', 'LineWidth', 2);
  hold on;
  plot(0, 'k--', 'LineWidth', 0.5);
  grid on;
  xlabel('x'); ylabel('f(x)');
  title(['f(x) en [', num2str(a), ', ', num2str(b), ']']);

  % 2. Buscar cambios de signo (ceros aproximados)
  cambios = [];
  for i = 1:length(y)-1
    if sign(y(i)) * sign(y(i+1)) < 0
      cambios = [cambios; x(i), x(i+1)];
    endif
  endfor

  n_ceros = size(cambios, 1);
  fprintf('Ceros encontrados: %d\n\n', n_ceros);

  if n_ceros == 0
    fprintf('No se encontraron ceros en el intervalo.\n');
    return;
  endif

  % 3. Analizar cada cero
  for i = 1:n_ceros
    xmin = cambios(i, 1);
    xmax = cambios(i, 2);
    x_aprox = (xmin + xmax) / 2;

    fprintf('--- Cero #%d ---\n', i);
    fprintf('Intervalo: [%.4f, %.4f]\n', xmin, xmax);
    fprintf('Aproximación inicial: %.4f\n', x_aprox);

    % Marcar en el gráfico
    plot(x_aprox, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    text(x_aprox, 0, sprintf('  x₀≈%.2f', x_aprox), 'FontSize', 10);

    % Evaluar derivada
    df_val = abs(df(x_aprox));
    fprintf('|f''(x₀)| ≈ %.4f\n', df_val);

    % Detectar raíz múltiple
    if abs(f(x_aprox)) < 1e-6 && df_val < 0.1
      fprintf('Posible RAÍZ MÚLTIPLE (f≈0 y f''≈0)\n');
      fprintf('Newton pierde convergencia cuadrática\n');
      raiz_multiple = true;
    else
      raiz_multiple = false;
    endif

    % Recomendar método
    fprintf('\nRecomendación:\n');

    if raiz_multiple
      fprintf('Bisección (confiable para raíces múltiples)\n');
      fprintf('Newton modificado: xₙ₊₁ = xₙ - m·f(xₙ)/f''(xₙ)\n');
    elseif df_val < 0.5
      fprintf('Derivada pequeña → curva casi horizontal\n');
      fprintf('Bisección (más confiable)\n');
      fprintf('Criterio de corte: |pₙ - pₙ₋₁| < tol\n');
    elseif df_val > 5
      fprintf('Derivada grande → curva casi vertical\n');
      fprintf('Newton-Raphson (rápido con buen x₀)\n');
      fprintf('Criterio de corte: |f(pₙ)| < tol\n');
    else
      fprintf('Newton-Raphson (convergencia rápida)\n');
      fprintf('Secante (si no hay derivada)\n');
      fprintf('Bisección como respaldo\n');
    endif

    % Evaluar ancho del intervalo
    ancho = xmax - xmin;
    if ancho > 1
      fprintf('Intervalo amplio → usar bisección para acotar primero\n');
    endif

    fprintf('\n');
  endfor

  % 4. Recomendación general
  fprintf('=== ESTRATEGIA GENERAL ===\n');
  if n_ceros > 1
    fprintf('️Múltiples ceros detectados\n');
    fprintf('Definir intervalos que encierren solo un cero\n');
    fprintf('Resolver cada uno por separado\n\n');
  endif

  hold off;
endfunction
