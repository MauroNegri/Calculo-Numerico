#function [x, it, r_h, t] = sor(A, b, x0, maxit, tol, w)
function [x, it, r_h] = sor(A, b, x0, maxit, tol, w)
  #tic();                   % Iniciar cronómetro
  n = length(b);           % Número de incógnitas (puede ser size(A,1) también)
  x = x0;                  % Inicializar x
  it = 0;                  % Contador de iteraciones
  r_h = [];                % Vector para guardar el residuo (error relativo)

  while (it < maxit)
    for i = 1:n
      % En SOR se mezcla el método de Gauss-Seidel con un factor de relajación w
      % Usa los valores actualizados de x(1:i-1) y los valores viejos de x0(i+1:n)
      x(i) = (1-w) * x0(i) + w *( b(i) - A(i,1:i-1)*x(1:i-1) ...
      - A(i,i+1:n)*x0(i+1:n) ) / A(i,i);
    endfor

    % Calcular el residuo (puede ser también el error relativo como en GS/Jacobi)
    % r_h(it + 1) = norm(x-x0, 'inf') / norm(x, 'inf'); % Norma Infinito
    % r_h(it+1) = norm(x-x0, 2) / norm(x, 2); % Norma 2 (Euclidiana)
    % r_h(it+1) = norm(x-x0, 1) / norm(x, 1); % Norma 1
    %it = it + 1;

        % Calcular error absoluto
     r_h(it+1) = norm(x-x0, 'inf');

    % Calcular Residuo
    % r_h(it+1) = norm(A*x-b, 'inf') / norm(b, 'inf'); % NORMALIZADO Error Relativo del residuo
    % r_h(it+1) = norm(A*x - b, 'inf'); % NO NORMALIZADO Error Absoluto del residuo

    % Verificar convergencia
    if r_h(it+1) < tol
      break;
    endif
    x0 = x;        % Actualizar x0 para la próxima iteración
    it = it + 1;   % Contador de iteraciones
  endwhile
  
  if it==maxit
    disp('Se ha llegado al Nro maximo de iteraciones')
  endif
  #t = toc();       % Finalizar cronómetro y devolver tiempo
 endfunction
