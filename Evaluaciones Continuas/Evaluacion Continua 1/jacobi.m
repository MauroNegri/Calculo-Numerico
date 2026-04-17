function [x, it, r_h, t] = jacobi(A, b, x0, maxit, tol)
##function [x, it, r_h] = jacobi(A, b, x0, maxit, tol)
  tic();                   % Iniciar cronómetro
  n = length(b);           % Tamaño del sistema
  x = x0;                  % Inicializar x
  it = 0;                  % Contador de iteraciones
  r_h = [];                % Vector para guardar los residuos

  while (it < maxit)
    for i = 1:n
      % En Jacobi se usa solo x0 (es decir, la solucion de la iteracion anterior)
      x(i) = (b(i) - A(i, 1:i-1)*x0(1:i-1)- A(i, i+1:n)*x0(i+1:n)) / A(i,i);
    endfor

    % Calcular el error relativo
    r_h(it+1) = norm(x-x0, 'inf') / norm(x, 'inf'); % Norma Infinito
    % r_h(it+1) = norm(x-x0, 2) / norm(x, 2); % Norma 2 (Euclidiana)
    % r_h(it+1) = norm(x-x0, 1) / norm(x, 1); % Norma 1
    %r_h(it+1) = norm(x-x0, 'inf');

    % Calcular absoluto
##     r_h(it+1) = norm(x-x0, 'inf')

    % Calcular Residuo
    % r_h(it+1) = norm(A*x-b, 'inf') / norm(b, 'inf'); % NORMALIZADO Error Relativo del residuo
    % r_h(it+1) = norm(A*x - b, 'inf'); % NO NORMALIZADO Error Absoluto del residuo

    % Verificar convergencia
    if r_h(it+1) < tol
      break
    endif
    x0 = x;        % Actualizar x0 para la proxima iteracion
    it = it + 1;   % Contador de iteraciones
  endwhile

  if it==maxit
    disp('Se ha llegado al Nro maximo de iteraciones')
  endif
  t = toc();       % Finalizar cronómetro y devolver tiempo
endfunction
