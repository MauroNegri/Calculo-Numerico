function [A, r] = doolittle_p(A)
  n = length(A);
  r = 1:n;
  epsilon = 1e-9;
  for k = 1:n-1
    % Buscar pivote mayor en la columna k
    [pmax, p] = max(abs(A(r(k:n), k)));

    if pmax < epsilon
      disp('Los posibles pivots son CERO')
      break
    endif

    p = p + k - 1;              % posición global
    if p ~= k
      r([p k]) = r([k p]);      % actualizar vector de permutaciones
    endif

    % Almacenar multiplicadores (como doolittle)
    A(r(k+1:n), k) = A(r(k+1:n), k) / A(r(k), k);

    % Operaciones elementales
    A(r(k+1:n), k+1:n) = A(r(k+1:n), k+1:n) - A(r(k+1:n), k) * A(r(k), k+1:n);
  endfor
endfunction

# Para extraer L y U en el main
# L = tril(A(r,:), -1) + eye(n);
# U = triu(A(r,:));
