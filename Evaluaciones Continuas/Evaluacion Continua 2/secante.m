function [p, it, r, t] = secante(f, p0, p1, maxit, tol)
  tic();  % comienza el cronómetro

  % Evaluaciones iniciales
  q0 = f(p0);
  q1 = f(p1);

  r = zeros(1, maxit);  % prealocamos vector de errores
  it = 2;  % ya tenemos dos puntos iniciales

  while it < maxit
    % Fórmula de la secante
     if (q1 - q0) == 0
      error("Division por cero en secante.");
    endif
    p = p1 - q1 * (p1 - p0) / (q1 - q0);

    % usamos una condicion de corte
    %   abs(f(p)) < tolerancia
    %   abs(p-p1) < tolerancia
    %   abs(p-p1)/abs(p) < tolerancia

    ## Criterio de corte 1: |f(p)| < tol
    ## r(it) = abs(f(p));
    ## if r(it) < tol
    ##   break;
    ## endif

    ## Criterio de corte 2: |p - p1| < tol
    ## r(it) = abs(p - p1);
    ## if r(it) < tol
    ##   break;
    ## endif

    ## Criterio de corte 3: |p - p1| / |p| < tol
    r(it) = abs(p - p1) / abs(p);
    if r(it) < tol
      break;
    endif


    % Actualizar variables para la próxima iteración
    p0 = p1;
    q0 = q1;
    p1 = p;
    q1 = f(p);

    it = it + 1;
  endwhile

  if it == maxit
    disp("No se encontró raíz en maxit iteraciones");
  endif

  t = toc();  % tiempo total
endfunction
