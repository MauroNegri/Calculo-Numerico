function [p,it,h,t] = puntofijo(f, p0, maxit, tol)
##function [x, h] = puntofijo(f, p0, maxit, tol)
  % puntofijo(f,g,p0,tol,maxit) si pide |f(pn)| <= tol, g la funcion original
  tic();

  it = 1;
  while (it <= maxit)
    p = f(p0);

##    r(it) = p-p0;
    % Criterio 1: |g(p)| < tol (si se provee f_original)
##    h(it) = abs(g(p));
##    if h(it) < tol
##      break;
##    endif

    % Criterio de corte 2: |pn - pn-1| <= tol Error absoluto
    h(it) = abs(p - p0);
    if h(it) < tol
      break;
    endif

    % Criterio de corte 3: |pn - pn-1| / |pn| <= tol Error relativo
    %if (abs(p - p0) / abs(p) < tol)
     % break;
    %endif
##    if it > 1
##      h(it) = abs(p - p0) / abs(p);
##      if h(it) < tol
##        break;
##      endif
##    else
##      h(it) = abs(p - p0);
##    endif



    p0 = p;
    it++;
  endwhile

  if it > maxit
    fprintf("No se encontro punto fijo en maxit iteraciones en: %s\n", func2str(f));
    %disp("no se encontro punto fijo en maxit iteraciones");
  endif

  t = toc();
endfunction

%function [x, h] = puntofijo(g, x0, kmax, tol)
%  tic();
%  h = [];  % Vector de errores
%  k = 0;
%
%  while k < kmax
%    x = g(x0);
%    h(k+1) = abs(x - x0);  % error actual%
%
%    if h(k+1) < tol
%      break;
%    endif
%
%    k++;
%    x0 = x;
%  endwhile
%
%  if k == kmax
%    disp('No se encontró el punto fijo en kmax iteraciones.');
%  endif
%
%  t = toc();
%endfunction
