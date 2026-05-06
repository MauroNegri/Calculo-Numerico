function [p, it, r, t] = newton(f, df, p0, maxit, tol)
  tic();
  r = zeros(1, maxit);  % Prealoca el vector de errores
  it = 1;

  while it < maxit
    fp0 = f(p0);
    dfp0 = df(p0);

    if dfp0 == 0
      error("Derivada cero en p0. No se puede continuar.");
    endif

    %usamos una condicion de corte
    %   abs(f(p)) < tolerancia
    %   abs(p-p0) < tolerancia
    %   abs(p-p0)/abs(p) < tolerancia

    p = p0 - fp0 / dfp0;

    ## Criterio de corte 1: |f(p)| < tol
    ## r(it) = abs(f(p));
    ## if r(it) < tol
    ##   break;
    ## endif

    ## Criterio de corte 2: |p - p0| < tol
     r(it) = abs(p - p0);
     if r(it) < tol
       break;
     endif

    ## Criterio de corte 3: |p - p0| / |p| < tol
##    r(it) = abs(p - p0) / abs(p);
##    if r(it) < tol
##      break;
##    endif

    p0 = p;
    it = it + 1;
  endwhile

  if it == maxit
    fprintf("No converge en %d iteraciones\n", maxit);
  endif

  t = toc();
  r = r(1:it);  % Recorta el vector a las iteraciones reales
endfunction

##function [x,it,rh,t] = newton(f,df,x0,kmax,tol)
##tic();
##it = 0;
##fx=f(x0);
##
##while it < kmax
##  it++;
##    x = x0 -(fx/df(x0));
##    #fx0 = f(x0);
##    fx = f(x);
##    %abs(fx);
##    %rh(it) = max([abs(fx), abs((x-x0)/x)])
##    rh(it) = abs((x-x0)/x);
##    %rh(it) = abs((x-x0));
##   %E_rela = 2err/(abs(x)+tol);
##
##    if rh(it) < tol
##      break;
##    endif
##
##    %if(E_rel(it)<tol)
##      %break;
##    %endif
##
##    %if(f(p) < tol)
##       %  break
##       %endif
##
## x0 = x;
##
##endwhile
##t =toc();
##endfunction
##
##%La funcion tiene un cero si y solo si f(p) = 0 y f'(p) != 0
