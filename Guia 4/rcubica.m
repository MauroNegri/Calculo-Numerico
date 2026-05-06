function x = rcubica(a)
  tol = 1e-12;
  maxit = 1000;
  if a == 0; x = 0; return; endif
  f = @(t) t.^3 - a;
  if a > 1
    xmin = 0; xmax = a;
  elseif a > 0
    xmin = 0; xmax = 1;
  else
    xmin = a; xmax = 0;
  endif
  [x, ~, ~, ~] = biseccion(f, xmin, xmax, maxit, tol);
endfunction
