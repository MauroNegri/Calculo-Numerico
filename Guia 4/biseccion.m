%function [p,h] = biseccion(f,a,b,maxit,tol)
function [p,it,h,t] = biseccion(f,a,b,maxit,tol)
  tic();
  fa = f(a);
  fb = f(b);

  if sign(fa) * sign(fb) > 0
    %fprintf("En la funcion: %s\n", func2str(f));
    error("No se cumple la regla de los signos");
  endif

  it = 0;
  prev_p = a;  % para calcular el error relativo
  h = [];  % Inicializar el vector de errores

  while (it < maxit)
    it++;
    p = a+(b-a)/2;
    fp = f(p);
## Criterio de corte 1: |f(pn)| <= tol
##    h(it) = abs(fp);
##    if (h(it) < tol)
##     break;
##   endif
##
## Criterio de corte 2: |pn - pn-1| <= tol
     h(it) = abs(p - prev_p);
     if (it > 1 && abs(p - prev_p) < tol)
       break;
     endif
##
## Criterio de corte 3: |pn - pn-1| / |pn| <= tol
##     if (it > 1 && abs(p - prev_p) / abs(p) < tol)
##       break;
##     endif
## Criterio de corte 3: |pn - pn-1| / |pn| <= tol
##    if it > 1
##      h(it) = abs(p - prev_p) / abs(p);
##      if h(it) < tol
##        break;
##      endif
##    else
##      h(it) = abs(p - prev_p);  % primera iteracion no hay prev_p valido
##    endif



% Actualizar los límites del intervalo
    if sign(fp)*sign(fb)<0
      fa = fp;
      a = p;
    else
      fb = fp;
      b = p;
    endif

    prev_p = p;  % Actualizar prev_p para la siguiente iteración
  endwhile

  t = toc();
 endfunction

 % posee convergencia lineal quiere decir que es muy lento, pero es muy robusto asique siempre converge
% Las raices son -1, 0.5, 1
% 1) Verificar los extremos [-2,1.5] -> f(-2) = -22.5 < 0 y f(1.5) = 3.75 > 0. Cambio de signo. Puedo usar Biseccion
% 2) iteracion 1 - Punto Medio del intervalo
% P1 =  (-2 + 1.5) / 2 = -0.25
% f(-0.25) = 2.109 > 0
% Nuevo intervalo: [-2, -0.25] f(-2) < 0 | f(-0.25) > 0
% iteracion 2
% P2 = (-2 - 0.25) / 2 = -1.125
% f(-1.125) = -1.2949 < 0
% Nuevo intervalo: [-1.125, -0.25] f(-1.125) < 0 | f(-0.25) > 0
% iteracion 3
% P3 = (-1.125 - 0.26) / 2 = -0.6875
% f(-0.6875) = 1.873 > 0
% Nuevo intervalo: [-1.125, -0.6875] f(-1.125 < 0 | f(-0.6875) > 0
% iteracion 4
% P4 = (-1.125 - 0.6875) / 2 = -0.90625
% f(-0.90625) = 0.755 > 0
% SI no logro identificar hacer una iteracion mas
% P5 = (-1.125 - 0.90625) / 2 = -1.0156
% La raiz a la que estoy convergiendo es a -1.
% E <= (b-a) / 2^(w * n^2 * it)
% E <= (1.5 - (-2))/ 2^4 -> E <= 0.21875
