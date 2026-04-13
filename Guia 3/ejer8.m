N_h = 10;           # Cantidad de nodos (puede ser 10, 20, 50, etc.)
n = N_h - 3;        # Tamaño del sistema (desconocidas eta_2 a eta_{N_h - 2})
B = [ 1*ones(n,1),  -4*ones(n,1), 6*ones(n,1), -4*ones(n,1), 1*ones(n,1)];
d = [-2, -1, 0, 1, 2];
A = spdiags(B,d,n,n);
A = full(A);
##disp(A);
##disp(eig(A)); definida positiva autovalores positivos
x0 = zeros(n,1);
tol = 1e-8;
##w_opt = 1.7;
w_opt = 1.6300;
maxit = 1000;

L = 5;             # Longitud de la viga en metros
E = 210e+3;        # Módulo de Young en Pascales (210 GPa)
w = 0.07;          # Ancho de la viga en metros
s = 0.14;          # Altura de la viga en metros
P = 1.0326e+4;     # Cte
J = w * s^3 / 12;  # Momento de inercia
h = L / N_h;       # Paso de discretización
f = P / (E*J);

b = h^4*f * (1:n)';

# Método de Gauss
[x] = gauss(A,b);
disp("Sistema Resuelto por Gauss");
disp(x);

##w_vals = linspace(0.01, 1.99, 100);
##iters = zeros(size(w_vals));
##
##for j = 1:length(w_vals)
##    [~, iters(j), ~] = sor(A, b, x0, maxit, tol, w_vals(j));
##end
##
##[~, idx] = min(iters);
##w_opt_num = w_vals(idx);
##disp("w óptimo numérico:"); disp(w_opt_num)
##
### Evaluar w segun el libro
##n = length(b);
##w_vals = linspace(1.0, 1.99, 100);
##tasa = zeros(size(w_vals));
##
##for j = 1:length(w_vals)
##    x = x0;
##    x0_ = x0;
##    for it = 1:15  % solo 15 iteraciones
##        % hacer una iteracion SOR manualmente
##        for i = 1:n
##            x(i) = (1-w_vals(j))*x0_(i) + w_vals(j)*(b(i) - A(i,1:i-1)*x(1:i-1) ...
##                   - A(i,i+1:n)*x0_(i+1:n)) / A(i,i);
##        end
##        if it == 14
##            diff_prev = norm(x - x0_, 'inf');
##        end
##        if it == 15
##            diff_curr = norm(x - x0_, 'inf');
##        end
##        x0_ = x;
##    end
##    tasa(j) = diff_curr / diff_prev;  % cociente de reduccion
##end
##
##[~, idx] = min(tasa);
##disp("w óptimo metodo del libro:"); disp(w_vals(idx))

# Radios espectrales
[rho_JA] = radio_espectral(A, 'ja', []); # radio espectral > 1, diverge
[rho_GS] = radio_espectral(A,'gs', []); # radio espectral cercano a 1, converge pero lento
w_opt2 = 2 / (1 + sqrt(1 - rho_JA^2)); # radio espectral cercano a 1 pero mas chico que gs, converge lento
disp("Valor óptimo de w:"); # Matriz definida positiva, trigiagonal y el rho_JA real
disp(w_opt2); # Da un valor imaginario, no se puede usar
[rho_SOR] = radio_espectral(A, 'sor', w_opt);

disp("Radio espectral Jacobi:");
disp(rho_JA);
disp("Radio espectral Gauss-Seidel:");
disp(rho_GS);
disp("Radio espectral SOR:");
disp(rho_SOR);

# Métodos iterativos
[x_JA, it_JA, rh_JA, t_JA] = jacobi(A, b, x0, maxit, tol);
[x_GS, it_GS, rh_GS, t_GS] = gauss_seidel(A, b, x0, maxit, tol);
[x_SOR, it_SOR, rh_SOR, t_SOR] = sor(A,b, x0, maxit, tol, w_opt);

# Mostrar resultados Jacobi
disp("Método Jacobi:");
disp(["Iteraciones: ", num2str(it_JA)]);
disp(["Error final: ", num2str(rh_JA(end))]);
disp(["Tiempo: ", num2str(t_JA), " segundos"]);

# Mostrar resultados Gauss-Seidel
disp("Método Gauss-Seidel:");
disp(["Iteraciones: ", num2str(it_GS)]);
disp(["Error final: ", num2str(rh_GS(end))]);
disp(["Tiempo: ", num2str(t_GS), " segundos"]);

# Mostrar resultados SOR
disp("Método SOR:");
disp(["Iteraciones: ", num2str(it_SOR)]);
disp(["Error final: ", num2str(rh_SOR(end))]);
disp(["Tiempo: ", num2str(t_SOR), " segundos"]);

## La matriz no es diagonalmente dominante
## diagonal = 6, suma off-diagonal = 4+4+1+1 = 10 > 6
## por eso Jacobi no converge. Gauss-Seidel y SOR convergen igual porque la
## matriz es simétrica definida positiva.
##
## Método directo de Gauss — pivoteo necesario?
## No es necesario aplicar pivoteo porque la matriz es simétrica definida positiva
## (todos sus autovalores son positivos, verificable con eig(A)). Esto garantiza
## que todos los pivotes durante la eliminación son estrictamente positivos,
## por lo que Gauss sin pivoteo funciona correctamente y produce la solución
## de referencia.
## x_gauss = [3046.3, 7257.4, 10944, 12800, 12096, 8870.2, 4121.5]
## x_JA = [-5.3507e+170, 9.0013e+170, -1.1542e+171, 1.2417e+171, -1.1542e+171, 9.0013e+170, -5.3507e+170]
## x_GS = [3046.3, 7257.4, 10944, 12800, 12096, 8870.2, 4121.5]
## x_SOR = [3046.3, 7257.4, 10944, 12800, 12096, 8870.2, 4121.5]
##
## Métodos iterativos
## - Jacobi diverge porque el radio espectral ρ = 1.481 > 1. La causa es que la
## matriz no es diagonalmente dominante, condición que Jacobi requiere para converger.
## - Gauss-Seidel converge a pesar de que la matriz no es diagonalmente dominante,
## porque la condición suficiente para GS es más débil: basta con que la matriz
## sea simétrica definida positiva. Sin embargo con ρ = 0.9743 muy cercano a 1,
## la convergencia es lenta (569 iteraciones).
## - SOR con w=1.63 es el método más eficiente: converge en solo 90 iteraciones,
## siendo 6 veces más rápido que Gauss-Seidel.
