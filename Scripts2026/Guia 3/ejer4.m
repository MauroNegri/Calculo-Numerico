A = [10  5  0  0;
      5 10 -4  0;
      0 -4  8 -1;
      0  0 -1  5];
b = [6; 25; -11; -11];
x0 = [0; 0; 0; 0];
maxit = 1000;
tol =1e-6;
##w = 1.1535;

## Parámetro óptimo w
rho_Tj = 0.6793;
w_opt = 2 / (1 + sqrt(1 - rho_Tj^2));
disp("Valor óptimo de w:");
disp(w_opt);
## Cuándo es válida la fórmula?
## Matriz definida positiva — todos sus autovalores son positivos
## Matriz consistentemente ordenada — TRIDIAGONAL
## Radio espectral de Jacobi real — rho_Tj debe ser real, lo que se garantiza con las condiciones anteriores

# si la matriz no cumple esas condiciones
% Estimación numérica de w_opt
w_vals = linspace(0.01, 1.99, 100);
iters = zeros(size(w_vals));

for j = 1:length(w_vals)
    [~, iters(j), ~] = sor(A, b, x0, maxit, tol, w_vals(j));
end

[~, idx] = min(iters);
w_opt_num = w_vals(idx);
disp("w óptimo numérico:"); disp(w_opt_num)

##  Rango       Nombre                Cuándo usarlo
##0 < w < 1   Subrelajación       Gauss-Seidel diverge o oscila
##  w = 1     Gauss-Seidel        caso particular de SOR
##1 < w < 2   Sobrelajación       caso de mayor interés, acelera convergencia
[X_SOR, it_SOR, rh_SOR] = sor(A,b, x0, maxit, tol, w_opt);
[X_SOR2, it_SOR2, rh_SOR2] = sor(A,b, x0, maxit, tol, w_opt_num);

% Radios espectrales
[rho_SOR_opt] = radio_espectral(A, 'sor', w_opt);
[rho_SOR_optnum] = radio_espectral(A, 'sor', w_opt_num);
[rho_GS] = radio_espectral(A,'gs', []);

disp("Radio espectral Gauss-Seidel:");
disp(rho_GS);
disp("Radio espectral SOR:");
disp(rho_SOR_opt);
disp("Radio espectral SOR:");
disp(rho_SOR_optnum);

# Mostrar resultados gauss seidel
[X_GS, it_GS, rh_GS] = gauss_seidel(A, b, x0, maxit, tol);
disp(["GS Iteraciones: ", num2str(it_GS)]);
disp(["Error final: ", num2str(rh_GS(end))]);

% Mostrar resultados SOR
disp("Método SOR:");
disp(["Iteraciones: ", num2str(it_SOR)]); % printf("Iteraciones: %d\n", itSOR); Más preciso y limpio
disp(["Error final: ", num2str(rh_SOR(end))]);

% Mostrar resultados SOR2
disp("Método SOR2:");
disp(["Iteraciones: ", num2str(it_SOR2)]); % printf("Iteraciones: %d\n", itSOR); Más preciso y limpio
disp(["Error final: ", num2str(rh_SOR2(end))]);

# Ambos SOR convergen:
# El w numérico(w_opt_num) 1.17 converge en 1 iteración menos, pero con un error final mayor
# El w analítico(w_opt) 1.1535 necesita 1 iteración más, pero llega a un error más chico
# Esto pasa porque el w numérico minimiza el número de iteraciones hasta
# cruzar la tolerancia 1e-6, pero no necesariamente el error en sí.
# El w analítico es el óptimo teórico real.

# Evaluar w segun el libro
n = length(b);
w_vals = linspace(1.0, 1.99, 100); # entre 10 y 20
tasa = zeros(size(w_vals));

for j = 1:length(w_vals)
    x = x0;
    x0_ = x0;
    for it = 1:15  % solo 15 iteraciones
        % hacer una iteracion SOR manualmente
        for i = 1:n
            x(i) = (1-w_vals(j))*x0_(i) + w_vals(j)*(b(i) - A(i,1:i-1)*x(1:i-1) ...
                   - A(i,i+1:n)*x0_(i+1:n)) / A(i,i);
        end
        if it == 14
            diff_prev = norm(x - x0_, 'inf');
        end
        if it == 15
            diff_curr = norm(x - x0_, 'inf');
        end
        x0_ = x;
    end
    tasa(j) = diff_curr / diff_prev;  % cociente de reduccion
end

[~, idx] = min(tasa);
disp("w óptimo metodo del libro:"); disp(w_vals(idx))
