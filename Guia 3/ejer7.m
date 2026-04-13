n = 10;
Bin =[ 0.25*ones(n,1), 0.5*ones(n,1), 2*ones(n,1), 0.5*ones(n,1), 0.25*ones(n,1)];
d = [-4, -2, 0, 2, 4];
A = spdiags(Bin, d, n, n);
A = full(A);
##disp(eig(A));

## Con diag()
##n = 10;
##A = zeros(n,n);
##main_diag = 2 * (1:n);   # Diagonal principal
##subdiag2 = 0.5 * (3:n);  # Segunda subdiagonal (empieza en i=3)
##subdiag4 = 0.25 * (5:n);
##updiag2 = 0.5 * (1:n-2);
##updiag4 = 0.25 * (1:n-4);
##
##A = diag(subdiag4, -4) + ...
##    diag(subdiag2, -2) + ...
##    diag(main_diag) + ...
##    diag(updiag2, 2) + ...
##    diag(updiag4, 4);

b = pi*ones(n,1);
x0 = zeros(n,1);
tol = 1e-5;
maxit = 100;
##w_opt = 1.00621;

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


% Radios espectrales
[rho_JA] = radio_espectral(A, 'ja', []);
[rho_GS] = radio_espectral(A,'gs', []);
w_opt = 2 / (1 + sqrt(1 - rho_JA^2));
disp("Valor óptimo de w:"); # Matriz definida positiva, trigiagonal y el rho_JA existe y es real
disp(w_opt);
[rho_SOR] = radio_espectral(A, 'sor', w_opt); # Calculo el Radio Espectral de Jacobi y con eso saco w

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
disp(["Iteraciones: ", num2str(it_SOR)]); % printf("Iteraciones: %d\n", itSOR); Más preciso y limpio
disp(["Error final: ", num2str(rh_SOR(end))]);
disp(["Tiempo: ", num2str(t_SOR), " segundos"]);

# Gráfico comparativo de errores
figure;
semilogy(1:length(rh_JA), rh_JA, 'm-^', 'DisplayName', 'Jacobi');
hold on;
semilogy(1:length(rh_GS), rh_GS, 'b-o', 'DisplayName', 'Gauss-Seidel');
hold on;
semilogy(1:length(rh_SOR), rh_SOR, 'r-s', 'DisplayName', 'SOR');

xlabel('Iteración');
ylabel('Error relativo');
title('Comparación de convergencia: Jacobi vs Gauss-Seidel vs SOR');
legend('show');
grid on;

## - El método que resultó más eficiente fue el [Gauss-Seidel / SOR], ya que alcanzó
##   la tolerancia con menos iteraciones y en menor tiempo.
## - Jacobi mostró una tasa de convergencia más lenta, como se esperaba por su
##   forma de actualización completamente desacoplada.
## - La convergencia de los métodos está justificada por el hecho de que la matriz A es
##   diagonalmente dominante y simétrica positiva definida, lo cual garantiza convergencia
##   para estos métodos.
