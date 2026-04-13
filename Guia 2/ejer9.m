A = [1 -2 3 0; 3 -6 9 3; 2 1 4 1; 1 -2 2 -2];
n = 4;

[Af, r] = doolittle_p(A);

L = tril(Af(r,:), -1) + eye(n);
U = triu(Af(r,:));

disp('L ='); disp(L)
disp('U ='); disp(U)
disp('r ='); disp(r)

# Construir P a partir de r
P = eye(n)(r,:);
disp('P ='); disp(P)

# Verificar PA = L*U
disp('PA ='); disp(P*A)
disp('L*U ='); disp(L*U)
disp('Residuo PA - L*U:'); disp(norm(P*A - L*U))
# norm(P*A - L*U) Mide qué tan parecidos son PA y LU.
# Si la factorización está bien, este número debe ser cero
# (o muy cercano a 0, en el orden de 1e-14 o menor por errores de redondeo).

