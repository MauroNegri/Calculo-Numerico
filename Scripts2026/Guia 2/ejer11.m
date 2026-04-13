A = [1 1+0.5e-15 3; 2 2 20; 3 6 4];
n = 3;

# SIN PIVOTEO
[L_a, U_a] = doolittle(A);
disp('Sin pivoteo')
disp('L ='); disp(L_a)
disp('U ='); disp(U_a)
disp('Residuo A - L*U:'); disp(norm(A - L_a*U_a))

# CON PIVOTEO
[Af, r] = doolittle_p(A);
L_p = tril(Af(r,:), -1) + eye(n);
U_p = triu(Af(r,:));
P = eye(n)(r,:);
disp('Con pivoteo')
disp('L ='); disp(L_p) # L: matriz triangular inferior con 1s en la diagonal, contiene los multiplicadores usados durante la eliminación
disp('U ='); disp(U_p) # U: matriz triangular superior, contiene los pivotes y los coeficientes eliminados
disp('r ='); disp(r) # r: vector de permutaciones, r(k) indica qué fila original de A quedó en la posición k tras el pivoteo, en este caso r = [3 2 1] significa: primero fila 3, luego fila 2, luego fila 1
disp('P ='); disp(P) # P: matriz de permutación construida a partir de r, PA reordena las filas de A según el pivoteo parcial
disp('Residuo A - L*U:');  disp(norm(A - L_p*U_p))
disp('Residuo PA - L*U:'); disp(norm(P*A - L_p*U_p))

# Si una implementación sin pivoteo te da un residual grande,
# pero la con pivoteo te da un residual cercano a cero, quiere decir que
# el pivoteo fue necesario para tener una factorización confiable.

# Cuando no usás pivoteo:
# Estás asumiendo que el elemento en la diagonal es el mejor candidato
# posible para dividir, lo cual no es cierto si hay números pequeños
# en esa posición.
# Esa división con un número muy cercano a cero amplifica los errores numéricos
# porque los multiplicadores se vuelven enormes (del orden de 1/pivote),
# distorsionando los elementos restantes de la matriz.

# Con pivoteo:
# Siempre buscás el mayor valor absoluto en la columna actual
# (desde la fila actual hacia abajo), lo que minimiza el error.
