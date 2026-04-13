b = [1; 2; 3];

# MATRIZ A1
# A1 no tiene problemas de pivote: usamos doolittle normal
A1 = [8 4 1; 2 6 2; 2 4 8];

[L1_a, U1_a] = doolittle(A1);
disp('=== A1 - doolittle ===')
disp('L1_a ='); disp(L1_a)
disp('U1_a ='); disp(U1_a)

# Resolver L1*y = b, luego U1*x = y
y1_a = sust_adelante(L1_a, b);
x1_a = sust_atras(U1_a, y1_a);
disp('x1_a ='); disp(x1_a)

r1_a = A1*x1_a - b;
disp('Residuo A1:'); disp(norm(r1_a))


[A1f, r1] = doolittle_p(A1);
n = 3;
L1 = tril(A1f(r1,:), -1) + eye(n);
U1 = triu(A1f(r1,:));
disp('=== A1 - doolittle_p ===')
disp('L1 ='); disp(L1)
disp('U1 ='); disp(U1)

# Resolver L2*y = b(r2), luego U2*x = y
y1 = sust_adelante(L1, b(r1));
x1 = sust_atras(U1, y1);
disp('x1 ='); disp(x1)

r1 = A1*x1 - b;
disp('Residuo A1:'); disp(norm(r1))


# MATRIZ A2
# A2 tiene valores de muy distinta magnitud: conviene doolittle_p
A2 = [5.00e-02 5.57e+02 -4.00e+02;
      1.98e+00 1.94e+02 -3.00e-03;
      2.74e+02 3.11e+00  7.50e-02];

[L2_a, U2_a] = doolittle(A2);
disp('=== A2 - doolittle ===')
disp('L2_a ='); disp(L2_a)
disp('U2_a ='); disp(U2_a)

# Resolver L1*y = b, luego U1*x = y
y2_a = sust_adelante(L2_a, b);
x2_a = sust_atras(U2_a, y2_a);
disp('x2_a ='); disp(x2_a)

r2_a = A2*x2_a - b;
disp('Residuo A2_a:'); disp(norm(r2_a))

[A2f, r2] = doolittle_p(A2);
n = 3;
L2 = tril(A2f(r2,:), -1) + eye(n);
U2 = triu(A2f(r2,:));
disp('=== A2 - doolittle_p ===')
disp('L2 ='); disp(L2)
disp('U2 ='); disp(U2)

# Resolver L2*y = b(r2), luego U2*x = y
y2 = sust_adelante(L2, b(r2));
x2 = sust_atras(U2, y2);
disp('x2 ='); disp(x2)

r2 = A2*x2 - b;
disp('Residuo A2:'); disp(norm(r2))

# MATRIZ A3
# A3 tiene pivote cero en paso 2 (sin pivoteo falla): usar doolittle_p
A3 = [1 2 -1; 2 4 0; 0 1 -1];

[A3f, r3] = doolittle_p(A3);
L3 = tril(A3f(r3,:), -1) + eye(n);
U3 = triu(A3f(r3,:));
disp('=== A3 - doolittle_p ===')
disp('L3 ='); disp(L3)
disp('U3 ='); disp(U3)

# Resolver L3*y = b(r3), luego U3*x = y
y3 = sust_adelante(L3, b(r3));
x3 = sust_atras(U3, y3);
disp('x3 ='); disp(x3)

r3 = A3*x3 - b;
disp('Residuo A3:'); disp(norm(r3))

[L3_a, U3_a] = doolittle(A3);
disp('=== A3 - doolittle ===')
disp('L3_a ='); disp(L3_a)
disp('U3_a ='); disp(U3_a)

# Resolver L1*y = b, luego U1*x = y
y3_a = sust_adelante(L3_a, b);
x3_a = sust_atras(U3_a, y3_a);
disp('x3_A ='); disp(x3_a)

r3_a = A3*x3_a - b;
disp('Residuo A3_a:'); disp(norm(r3_a))
