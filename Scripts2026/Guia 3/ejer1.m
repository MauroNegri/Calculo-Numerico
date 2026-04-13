% Ejercicio 1.
% i)- Calcule los autovalores y el radio espectral de la matriz A,
% ii)- defina matriz convergente e indique si lo es la matriz A.

A = [ 3/4 1/6; 1/4 1/4];

% Autovalores
autovalores = eig(A);

% Radio espectral: el mayor valor absoluto de los autovalores

radio_espectral = max(abs(autovalores));

% Una matriz A se dice convergente si su radio espectral es menor que 1, es decir:
% p(A) = max|lambda| < 1
% La matriz A es convergente
