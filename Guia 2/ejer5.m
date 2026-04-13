S = [ 16.87 0.1650 0.2019 0.3170 0.2340   0.1820  0.1100;
         0   27.70 0.8620 0.0620 0.0730   0.1310  0.1200;
         0       0  22.35  13.05  4.420    6.001  0.3710;
         0       0      0  11.28      0    1.110  0.3710;
         0       0      0      0  9.850   1.1684   2.108;
         0       0      0      0  0.2990   15.98   2.107;
         0       0      0      0       0       0   4.670];
         # Matriz de coeficientes de sensibilidad (de la tabla)
h = [17.1; 65.1; 186.0; 82.7; 84.2; 63.7; 119.7]; # Vector de alturas de picos (dado en el enunciado)

[p] = gauss(S,h);

P_total = sum(p);
#P_total = 41.738
# Comparar con la presión medida
disp("Presión total calculada: "), disp(P_total);
disp("Presión medida en el experimento: 38.78 μm Hg");

# Error relativo = (| Valor Calculado - Valor Real  / Valor real |) * 100
# Error relativo = 7.63%

error_rel = abs(P_total - 38.78) / 38.78 * 100;
disp("Error relativo (%): "), disp(error_rel)

t = linspace(0, 1, 7);
plot(t, p);
xlabel('t'); ylabel('p(t)');
title('Solución del sistema');
