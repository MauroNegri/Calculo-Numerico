function [rhoA] = radio_espectral(A, metodo,w)
  % Descomponemos la matriz A para calcular la matriz
  #[L D U] = DescomponerMatriz(A);
  n=length(A);
  v=diag(A);
  L=tril(A,-1);
  U=triu(A,1);
  D=diag(v);
  % Seleccion del metodo
  if strcmp(metodo, "ja")
    rhoA = max(abs(eig(-inv(D)*(L+U))));
 # rhoA = abs(eigs(-inv(D)*(L+U)), 1);
   elseif strcmp(metodo, "gs")
    rhoA = max(abs(eig(-inv(D+L)*U)));
 # rhoA = abs(eigs(-inv(D+L)*U), 1);

   elseif strcmp(metodo, "sor")
    rhoA = max(abs(eig(-inv(D+w*L)*((w-1)*D+w*U))));
 %   rhoA = (abs(eigs(-inv(D+w*L)*((w-1)*D+w
  endif
  # Para matrices muy grandes N >= 100 usar la linea comentada eigs(M,1)
  # Para matrices mas chicas N < 100 usar eig(M)
  # eig(M) -> devuelve todos los autovalores
  # eigs(M,1) -> devuelve solo el mayor autovalor
endfunction
