function [L D U]=descomponer_matriz(A)
  n=length(A);
  v=diag(A);
  L=tril(A,-1);
  U=triu(A,1);
  D=diag(v);
endfunction
