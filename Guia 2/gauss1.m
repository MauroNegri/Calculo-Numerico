function [x] = gauss1(A,b)
  n=length(b);
  A=[A b];
  % Eliminacion
    for k=1:n-1
      for i=k+1:n
        m = A(i,k)/A(k,k);
        A(i,k)=0;
        for j=k+1:n+1
          A(i,j) = A(i,j)-m*A(k,j);
        endfor
      endfor
    endfor
    if (A(n,n)==0)
      disp(’no hay sol. unica’)
    endif
  x=sust_atras1(A); %retrosustitucion
endfunction
