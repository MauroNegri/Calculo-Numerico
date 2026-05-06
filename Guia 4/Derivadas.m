%clc;
pkg load symbolic
syms x
%f = @(x) 3*(x+0.5)*sin((x-2.7)/2)^4;
%f = 3*(x + 1/2)*sin((x - 27/10)/2)^4;
%f = 3*(x + sym(1)/2)*sin((x - sym(27)/10)/2)^4;
%f = 3*(x+0.5)*sin((x-2.7)/2)^4;
%df = diff(f, x)

%f = 1 / (3*exp(x) + 4*x + 5);

f =  x + exp(-10 *x^2) * cos(x);

df = diff(f, x)
%syms x y
%f = x^2 * y + sin(y);
%df_dx = diff(f, x)  % derivada parcial respecto a x
%df_dy = diff(f, y)  % derivada parcial respecto a y

%d2f = diff(f, x, 2);  % segunda derivada respecto a x

