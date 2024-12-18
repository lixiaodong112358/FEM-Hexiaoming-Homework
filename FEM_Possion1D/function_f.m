function [result]=function_f(x)
%%%Xiaodong LI 2021/4/21
result=-exp(x).*(cos(x)-2.*sin(x)-x.*cos(x)-x.*sin(x));
end