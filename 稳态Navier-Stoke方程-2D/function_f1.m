function result=function_f1(x,y)
%%ÐÎ³Évector b1µÄunction
result=-2*function_v(x,y)*x^2-2*function_v(x,y)*y^2-function_v(x,y)*exp(-y)+pi^2*cos(pi*x)*cos(2*pi*y)...
       +2*x*y^2*(x^2*y^2+exp(-y))+(-2*x*y^3/3+2-pi*sin(pi*x))*(2*x^2*y-exp(-y));
end