function result=function_f2(x,y,t)
%%ÐÎ³Évector b2µÄfunction
result=-2*pi*(-2/3*x*y^3+2-pi*sin(pi*x))*sin(2*pi*t)+...
    (4*function_v(x,y)*x*y-function_v(x,y)*pi^3*sin(pi*x)+2*pi*(2-pi*sin(pi*x))*sin(2*pi*y))*cos(2*pi*t);
end