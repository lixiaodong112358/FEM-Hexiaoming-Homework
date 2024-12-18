function result=function_f2(x,y,t)
%%ÐÎ³Évector b2µÄfunction
lambda=function_lambda(x,y);
mu=function_nu(x,y);
result=cos(t)*(2*mu*(-2*(x-1)*x-y^2+y)-(x-1)*x*(2*lambda+(y-1)*y))...
       -pi^2*(lambda+mu)*sin(pi*t)*cos(pi*x)*cos(pi*y);
end