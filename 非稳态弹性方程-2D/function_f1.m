function result=function_f1(x,y,t)
%%ÐÎ³Évector b1µÄfunction
lambda=function_lambda(x,y);
mu=function_nu(x,y);
result=pi^2.*(lambda+3*mu-1)*sin(pi*t)*sin(pi*x)*sin(pi*y)-(2*x-1)*(2*y-1)*(lambda+mu)*cos(t);
end