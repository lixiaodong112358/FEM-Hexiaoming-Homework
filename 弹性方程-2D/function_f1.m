function result=function_f1(x,y)
%%ÐÎ³Évector b1µÄunction
result=-(function_lambda(x,y)+2*function_nu(x,y))*(-pi^2*sin(pi*x)*sin(pi*y))...
-(function_lambda(x,y)+function_nu(x,y))*((2*x-1)*(2*y-1))-function_nu(x,y)*(-pi^2*sin(pi*x)*sin(pi*y));
end