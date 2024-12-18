function result=function_f2(x,y)
%%ÐÎ³Évector b2µÄfunction
result=-(function_lambda(x,y)+2*function_nu(x,y))*(2*x*(x-1))-(function_lambda(x,y)+function_nu(x,y))*(pi^2*cos(pi*x)*cos(pi*y))...
    -function_nu(x,y)*(2*y*(y-1));
end