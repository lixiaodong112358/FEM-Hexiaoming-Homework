function result=function_f1(x,y,t)
%%非稳态NS方程2021/9/13
%%形成vector b1的unction
result=cos(2*pi*t)^2*(2*x^2*y-exp(-y))*(-2/3*x*y^3-pi*sin(pi*x)+2)...
    +function_v(x,y)*(-cos(2*pi*t)*(2*x^2+exp(-y))-2*y^2*cos(2*pi*t))-2*pi*sin(2*pi*t)*(x^2*y^2+exp(-y))...
    +2*x*y^2*cos(2*pi*t)^2*(x^2*y^2+exp(-y))+pi^2*cos(2*pi*t)*cos(pi*x)*cos(2*pi*y);
end