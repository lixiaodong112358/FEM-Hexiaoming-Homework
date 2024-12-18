function result=function_f2(x,y,t)
%%非稳态NS方程2021/9/13
%%形成vector b2的function
result=cos(2*pi*t)^2*(x^2*y^2+exp(-y))*(-pi^2*cos(pi*x)-2/3*y^3)-...
    2*pi*sin(2*pi*t)*(-2/3*x*y^3-pi*sin(pi*x)+2)-2*x*y^2*cos(2*pi*t)^2*(-2/3*x*y^3-pi*sin(pi*x)+2)...
    +function_v(x,y)*(4*x*y*cos(2*pi*t)-pi^3*cos(2*pi*t)*sin(pi*x))-2*pi*cos(2*pi*t)*(pi*sin(pi*x)-2)*sin(2*pi*y);
end