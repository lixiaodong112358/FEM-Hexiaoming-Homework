function result=function_f1(x,y,t)
%%ÐÎ³Évector b1µÄunction
v=function_v(x,y);
result=-2*pi*(x^2*y^2+exp(-y))*sin(2*pi*t)+(-2*v*x^2-2*v*y^2-v*exp(-y)+pi^2*cos(pi*x)*cos(2*pi*y))*cos(2*pi*t);
end