%%%用于画图的脚本
close all
[xx,yy]=meshgrid(left:h1:right,top:-h2:bottom);

zz=flip(reshape(solution,[N2+1,N1+1]));

figure

surf(xx,yy,zz)
% figure
% hold on
% plot3(Pb(1,:),Pb(2,:),solution,'ro')