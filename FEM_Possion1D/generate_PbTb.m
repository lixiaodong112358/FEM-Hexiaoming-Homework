function [Pb,Tb]=generate_PbTb(a,b,h,ID)
N=(b-a)/h;
if ID==101
    Pb=a:h:b;
    Nm=N+1;
    Tb=[1:Nm-1;2:Nm];
elseif ID==102
    Pb=a:h/2:b;
    Tb=[1:2:2*N-1;
        3:2:2*N+1;
        2:2:2*N];
    %%
    % _ingbu_
end
end