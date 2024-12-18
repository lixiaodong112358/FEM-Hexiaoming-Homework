function result=FE_basis_local_fun_1D(xi,vertices,basis_type,basis_index,basis_der_x)
% trial 和test 有可能不一样，届时可能需要改动程序
% xi: 高斯节点
%101: 1D linear
%102: 1D quadratic
h=vertices(2)-vertices(1);
    if basis_type==101
        if basis_index==1
            if basis_der_x==0
                result=(vertices(2)-xi)/h;
            elseif basis_der_x==1
                result=-1/h;
            elseif basis_der_x>=2&&rem(basis_der,1)==0%%%整数
                result=0;
            else
                warning='wrong input for basis derivative order'
            end
        elseif basis_index==2
            if basis_der_x==0
                result=(xi-vertices(1))/h;
            elseif basis_der_x==1
                result=1/h;
            elseif basis_der_x>=2&&rem(basis_der,1)==0%%%整数
                result=0;
            else
                warning='wrong input for basis derivative order'
            end
        else
            warning='wrong input for basis index'
        end
    elseif basis_type==102
        %%%基函数是非线性基函数的情况
        if basis_index==1
            if basis_der_x==0
                xx=(xi-vertices(1))/h;
                result=2*xx^2-3*xx+1;
            elseif basis_der_x==1
                result=-3/h+4*(xi-vertices(1))/h/h;
            elseif basis_der_x==2
                result=4/h/h;
            elseif basis_der_x>2&&rem(basis_der_x,1)==0
                result=0;
            else
                warning='基函数导数阶输入错误，请检查'
            end
        elseif basis_index==2
            if basis_der_x==0
                xx=(xi-vertices(1))/h;
                result=2*xx^2-xx;
            elseif basis_der_x==1
                result=-1/h+4*(xi-vertices(1))/h/h;
            elseif basis_der_x==2
                result=4/h/h;
            elseif basis_der_x>2&&rem(basis_der_x,1)==0
                result=0;
            else
                warning='基函数导数阶输入错误，请检查'
            end
        elseif basis_index==3
            if basis_der_x==0
                xx=(xi-vertices(1))/h;
                result=-4*xx^2+4*xx;
            elseif basis_der_x==1
                result=4/h-8*(xi-vertices(1))/h/h;
            elseif basis_der_x==2
                result=-8/h/h;
            elseif basis_der_x>2&&rem(basis_der_x,1)==0
                result=0;
            else
                warning='基函数导数阶输入错误，请检查'
            end
        else
            warning='非线性基函数参数输入错误，请检查'
        end
    end
    
end