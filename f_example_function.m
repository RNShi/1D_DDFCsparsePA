%% Test Functions%%
% This script gives some test functions:
% 
% f1: a smooth function sin(pi*x) (no jump);
% f2: a function has one jump;
% f3: a function has four jumps with similar maginitude;
% f4: a function has two large jumps and two small jumps;
% f5: a function has two large jumps and four small jumps;
% f5: a function has two large jumps and four small jumps;
% f6: a piecewise smooth function has more jumps(both large jumps and
% small jumps.

% f_num = no. of function;
% n     = no. of points on the function
% mres  = change resolution of Fourier coefficient calculation.
% x     = n+1 uniform points on [-1, 1]
% quadx = n*mres uniform points on [-1, 1]

function [fx, f_hat] = f_example_function(f_num,n,mres,x,quadx)

fx = zeros(n+1,1);
pp=10; 
qq=1;
switch f_num
    case 'f1'
        for j=1:n+1
            fx(j)=sin(pi*x(j));
        end
    case 'f2'
        for j=1:n+1
            if x(j)<0
                fx(j)=-x(j)-1;
            else if x(j)>=0
                    fx(j)=-x(j)+1;
                end
            end
        end       
    case 'f3'    
        for j=1:n+1
            if (x(j)> -3/5) && (x(j)<=-1/5)
                fx(j)= 1;
            elseif (x(j)> 1/5) && (x(j)<=3/5)
                fx(j) = 0.9;
            else
                fx(j) = 0;
            end
        end
     case 'f4'
         for j=1:n+1
            if (x(j)> -3/5) && (x(j)<=-1/5)
                fx(j)= 1;
            elseif (x(j)> 1/5) && (x(j)<=3/5)
                fx(j) = 0.1;
            else
                fx(j) = 0;
            end
        end
     case 'f5'
        for j=1:n+1
            if x(j)<-3/4
                fx(j)=0;
            elseif -3/4<=x(j) && x(j)<-1/2
                    fx(j)=3/2*pp;
            elseif -1/4<=x(j) && x(j)<1/8
                 fx(j)=7/4-x(j)*pi/(2)+sin(x(j)*pi-1/4);
            elseif 3/8<=x(j) && x(j)<3/4
                 fx(j)=qq*(11*x(j)*pi/(4)-5);
            else 
                 fx(j)=0;
            end
        end
    case 'f6'
        for j=1:n+1
            if -1<=x(j) && x(j)<-5/6
                 fx(j)=0;
            elseif -5/6<=x(j) && x(j)<-1/2
                 fx(j)=2*exp(3*sin(2.5*pi*(x(j)+0.5)));
%                 elseif  -1/4<=x(j) && x(j)<0
            elseif  -7/20<=x(j) && x(j)<0
                 fx(j)=exp(3*sin(2.5*pi*(x(j)+0.5)))/4;
            elseif 2/5<=x(j) && x(j)<5/6
                 fx(j)=sin(pi*x(j))*(11*x(j)*x(j)*pi/4-1);
            elseif 5/6<=x(j) && x(j)<=1
                 fx(j)=-4*sin(pi*x(j));
            else 
                 fx(j)=1-x(j);
            end
        end 
end

%% Test function defined on mres*n grid points for numerical integration
f_quadx=zeros(mres*n,1);   
h=quadx(2)-quadx(1);
switch f_num
    case 'f1'
        for j=1:mres*n
%             if quadx(j)<0
                f_quadx(j)=sin(pi*quadx(j));
%             else if quadx(j)>=0
%                     f_quadx(j)=exp(-3*pi*quadx(j));
%                 end
%             end
        end
    case 'f2'
        for j=1:mres*n
            if quadx(j)<0
                f_quadx(j)=-quadx(j)-1;
            else if quadx(j)>=0
                    f_quadx(j)=-quadx(j)+1;
                end
            end
        end
    case 'f3'    
        for j=1:mres*n
            if (quadx(j)> -3/5) && (quadx(j)<=-1/5)
                f_quadx(j)= 1;
            elseif (quadx(j)> 1/5) && (quadx(j)<=3/5)
                f_quadx(j) = 0.9;
            else
                f_quadx(j) = 0;
            end
        end
     case 'f4'
         for j=1:mres*n
            if (quadx(j)> -3/5) && (quadx(j)<=-1/5)
                f_quadx(j)= 1;
            elseif (quadx(j)> 1/5) && (quadx(j)<=3/5)
                f_quadx(j) = 0.1;
            else
                f_quadx(j) = 0;
            end
        end
    case 'f5'
        for j=1:mres*n
            if quadx(j)<-3/4
                f_quadx(j)=0;
            elseif -3/4<=quadx(j) && quadx(j)<-1/2
                 f_quadx(j)=3/2*pp;
            elseif -1/4<=quadx(j) && quadx(j)<1/8
                 f_quadx(j)=7/4-quadx(j)*pi/(2)+sin(quadx(j)*pi-1/4);
            elseif 3/8<=quadx(j) && quadx(j)<3/4
                 f_quadx(j)=qq*(11*quadx(j)*pi/(4)-5);
            else 
                 f_quadx(j)=0;
            end
        end
    case 'f6'
        for j=1:mres*n
            if -1<=quadx(j) && quadx(j)<-5/6
                 f_quadx(j)=0;
            elseif -5/6<=quadx(j) && quadx(j)<-1/2
                 f_quadx(j)=2*exp(3*sin(2.5*pi*(quadx(j)+0.5)));
%                 elseif -1/4<=quadx(j) && quadx(j)<0 
            elseif-7/20<=quadx(j) && quadx(j)<0
                 f_quadx(j)=exp(3*sin(2.5*pi*(quadx(j)+0.5)))/4;
            elseif 2/5<=quadx(j) && quadx(j)<5/6
                 f_quadx(j)=sin(pi*quadx(j))*(11*quadx(j)*quadx(j)*pi/4-1);
            elseif 5/6<=quadx(j) && quadx(j)<=1
                 f_quadx(j)=-4*sin(pi*quadx(j));
            else 
                 f_quadx(j)=1-quadx(j);
            end
        end
end

%% Calculate Fourier Coefficients using quadrature formula
k = -n/2:n/2; k=k';  
f_hat=zeros(n+1,1);

for i=1:n+1
    for j=1:mres*n
        f_hat(i)=f_hat(i)+f_quadx(j)*exp(-1i*k(i)*pi*quadx(j))*h;
    end
end
f_hat=f_hat/2;

return




