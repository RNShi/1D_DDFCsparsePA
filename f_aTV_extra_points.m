% function: Using average total variation to find the extra points
% n_more_l = no. of extra points on the left
% n_more_r = no. of extra points on the right


function [n_more_l, n_more_r] = f_aTV_extra_points(fn_partial, r, beta, n1, n2)
%right
n_more_l=zeros(r,1);
n_more_r=zeros(r,1);
K=1.5;
n_more_max=2*beta;

for i=1:r-1
    nTV_r(i)=abs(fn_partial(2,i+1)-fn_partial(1,i+1))+abs(fn_partial(3,i+1)-fn_partial(2,i+1));
    n_more_r(i)=2;
    for j=3:n_more_max
        nTV_r(i)=nTV_r(i)+abs(fn_partial(j+1,i+1)-fn_partial(j,i+1));
        if(abs(fn_partial(j+1,i+1)-fn_partial(j,i+1))>nTV_r(i)/(n_more_r(i)+1)*K)
            break;
        end
        n_more_r(i)=n_more_r(i)+1;
    end
end
nTV_r(r)=abs(fn_partial(2,1)-fn_partial(1,1))+abs(fn_partial(3,1)-fn_partial(2,1));
n_more_r(r)=2;
for j=3:n_more_max
    nTV_r(r)=nTV_r(r)+abs(fn_partial(j+1,1)-fn_partial(j,1));
    if(abs(fn_partial(j+1,1)-fn_partial(j,1))>nTV_r(r)/(n_more_r(r)+1)*K)
        break;
    end
    n_more_r(r)=n_more_r(r)+1;
end

%left
for i=2:r
    nTV_l(i)=abs(fn_partial(n1+1,i-1)-fn_partial(n1,i-1))+abs(fn_partial(n1,i-1)-fn_partial(n1-1,i-1));
    n_more_l(i)=2;
    for j=1:n_more_max-2
        nTV_l(i)=nTV_l(i)+abs(fn_partial(n1-j,i-1)-fn_partial(n1-j-1,i-1));
        if(abs(fn_partial(n1-j,i-1)-fn_partial(n1-j-1,i-1))>nTV_l(i)/(n_more_l(i)+1)*K)
            break;
        end
        n_more_l(i)=n_more_l(i)+1;
    end
end
nTV_l(1)=abs(fn_partial(n2+1,r)-fn_partial(n2,r))+abs(fn_partial(n2,r)-fn_partial(n2-1,r));
n_more_l(1)=2;
for j=1:n_more_max-2
    nTV_l(1)=nTV_l(1)+abs(fn_partial(n2-j,r)-fn_partial(n2-j-1,r));
    if(abs(fn_partial(n2-j,r)-fn_partial(n2-j-1,r))>nTV_l(1)/(n_more_l(1)+1)*K)
        break;
    end
    n_more_l(1)=n_more_l(1)+1;
end

return
