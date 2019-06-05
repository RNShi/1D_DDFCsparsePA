% function: Using advanced average total variation to find the extra points
% n_more_l_new = no. of extra points on the left
% n_more_r_new = no. of extra points on the right

function [n_more_l_new, n_more_r_new] = f_advanced_aTV_extra_points(f_partial, n_more_l, n_more_r, r, beta, n, n1)
%right
n_more_min=2;
n_more_l_new=zeros(r,1);
n_more_r_new=zeros(r,1);
nTV_r_min=zeros(r,1);
nTV_r_temp=zeros(r,1);
nTV_l_min=zeros(r,1);
nTV_l_temp=zeros(r,1);
for i=1:r-1
    n_more_r_new(i)=n_more_min;
    for j=1:beta
        nTV_r_min(i)= nTV_r_min(i)+abs(f_partial(i*n1-beta+j+1+n_more_min)-f_partial(i*n1-beta+j+n_more_min));
    end
    for j=1:n_more_r(i)-n_more_min
        nTV_r_temp(i)=0;
        for l=1:beta
            nTV_r_temp(i)= nTV_r_temp(i)+abs(f_partial(i*n1-beta+l-j+1+n_more_min)-f_partial(i*n1-beta+l-j+n_more_min));
        end
        if nTV_r_temp(i)<=nTV_r_min(i)
            nTV_r_min(i)=nTV_r_temp(i);
            n_more_r_new(i)=j+n_more_min;
        end
    end
end
temp_v=zeros(1,n_more_r(r)+beta+1);
for i=1:beta-n_more_min
    temp_v(i)=f_partial(n-beta+i+n_more_min);
end
for i=beta-n_more_min+1:n_more_r(r)+beta+1
    temp_v(i)=f_partial(i-(beta-n_more_min));
end
for j=1:beta
    nTV_r_min(r)= nTV_r_min(r)+abs(temp_v(j+1)-temp_v(j));
end
for j=1:n_more_r(r)-n_more_min
    n_more_r_new(r)=n_more_min;
    nTV_r_temp(r)=0;
    for l=1:beta
        nTV_r_temp(r)= nTV_r_temp(r)+abs(temp_v(l+j+1)-temp_v(l+j));
    end
    if nTV_r_temp(r)<=nTV_r_min(r)
        nTV_r_min(r)=nTV_r_temp(r);
        n_more_r_new(r)=j+n_more_min;
    end
end

%left
for i=2:r
    n_more_l_new(i)=n_more_min;
    for j=1:beta
        nTV_l_min(i)= nTV_l_min(i)+abs(f_partial((i-1)*n1+j+1-n_more_min)-f_partial((i-1)*n1+j-n_more_min));
    end
    for j=1:n_more_l(i)-n_more_min
        nTV_l_temp(i)=0;
        for l=1:beta
            nTV_l_temp(i)= nTV_l_temp(i)+abs(f_partial((i-1)*n1+l-j+1-n_more_min)-f_partial((i-1)*n1+l-j-n_more_min));
        end
        if nTV_l_temp(i)<=nTV_l_min(i)
            nTV_l_min(i)=nTV_l_temp(i);
            n_more_l_new(i)=j+n_more_min;
        end
    end
end
temp_v=zeros(1,n_more_l(1)+beta+1);
for i=1:n_more_l(1)+n_more_min
    temp_v(i)=f_partial(n-(n_more_l(1)+n_more_min)+i);
end
for i=n_more_l(1)+n_more_min+1:n_more_l(1)+beta+1
    temp_v(i)=f_partial(i-(n_more_l(1)+n_more_min));
end
for j=1:beta
    nTV_l_min(1)= nTV_l_min(1)+abs(temp_v(n_more_l(1)+beta+1+1-j)-temp_v(n_more_l(1)+beta+1-j));
end
for j=1:n_more_l(1)-n_more_min
    n_more_l_new(1)=n_more_min;
    nTV_l_temp(1)=0;
    for l=1:beta
        nTV_l_temp(1)= nTV_l_temp(1)+abs(temp_v(n_more_l(1)+l-j+1)-temp_v(n_more_l(1)+l-j));
    end
    if nTV_l_temp(1)<=nTV_l_min(1)
        nTV_l_min(1)=nTV_l_temp(1);
        n_more_l_new(1)=j+n_more_min;
    end
end

return
