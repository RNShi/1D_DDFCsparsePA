% function: Fourier coefficients
% input: a periodic function
% output: its Fourier coefficients

function fc_hat = f_fourier_coefficients(fc)

N=size(fc,1);
if rem(N,2)==0  %even expansion:N1=evem number
    j=0:N;
    z1=2*pi*j/N;
    l1=-N/2:N/2;l1=l1';
    fc_hat=zeros(N+1,1);
    for i=1:N+1
        for j=1:N
            fc_hat(i,1)=fc_hat(i,1)+fc(j)*exp(-1i*l1(i)*z1(j));
        end
    end
    fc_hat(1,1)=fc_hat(1,1)/(2*N);
    fc_hat(N+1,1)=fc_hat(N+1,1)/(2*N);
    fc_hat(2:N,1)=fc_hat(2:N,1)/N; 
else        %odd expansion:new N1=odd number
    j=0:N-1;
    z1=2*pi*j/(N-1);
    l1=-(N-1)/2:(N-1)/2;l1=l1';
    fc_hat=zeros(N,1);
    for i=1:N
        for j=1:N
            fc_hat(i,1)=fc_hat(i,1)+fc(j)*exp(-1i*l1(i)*z1(j));
        end
    end
    fc_hat(:,1)=fc_hat(:,1)/N;
end

return