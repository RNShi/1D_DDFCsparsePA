% function: discrete Fourier transform matrix
% input: size N
% output: discrete Fourier transform matrix

function A = f_DFT(N)
if rem(N,2)==0  %even expansion:N1=evem number
    j=0:N;
    z1=2*pi*j/N;
    l1=-N/2:N/2;l1=l1';
    A=(1/N)*exp(-1i*l1*z1);
    A(1,:)=A(1,:)/2;
    A(N+1,:)=A(N+1,:)/2;
else        %odd expansion:new N1=odd number
    j=0:N-1;
    z1=2*pi*j/(N-1);
    l1=-(N-1)/2:(N-1)/2;l1=l1';
    A=(1/N)*exp(-1i*l1*z1);
end
return