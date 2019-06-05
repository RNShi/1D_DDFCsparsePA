clear all 
% close all
format long

%% input data
f_num = 'f6';
% no. of input Fourier modes is n+1
n = 512; 
%no. of seperate intevals
r = 3; 
method = input('Please type the number of method (1-non-overlapping; 2-static overlapping; 3-overlapping \n with average total variation; 4-overlapping with advanced average total variation): ');

%% Creates Uniform Grid for Reconstruction
lambda = 0.015;
mres = 12;  %% Change Resolution of Fourier Coefficient Calculation.
j = 0:n;
x = -1+2*j/n;
dx = x(2)-x(1);
fx = zeros(n+1,1); 
j = 0:mres*n-1;
quadx = -1+j/(mres*n/2);
% f_quadx = zeros(mres*n,1); 
k = -n/2:n/2; k=k'; 

%% Uses example_function.m to get example function and the 
%% correspond Fourier coefficient
f_hat = zeros(n+1,1);
[fx,f_hat] = f_example_function(f_num,n,mres,x,quadx);

%% Fourier partial sum
% (with Gibbs oscillations if function is piecewise smooth)
S_n = zeros(n,n+1); %The inverse discrete Fourier transform matrix
for i = 1:n
    for j = 1:n+1
        S_n(i,j) = exp(1i*k(j)*pi*x(i));
    end
end
f_partial = S_n*f_hat;

%% Filter Fourier Reconstruction (for comparison)
alpha = 32;
p = 5;
G = zeros(n,n+1);
for i = 1:n+1
    for j = 1:n
    G(j,i) = exp(-alpha*(abs(k(i)/(n/2))).^p);
    end
end
F = G.*S_n;
f_filter = F*f_hat;
f_filter(n+1,1) = f_filter(1,1);

%% some constant
n1 = ceil(n/r);   % no. of points in first r-1 intervals
re = r*n1-n;
n2 = n1-re;     % no. of points in the last intervals
% beta = the size of the data near the boundaries used for the FC 
beta = 10;
% gamma = the size of the extension data for the FC 
gamma = 20;

N1 = n1+gamma;
N2 = n2+gamma;

%% non-overlapping domain decomposition 
% split f_partial to r non-overlapping parts
fn_partial = zeros(n1+1,r);  
for i = 1:r-1
    for j = 1:n1+1
        fn_partial(j,i) = f_partial(j+(i-1)*n1);
    end
end
for j = 1:n2
    fn_partial(j,r) = f_partial(j+(r-1)*n1);
end
fn_partial(n2+1,r) = f_partial(1);

%% %%%%%%%%%%%%%%%%%%%%%%begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%non-overlapping%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if method == 1

%%  Finds Fourier continuation (periodic)
fc = zeros(N1,r);
for i = 1:r-1
    fc(:,i) = f_FC_method(fn_partial(:,i),beta,gamma);
end
fc(1:N2,r) = f_FC_method(fn_partial(1:n2+1,r),beta,gamma);

%% Fourier Coefficients for periodic extension
if rem(N1,2) == 0
    fc_hat = zeros(N1+1,r);
else
    fc_hat = zeros(N1,r);
end

for i = 1:r-1
    fc_hat(:,i) = f_fourier_coefficients(fc(:,i));
end

if rem(N2,2) == 0
    fc_hat(1:N2+1,r) = f_fourier_coefficients(fc(1:N2,r));
else
    fc_hat(1:N2,r) = f_fourier_coefficients(fc(1:N2,r));
end

%% discrete Fourier transform matrix
A1 = f_DFT(N1);
A2 = f_DFT(N2);

%% Uses L1 regularization with sparse PA to find the reconstruction
f1 = 0*fx;
if rem(N1,2) == 0
    frecon = zeros(N1+1,r);
else
    frecon = zeros(N1,r);
end
for i = 1:r-1
    frecon(:,i) = f_L1_regularization(fc_hat(:,i), A1, lambda);
    if i == 1
        f1(1,1) = frecon(1,1);
    else
        f1(1+(i-1)*n1,1) = (f1(1+(i-1)*n1,1)+frecon(1,i))/2;
    end
    for j = 2:n1+1
        f1(j+(i-1)*n1,1) = frecon(j,i);
    end
end

if rem(N2,2) == 0 
    frecon(1:N2+1,r) = f_L1_regularization(fc_hat(1:N2+1,r), A2, lambda);
else
    frecon(1:N2,r) = f_L1_regularization(fc_hat(1:N2,r), A2, lambda);
end
f1(1+(r-1)*n1,1) = (f1(1+(r-1)*n1,1)+frecon(1,r))/2;
for i = 2:n2
    f1(i+(r-1)*n1,1) = frecon(i,r);
end
f1(n+1,1) = frecon(n1+1-re+1,r);

end
%%%%%%%%%%%%%%%%%%%%%% non-overlapping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %%%%%%%%%%%%%%%%%%%%%%%% begin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% static overlapping%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if method == 2

n_more = 10;
n1_new = n1+2*n_more;
n2_new = n2+2*n_more;
N1 = n1_new+gamma;
N2 = n2_new+gamma;

%% static overlapping domain decomposition
% split f_partial to r static overlapping parts
fn_partial_new = zeros(n1_new+1,r); 
fn_partial_new(1:n_more,1) = f_partial(n-n_more+1:n);
fn_partial_new(n_more+1:n1_new+1,1) = f_partial(1:n1_new+1-n_more);
for i = 2:r-1
    for j = 1:n1_new+1
        fn_partial_new(j,i) = f_partial(j+(i-1)*n1-n_more);
    end
end
for j = 1:n2_new-n_more
    fn_partial_new(j,r) = f_partial(j+(r-1)*n1-n_more);
end
fn_partial_new(n2_new-n_more+1:n2_new+1,r) = f_partial(1:n_more+1);

%%  Finds Fourier continuation (periodic)
fc = zeros(N1,r);
for i = 1:r-1
    fc(:,i) = f_FC_method(fn_partial_new(:,i),beta,gamma);
end
fc(1:N2,r) = f_FC_method(fn_partial_new(1:n2_new+1,r),beta,gamma);

%% Fourier Coefficients for periodic extension
if rem(N1,2) == 0
    fc_hat = zeros(N1+1,r);
else
    fc_hat = zeros(N1,r);
end

for i = 1:r-1
    fc_hat(:,i) = f_fourier_coefficients(fc(:,i));
end

if rem(N2,2) == 0
    fc_hat(1:N2+1,r) = f_fourier_coefficients(fc(1:N2,r));
else
    fc_hat(1:N2,r) = f_fourier_coefficients(fc(1:N2,r));
end
%% discrete Fourier transform matrix
A1 = f_DFT(N1);
A2 = f_DFT(N2);

%% Uses L1 regularization with sparse PA to find the reconstruction
f1 = 0*fx;
if rem(N1,2) == 0
    frecon = zeros(N1+1,r);
else
    frecon = zeros(N1,r);
end
for i = 1:r-1
    frecon(:,i) = f_L1_regularization(fc_hat(:,i), A1, lambda);
    if i == 1
        f1(1,1) = frecon(i+n_more,1);
    else
        f1(1+(i-1)*n1,1) = (f1(1+(i-1)*n1,1)+frecon(1+n_more,i))/2;
    end
    for j = 2:n1+1
        f1(j+(i-1)*n1,1) = frecon(j+n_more,i);
    end
end

if rem(N2,2) == 0 
    frecon(1:N2+1,r) = f_L1_regularization(fc_hat(1:N2+1,r), A2, lambda);
else
    frecon(1:N2,r) = f_L1_regularization(fc_hat(1:N2,r), A2, lambda);
end
f1(1+(r-1)*n1,1) = (f1(1+(r-1)*n1,1)+frecon(1+n_more,r))/2;
for i = 2:n2
    f1(i+(r-1)*n1,1) = frecon(i+n_more,r);
end
f1(n+1,1) = frecon(n1+1-re+1+n_more,r);

f1(1,1) = (f1(1,1)+f1(n+1,1))/2;
f1(n+1,1) = f1(1,1);

end
%%%%%%%%%%%%%%%%%%%%% static overlapping %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %%%%%%%%%%%%%%%%%%%%%%%% begin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% overlapping using average total variation %%%%%%%%%%%%%%%

if method == 3
    
%% Uses average total variation to find the extra points

[n_more_l, n_more_r] = f_aTV_extra_points(fn_partial, r, beta, n1, n2);
n1_new = zeros(r,1);
for i = 1:r-1
    n1_new(i) = n1+n_more_l(i)+n_more_r(i);
end
n1_new(r) = n2+n_more_l(r)+n_more_r(r);

N1 = n1_new+gamma;

%% overlapping domain decomposition using aTV
% split f_partial to r overlapping parts using aTV
fn_partial_new = zeros(max(n1_new)+1,r);
fn_partial_new(1:n_more_l(1),1) = f_partial(n-n_more_l+1:n);
fn_partial_new(n_more_l(1)+1:n1_new(1)+1,1) = f_partial(1:n1_new(1)+1-n_more_l(1));
for i = 2:r-1
    for j = 1:n1_new(i)+1
        fn_partial_new(j,i) = f_partial(j+(i-1)*n1-n_more_l(i));
    end
end
for j = 1:n1_new(r)-n_more_r(r)
    fn_partial_new(j,r) = f_partial(j+(r-1)*n1-n_more_l(r));
end
fn_partial_new(n1_new(r)-n_more_r(r)+1:n1_new(r)+1,r) = f_partial(1:n_more_r(r)+1);

%% Finds Fourier continuation (periodic)
fc = zeros(max(N1),r);
for i = 1:r
    fc(1:N1(i),i) = f_FC_method(fn_partial_new(1:n1_new(i)+1,i),beta,gamma);
end

%% Fourier Coefficients for periodic extension
fc_hat = zeros(max(N1)+1,r);
for i = 1:r
    if rem(N1(i),2) == 0
        fc_hat(1:N1(i)+1,i) = f_fourier_coefficients(fc(1:N1(i),i));
    else
        fc_hat(1:N1(i),i) = f_fourier_coefficients(fc(1:N1(i),i));
    end
end

%% discrete Fourier transform matrix
A1 = zeros(max(N1)+1,max(N1),r);
for i = 1:r
    if rem(N1(i),2) == 0
        A1(1:N1(i)+1,1:N1(i)+1,i) = f_DFT(N1(i));
    else
        A1(1:N1(i),1:N1(i),i) = f_DFT(N1(i));
    end
end

%% Uses L1 regularization with sparse PA to find the reconstruction
f1 = 0*fx;
frecon = zeros(max(N1)+1,r);
for i = 1:r-1
    if rem(N1(i),2) == 0 
        A = A1(1:N1(i)+1,1:N1(i)+1,i);
        frecon(1:N1(i)+1,i) = f_L1_regularization(fc_hat(1:N1(i)+1,i), A, lambda);
    else
        A = A1(1:N1(i),1:N1(i),i);
        frecon(1:N1(i),i) = f_L1_regularization(fc_hat(1:N1(i),i), A, lambda);
    end
    if i == 1
        f1(1,1) = frecon(i+n_more_l(i),1);
    else
        f1(1+(i-1)*n1,1) = (f1(1+(i-1)*n1,1)+frecon(1+n_more_l(i),i))/2;
    end
    for j = 2:n1+1
        f1(j+(i-1)*n1,1) = frecon(j+n_more_l(i),i);
    end
end

if rem(N1(r),2) == 0 
    A = A1(1:N1(r)+1,1:N1(r)+1,r);
    frecon(1:N1(r)+1,r) = f_L1_regularization(fc_hat(1:N1(r)+1,r), A, lambda);
else
    A = A1(1:N1(r),1:N1(r),r);
    frecon(1:N1(r),r) = f_L1_regularization(fc_hat(1:N1(r),r), A, lambda);
end
f1(1+(r-1)*n1,1) = (f1(1+(r-1)*n1,1)+frecon(1+n_more_l(r),r))/2;
for i = 2:n2
    f1(i+(r-1)*n1,1) = frecon(i+n_more_l(r),r);
end
f1(n+1,1) = frecon(n1+1-re+1+n_more_l(r),r);

f1(1,1) = (f1(1,1)+f1(n+1,1))/2;
f1(n+1,1) = f1(1,1);

end
%%%%%%%%%%%% overlapping using average total variation %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% overlapping using average total variation (advanced)%%%%%%%%%
if method == 4

%% Uses advanced average total variation to find the extra points
[n_more_l, n_more_r] = f_aTV_extra_points(fn_partial, r, beta, n1, n2);
[n_more_l_new, n_more_r_new] = f_advanced_aTV_extra_points(f_partial, n_more_l, n_more_r, r, beta, n, n1);

n1_new = zeros(r,1);
for i = 1:r-1
    n1_new(i) = n1+n_more_l_new(i)+n_more_r_new(i);
end
n1_new(r) = n2+n_more_l_new(r)+n_more_r_new(r);

N1 = n1_new+gamma;

%% overlapping domain decomposition using aTV
% split f_partial to r overlapping parts using aTV
fn_partial_new = zeros(max(n1_new)+1,r);
fn_partial_new(1:n_more_l_new(1),1) = f_partial(n-n_more_l_new+1:n);
fn_partial_new(n_more_l_new(1)+1:n1_new(1)+1,1) = f_partial(1:n1_new(1)+1-n_more_l_new(1));
for i = 2:r-1
    for j = 1:n1_new(i)+1
        fn_partial_new(j,i) = f_partial(j+(i-1)*n1-n_more_l_new(i));
    end
end
for j = 1:n1_new(r)-n_more_r_new(r)
    fn_partial_new(j,r) = f_partial(j+(r-1)*n1-n_more_l_new(r));
end
fn_partial_new(n1_new(r)-n_more_r_new(r)+1:n1_new(r)+1,r) = f_partial(1:n_more_r_new(r)+1);

%% Finds Fourier continuation (periodic)
fc = zeros(max(N1),r);
for i = 1:r
    fc(1:N1(i),i) = f_FC_method(fn_partial_new(1:n1_new(i)+1,i),beta,gamma);
end

%% Fourier Coefficients for periodic extension
fc_hat = zeros(max(N1)+1,r);
for i = 1:r
    if rem(N1(i),2) == 0
        fc_hat(1:N1(i)+1,i) = f_fourier_coefficients(fc(1:N1(i),i));
    else
        fc_hat(1:N1(i),i) = f_fourier_coefficients(fc(1:N1(i),i));
    end
end

%% discrete Fourier transform matrix
A1 = zeros(max(N1)+1,max(N1),r);
for i = 1:r
    if rem(N1(i),2) == 0
        A1(1:N1(i)+1,1:N1(i)+1,i) = f_DFT(N1(i));
    else
        A1(1:N1(i),1:N1(i),i) = f_DFT(N1(i));
    end
end

%% Uses L1 regularization with sparse PA to find the reconstruction
f1 = 0*fx;
frecon = zeros(max(N1)+1,r);
for i = 1:r-1
    if rem(N1(i),2) == 0 
        A = A1(1:N1(i)+1,1:N1(i)+1,i);
        frecon(1:N1(i)+1,i) = f_L1_regularization(fc_hat(1:N1(i)+1,i), A, lambda);
    else
        A = A1(1:N1(i),1:N1(i),i);
        frecon(1:N1(i),i) = f_L1_regularization(fc_hat(1:N1(i),i), A, lambda);
    end
    if i == 1
        f1(1,1) = frecon(i+n_more_l_new(i),1);
    else
        f1(1+(i-1)*n1,1) = (f1(1+(i-1)*n1,1)+frecon(1+n_more_l_new(i),i))/2;
    end
    for j = 2:n1+1
        f1(j+(i-1)*n1,1) = frecon(j+n_more_l_new(i),i);
    end
end

if rem(N1(r),2) == 0 
    A = A1(1:N1(r)+1,1:N1(r)+1,r);
    frecon(1:N1(r)+1,r) = f_L1_regularization(fc_hat(1:N1(r)+1,r), A, lambda);
else
    A = A1(1:N1(r),1:N1(r),r);
    frecon(1:N1(r),r) = f_L1_regularization(fc_hat(1:N1(r),r), A, lambda);
end
f1(1+(r-1)*n1,1) = (f1(1+(r-1)*n1,1)+frecon(1+n_more_l_new(r),r))/2;
for i = 2:n2
    f1(i+(r-1)*n1,1) = frecon(i+n_more_l_new(r),r);
end
f1(n+1,1) = frecon(n1+1-re+1+n_more_l_new(r),r);

f1(1,1) = (f1(1,1)+f1(n+1,1))/2;
f1(n+1,1) = f1(1,1);

end
%%%%%%%%%%%% overlapping using average total variation (advanced)%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Error Analysis
f_exact = fx;
f_partial = S_n*f_hat;
f_partial(n+1,1) = f_partial(1,1);
f_DDFCsparsePA = f1;

err_partial = abs(f_exact-f_partial);
err_filter = abs(f_exact-f_filter);
err_DDFCsparsePA = abs(f_exact-f_DDFCsparsePA);

%% Reconstruction Plots
figure, h = plot(x,fx, '-b',x,real(f_partial),'or');
set(h(1),'LineWidth',1.5),set(h(2),'MarkerSize',4)
legend('Original Function', 'Fourier Partial Sum Reconstruction')
set(gca,'FontSize',13);
xlabel('x','fontsize',13),ylabel('y','fontsize',13)

figure, h = plot(x,fx, '-b',x,real(f_filter),'or');
set(h(1),'LineWidth',1.5),set(h(2),'MarkerSize',4)
legend('Original Function', 'filtered Fourier Partial Sum Reconstruction')
set(gca,'FontSize',13);
xlabel('x','fontsize',13),ylabel('y','fontsize',13)

figure, h = plot(x,fx,'-b',x,real(f_DDFCsparsePA),'ro');
set(h(1),'LineWidth',1.5),set(h(2),'MarkerSize',4)
legend('Original Function','DDFCsparsePA Reconstruction')
set(gca,'FontSize',13);
xlabel('x','fontsize',13),ylabel('y','fontsize',13)

%% Error plots
figure, h = plot(x,log(err_partial),'ksquare-',...
                 x,log(err_filter),'-ro',...
                 x,log(err_DDFCsparsePA),'-bo',...
                 'LineWidth',1,'MarkerSize',4);
xlabel('x','FontSize',14),ylabel('log Error','FontSize',14)
set(gca,'FontSize',14)

legend('Fourier Partial Sum Reconstruction', 'filtered Fourier Partial Sum Reconstruction',...
    'DDFCsparsePA Reconstruction')
