% function: 1D Fourier continuation method using boundary data
% input: non-periodic function
% output: periodic function on the extended region

% f_input = non-periodic function
% beta    = the size of the data near the boundaries used for the FC 
% gamma   = the size of the extension data for the FC 

function fc = f_FC_method(f_input,beta,gamma)

N = size(f_input,1); % size of input data
x0 = 0; x1 = 1;

% -------------------------- some parameters ---------------------------     
if beta>5
    m = 5; % % m=M<beta
else
    m = beta-1;
end
Q     = 100;     
mm    = 30 ;     % mm = g
dx    = (x1-x0)/(N-1);  % dx=dx
delta = beta*dx;     % delta=beta*dx
d     = (gamma-1)*dx;   % d=gamma*dx
d_Q   = delta/Q;    % d_Q=delta/Q
x = zeros(N,1);
for ii = 1:N
    x(ii) = x0 + (ii-1)*dx;    
end
sN = N+gamma-1; % size of the periodic function reconstructed by FC
bx = (x1-x0)/(N-1)*(sN-1); 

%-------------------- The points near the boundary -------------------
x_left (1:beta+1) = x(N-beta:N);
x_right(1:beta+1) = x(1:beta+1);

%----------------------- recurrence Gram polynomial ---------------------
a0 = x0;  a1 = x0+delta;
xi = 2*(x_right-a0)/(a1-a0) - 1;   %xi in [-1,1] beta+1 points

alpha = zeros(m,1);
for n = 1: m
    alpha(n) = beta/n*sqrt((4*n^2-1)/((beta+1)^2-n^2));
end  
P = zeros(beta+1, m+1);  
P(:,1) = 1/sqrt(beta+1);  P(:,2) = alpha(1)*xi'.*P(:,1);
for ii = 3: m+1
    P(:,ii) = alpha(ii-1)*xi'.*P(:,ii-1) - alpha(ii-1)/alpha(ii-2)*P(:,ii-2);
end

%------------------- function values near the boundaries ----------------
f_left  = f_input(N-beta:N,1);
f_right = f_input(1:beta+1,1);
a_left = zeros(m+1,1);
a_right = zeros(m+1,1);
for ii = 1:m+1
    a_left(ii,1)  = sum(f_left(:,1) .*(P(:,ii))); %a_left=b_l
    a_right(ii,1) = sum(f_right(:,1).*(P(:,ii)));
end

x_d = zeros(2*(beta+gamma-1),1);
for ii = 1:2*(beta+gamma-1)
    x_d(ii) = (x1-delta) + (ii-1)*dx;   % x in [b-delta,b+2*d+delta]  2*(beta+gamma) points
end
x_d_left = zeros(Q+1,1);
for ii = 1:Q+1
    x_d_left(ii) = (x1-delta) + (ii-1)*d_Q;    % x in [b-delta,b] Q+1  points
end

% ---------------- recurrence Gram polynomial for Q+1 points---------------
a0 = x1 - delta;  a1 = x1;
xi_Q = 2*(x_d_left-a0)/(a1-a0) - 1;   % xi_Q in [-1,1] Q+1 points

P_Q = zeros(Q+1, m+1);
P_Q(:,1) = 1/sqrt(beta+1); 
P_Q(:,2) = alpha(1)*xi_Q.*P_Q(:,1);
for ii = 3: m+1
    P_Q(:,ii) = alpha(ii-1)*xi_Q.*P_Q(:,ii-1) - alpha(ii-1)/alpha(ii-2)*P_Q(:,ii-2);
end

%------------------------ Get the even and odd k --------------------------
if mod(mm, 2) == 0
    t_g = -mm/2+1:mm/2;   %t_g=g_mm
else
    t_g = -(mm-1)/2:(mm-1)/2;
end

if mod(t_g(1), 2) == 0
    t_g_even = t_g(1:2:mm);
    t_g_odd  = t_g(2:2:mm);
else
    t_g_even = t_g(2:2:mm);
    t_g_odd  = t_g(1:2:mm);
end

%---------------- matrix A/B (for calculating a_hat/b_hat) ----------------
for ii = 1:Q+1
    angle = pi/(d+delta)*x_d_left(ii);
    for jj = 1:length(t_g_even)
        A(ii,jj) = exp(1i*angle*t_g_even(jj));  
    end
    for jj = 1:length(t_g_odd)
        B(ii,jj) = exp(1i*angle*t_g_odd(jj));
    end
end

[U_even, S_even, V_even] = svd(A,0);
[U_odd,  S_odd,  V_odd ] = svd(B,0);
S_even_inv = inv(S_even);
S_odd_inv = inv(S_odd);
svd_tol = 1e-11; 
for ii = 1:length(t_g_even)
    if(abs(S_even_inv(ii, ii)) > 1/svd_tol)
        S_even_inv(ii, ii) = 0;
    end
end
for ii = 1:length(t_g_odd)
    if(abs(S_odd_inv(ii, ii)) > 1/svd_tol)
        S_odd_inv(ii, ii) = 0;
    end
end
a_hat = zeros(length(t_g_odd), m+1);
b_hat = zeros(length(t_g_odd), m+1);
for ii=1:m+1
    a_hat(:,ii) = V_even*S_even_inv*(U_even')*P_Q(:,ii);
    b_hat(:,ii) = V_odd *S_odd_inv *(U_odd' )*P_Q(:,ii);
end
a_hat = conj(a_hat);
b_hat = conj(b_hat);

% -------------------- Compute f_even and f_odd -------------------------
f_even = zeros(gamma-2, m+1);
f_odd  = zeros(gamma-2, m+1);
for ii = 1:m+1
    for jj = 1:2*(beta+gamma-1)
        angle = pi/(d+delta)*x_d(jj);
        f_even(jj,ii) = sum(a_hat(:,ii)'.*exp(1i*angle*t_g_even));
        f_odd (jj,ii) = sum(b_hat(:,ii)'.*exp(1i*angle*t_g_odd ));
    end
end

%----------------------- find f_match function ----------------------
f_match = zeros(2*(beta+gamma-1),1);
for ii = 1:2*(beta+gamma-1)
    f_match(ii,1) = sum((a_left(:,1)+a_right(:,1))'.*f_even(ii,:) +...
                         (a_left(:,1)-a_right(:,1))'.*f_odd(ii,:))/2;
end

%---------------------- periodic function: fc ---------------------------
M = N+gamma-1-1;  
fc = zeros(M+1,1);
fc(1:N,1) = f_input(1:N,1);  
fc(N+1:M,1) = real(f_match(beta+2:beta+gamma-1,1));
fc(M+1,1) = f_input(1,1); 

return
