% function: solving the convex optimization problem using sparse PA L1
% regularization

function f = f_L1_regularization(fc_hat, A, lambda)
N = size(fc_hat,1);
m_1=2;

C_sub = ones(m_1+1,1);
for j = 1:m_1+1
    for js =1:m_1+1
        if js ~= j
            C_sub(j) = C_sub(j)/(j-js);
        end
    end
end
C_sub = C_sub/sum(C_sub(1:floor(m_1/2)+1)); % Normalize Coefficients

%% edge detection matrix
if rem(N,2)==0 
    C_new= zeros(N,N+1); 
    for j=1:N
        inds = j-floor(m_1/2):j+m_1-floor(m_1/2);
        if min(inds)<1
            gap=1-min(inds);
            inds = inds+gap;  
        elseif max(inds)>N+1
            gap=max(inds)-(N+1);
            inds = inds-gap;
        end
        C_new(j,inds) = C_sub;
    end
else
    C_new= zeros(N,N); 
    for j=1:N
        inds = j-floor(m_1/2):j+m_1-floor(m_1/2);
        if min(inds)<1
            gap=1-min(inds);
            inds = inds+gap;  
        elseif max(inds)>N
            gap=max(inds)-N;
            inds = inds-gap;
        end
        C_new(j,inds) = C_sub;
    end
end

%% Begin CVX Program
if rem(N,2)==0
    cvx_begin quiet
        variable f(N+1);
        minimize(norm(A*f-fc_hat,2)+ lambda*norm(C_new*f,1));  
    cvx_end
else
    cvx_begin quiet
        variable f(N);
        minimize(norm(A*f-fc_hat,2)+ lambda*norm(C_new*f,1));  
    cvx_end
end

return