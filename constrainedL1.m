function [x,p] = constrainedL1(A,b,pm)
%initialize variables for inner problem
[M,N] = size(A);

%% parameters
lambda = 2; detla = 10*lambda;
maxit = 2*N; tol = 1e-3; x0 = zeros(N,1);
abstol = 1e-7; reltol = 1e-5;
eps = 1e-16;


if isfield(pm,'delta'); delta = pm.delta; end
if isfield(pm,'lambda'); lambda = pm.lambda; end
if isfield(pm,'maxit'); maxit = pm.maxit; end
if isfield(pm,'tol'); tol = pm.tol; end
if isfield(pm,'x0'); x0 = pm.x0; end
if isfield(pm,'abstol'); abstol = pm.abstol; end
if isfield(pm,'reltol'); reltol = pm.reltol; end


%% pre-computing/initialize
AAt = A*A';
L = chol( speye(M) + 1/delta*AAt, 'lower' );
L = sparse(L);
U = sparse(L');

x = x0;
Atb = A'*b;
y = zeros(N,1); u = y;
v = b; 



for it = 1:maxit
    
     %update x
     rhs = A'*(v+b) + delta*(y+u);
     x = rhs/delta - (A'*(U\(L\(A*rhs))))/delta^2;
 
     %update y
     yold = y;
     y =shrink(x-u, lambda/delta);
  
     %update u
     u = u + y -x;
    
     %update v
     v = v + b - A*x;
        
     % stopping condition for ADMM
     r = norm(x-y);
     s = norm(delta*(y-yold));
            
     eps_pri = sqrt(N)*abstol + reltol*max(norm(x),norm(y));
     eps_dual = sqrt(N)*abstol + reltol*norm(u);
            
     if (r < eps_pri && s < eps_dual)
         break;
     end    
end

end

function z = shrink(x,r)
    z = max(0,x - r) - max(0,-x - r);
end