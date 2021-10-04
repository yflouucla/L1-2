function [x,output] = constrainedL1L2(A,b,pm)
%min_x |x|_1-|x|_2, s.t. Ax = b

%Input: dictionary A, data b, parameters set pm
%       pm.lambda: regularization paramter
%       pm.delta: penalty parameter for ADMM, default value: 10*lambda
%       pm.maxoit: max outer iterations, default value: 10
%       pm.maxit: max inner iterations: default value: 5000
%       pm.tol: outer tolerace, default value: 1e-3
%       pm.abstol: abs tolerance for ADMM: default value: 1e-7
%       pm.reltol: rel tolerance for ADMM: default value: 1e-5
%Output: computed coefficients x
%       output.relerr: relative error of x_n and x_{n-1}
%       output.obj: objective function of x_n

%% DCA for outer interation
[M,N] = size(A);

%% parameters
lambda = 2; detla = 10*lambda;
maxit = 2*N; maxoit = 10; 
tol = 1e-3; x0 = zeros(N,1);
abstol = 1e-7; reltol = 1e-5;
eps = 1e-16;


if isfield(pm,'delta'); delta = pm.delta; end
if isfield(pm,'lambda'); lambda = pm.lambda; end
if isfield(pm,'maxit'); maxit = pm.maxit; end
if isfield(pm,'maxoit'); maxoit = pm.maxoit; end
if isfield(pm,'x0'); x0 = pm.x0; end
if isfield(pm,'tol'); tol = pm.tol; end
if isfield(pm,'abstol'); abstol = pm.abstol; end
if isfield(pm,'reltol'); reltol = pm.reltol; end


%% pre-computing/initialize
AAt = A*A';
L = chol( speye(M) + 1/delta*AAt, 'lower' );
L = sparse(L);
U = sparse(L');
x = x0;
y = zeros(N,1); 
u = zeros(N,1); 
v = b; 

obj = @(x) (norm(x,1)-norm(x));
output.obj(1) = obj(x);


for oit = 1: maxoit
    
    c = x/(norm(x,2)+eps);
    xold = x;    
    
    %ADMM method for solving the sub-problem
    for it = 1:maxit
  
        %update x
        rhs = A'*(v+b) + lambda*c + delta*(y+u);
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
    
    % Stopping condition for DCA
    relerr = sqrt(sum((x-xold).^2))/max(sqrt(sum(x.^2)),1);
    
    output.relerr(oit) = relerr;
    output.obj(oit+1) = obj(x);
    
    if relerr < tol
        disp(['tolerance met after ' num2str(oit) ' iterations']);
        break;
    end
end

if (oit == maxoit)
    disp('Max outer iteration reached');
end

end

function z = shrink(x,r)
   z = max(0,x - r) - max(0,-x - r);
end

