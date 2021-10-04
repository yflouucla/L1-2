clc, close all,clear all

M = 100; N = 2000; F = 20;
s = 25; % sparsity of signal
%% data
A = zeros(M,N);
r = rand(M,1);
l = 1:N;
% randomly oversampled DCT matrix
for k = 1:M
    A(k,:) = sqrt(2/M)*cos(2*pi*r(k)*(l-1)/F);
end
% compute coherence(A)
fprintf(['Coherence of A is ' num2str(coherence(A)) '\n\n'])

supp = randsample_separated(N,s,2*F);
x = zeros(N,1);
xs = randn(s,1);
x(supp) = xs;
b = A*x;


% unconstrained versions for L1 and L1-L2
pm.lambda = 1e-6; pm.delta = 1e-5;
x1uncon = lasso(A,b, pm.lambda, pm.delta,1);
x12uncon = unconstrainedL1L2(A,b,pm);

fprintf(['\nRelative error of unconstrained L1 is ' num2str(norm(x1uncon-x)/norm(x)) '\n'])
fprintf(['\nRelative error of unconstrained L1-L2 is ' num2str(norm(x12uncon-x)/norm(x)) '\n\n\n'])


%% constrained versions of L1 and L1-L2
pm.lambda =2;  
pm.delta = 10;

x1con = constrainedL1(A,b,pm);
x12con = constrainedL1L2(A,b,pm);

fprintf(['\nRelative error of constrained L1 is ' num2str(norm(x1con-x)/norm(x)) '\n'])
fprintf(['\nRelative error of constrained L1-L2 is ' num2str(norm(x12con-x)/norm(x)) '\n'])
    