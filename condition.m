clc; clear; close all

%% a
M = ceil(logspace(log10(4), 3,11)); N = M-2;
Cn = zeros(size(M));

figure(1)

for i_dim = 1:length(M)
    m = M(i_dim); n = N(i_dim);
    
    A = randn(m, n);
    
    Cn(i_dim) = cond(A);
end

loglog(M, Cn)
xlim([M(1), M(end)])
xlabel('m (n = m-2)')
ylabel('condition number')


%% b
A = randn(5,3);
A_ = [A, A(:,1)];
cn = cond(A);
cn_ = cond(A_);

%% c
eps = logspace(-7,0, 17);
Cn_ = zeros(1,length(eps));

for ee = 1:length(eps)
    An = [A, A(:,1)+eps.*randn(5,1)];
    Cn_(ee) = cond(An);
end

figure()
loglog(eps, Cn_)
xlabel('$\epsilon$', 'Interpreter', 'latex')
ylabel('condition number')
