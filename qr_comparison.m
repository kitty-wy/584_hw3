clc; clear; close all

A1 = eye(5,5);
A2 = randn(5,3);
A3 = randn(1000,1000);
A4 = randn(1000,3);

%% perform qr: 3 methods
[Q1,R1] = qr(A1); [Q1m, R1m] = mGS(A1); [Q1hh, R1hh] = qrfactor(A1);
[Q2,R2] = qr(A2); [Q2m, R2m] = mGS(A2); [Q2hh, R2hh] = qrfactor(A2);
[Q3,R3] = qr(A3);[Q3m, R3m] = mGS(A3); [Q3hh, R3hh] = qrfactor(A3);
[Q4,R4] = qr(A4); [Q4m, R4m] = mGS(A4); [Q4hh, R4hh] = qrfactor(A4);

%% check R upper triangular
D = [istriu(R1),istriu(R2),istriu(R3),istriu(R4)];
Dm = [istriu(R1m),istriu(R2m),istriu(R3m),istriu(R4m)];
Dhh = [istriu(R1hh),istriu(R2hh),istriu(R3hh),istriu(R4hh)];

%% check Q'Q = I
U = [inf_norm(Q1'*Q1-eye(size(Q1,2))), inf_norm(Q2'*Q2-eye(size(Q2,2))), ...
    inf_norm(Q3'*Q3-eye(size(Q3,2))), inf_norm(Q4'*Q4-eye(size(Q4,2)))];
Um = [inf_norm(Q1m'*Q1m-eye(size(Q1m,2))), inf_norm(Q2m'*Q2m-eye(size(Q2m,2))), ...
    inf_norm(Q3m'*Q3m-eye(size(Q3m,2))), inf_norm(Q4m'*Q4m-eye(size(Q4m,2)))];
Uhh = [inf_norm(Q1hh'*Q1hh-eye(size(Q1hh,2))), inf_norm(Q2hh'*Q2hh-eye(size(Q2hh,2))), ...
    inf_norm(Q3hh'*Q3hh-eye(size(Q3hh,2))), inf_norm(Q4hh'*Q4hh-eye(size(Q4hh,1)))];

%% error in reproducing A
E = [inf_norm(Q1*R1-A1)/inf_norm(A1), inf_norm(Q2*R2-A2)/inf_norm(A2), ...
    inf_norm(Q3*R3-A3)/inf_norm(A3), inf_norm(Q4*R4-A4)/inf_norm(A4)];
Em = [inf_norm(Q1m*R1m-A1)/inf_norm(A1), inf_norm(Q2m*R2m-A2)/inf_norm(A2), ...
    inf_norm(Q3m*R3m-A3)/inf_norm(A3), inf_norm(Q4m*R4m-A4)/inf_norm(A4)];
Ehh = [inf_norm(Q1hh*R1hh-A1)/inf_norm(A1), inf_norm(Q2hh*R2hh-A2)/inf_norm(A2), ...
    inf_norm(Q3hh*R3hh-A3)/inf_norm(A3), inf_norm(Q4hh*R4hh-A4)/inf_norm(A4)];

%% comp. time
M = ceil(logspace(1,log10(2000),11));
T = zeros(1,length(M));
Tm = zeros(1,length(M));
Thh = zeros(1,length(M));

for i_size = 1:length(M)
    disp(M(i_size))
    A = rand(M(i_size));
    
    tic
    [Q,R] = qr(A);
    T(i_size) = toc;
    
    tic
    [Q,R] = mGS(A);
    Tm(i_size) = toc;
    
    tic
    [Q,R] = qrfactor(A);
    Thh(i_size) = toc;
end

figure()
plot(M,T,M,Tm,M,Thh)
xlim([M(1), M(end)])
legend({'MATLAB qr()','Gram-Schmidt','Householder'},'location','NW')
xlabel('matrix size')
ylabel('computation time (s)')
title('Time to QR factorization solution')
