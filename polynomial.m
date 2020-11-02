clc; clear; close all

figure()
%% RHS plot
A = [-512,2304,-4608,5376,-4032,2016,-672,144,-18,1];
X = 1.920:1e-3:2.080;
P = A*[ones(size(X)); X; X.^2; X.^3; X.^4; X.^5; X.^6; X.^7; X.^8; X.^9];
plot(X, P)

%% LHS plot
hold on
plot(X, (X-2).^9)


%% plot appearance
xlim([X(1), X(end)])
legend('RHS','LHS')
xlabel('x')
ylabel('p(x)')
