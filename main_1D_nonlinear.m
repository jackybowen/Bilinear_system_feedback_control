close all;clear all;clc
% Stablizing controller design for nonlinear system(single control input)
% Assuming#1 monomial basis functions being used
% Assuming#2 A matrix corresponding to f(x)=x(1-x^2)
% tic
%% dynamic system formulation
f = @(t,x) x*(1-x^2);% +B*u
N = 5; % Number of monomial basis functions
A = diag([0:1:N-1]);
A0 = diag([0:1:N-3],2);A = A - A0;
B = 1;

B1 = diag(B*[1:N-1],-1);

beta = 10;
alpha = 10;
x = sym('x');
for i = 1:N
    Psi(i) = x^(i-1);
end
Dpsi = jacobian(Psi,x);
for i = 1:N
    B0 = Dpsi*B*Psi(i);
    for j = 1:N
        Bcoeff = fliplr(sym2poly(B0(j)));
        m = min(length(Bcoeff),N);
        Bi(j,:,i) = zeros(1,N);
        Bi(j,1:m,i) = Bcoeff(1:m);
    end
end


% pause
%% Test for stability conditions
% P = eye(N); % P should be symmetric, positive definite
cvx_begin sdp

variable P(N,N) symmetric
variable t
minimize(t)%-trace(B1'*P+P*B1))
subject to
t*eye(N) - (A'*P + P*A) == semidefinite(N)
P - 1e-3*eye(N) == semidefinite(N)
alpha*eye(N) - P == semidefinite(N)
cvx_end
% cvx_begin sdp
% variable P(N,N) symmetric
% maximize(trace(B1'*P+P*B1))
% subject to
% % -(A'*P + P*A) == semidefinite(N)
% P - 1e-3*eye(N) == semidefinite(N)
% alpha*eye(N) - P == semidefinite(N)
% % trace(B'*P+P*B)<=100
% cvx_end
pause
x = -10:0.1:10;
for i = 1:length(x)
    psi_x = subs(Psi,'x',x(i))';
    d = zeros(N,1);
    diBi = 0;
    diBiT= 0;
    for j = 1:size(Bi,3)
        d(j) = -alpha*(Bi(:,:,j)*psi_x)'*P*psi_x;
        diBi = diBi + d(j)*Bi(:,:,j);
        diBiT = diBiT + d(j)*Bi(:,:,j)';
    end
    dV(i) = psi_x'*(A'*P+P*A + diBiT*P + P*diBi)*psi_x;
end
figure
plot(x,dV)
xlabel('x')
ylabel('Vdot')
max(dV)
%% Closed-loop simulation
x0 = 4*rand(1,20)-2;
clearvars d
for j = 1:size(Bi,3)
    d(j) = -beta*(Bi(:,:,j)*Psi.').'*P*Psi.';
end
u = d*Psi.';
syms t x
f_c = matlabFunction(f+B*u,'Vars',[t,x]);
figure
for i = 1:length(x0)
    [t,x] = ode45(f_c,[0 0.2],x0(i));
    plot(t,x)
    hold on
end
xlabel('t')
ylabel('x')
title(['N=' num2str(N) ',\beta=' num2str(beta) ',\alpha=' num2str(alpha)])
% title(['\alpha=' num2str(alpha)]) 