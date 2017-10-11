close all;clear all;clc
% Stablizing controller design for 3D nonlinear system(single control input)
% Assuming#1 monomial basis functions being used
% Assuming#2 A matrix corresponding to Lorentz attractor
% tic
n = 3; % dimension of system
x = sym('x',[n,1]);
%% dynamic system formulation
sigma = 10;
beta = 8/3;
rho = 28;
f = [sigma*(x(2)-x(1));...
    x(1)*(rho -x(3))-x(2);...
    x(1)*x(2)-beta*x(3)];% +B*u
D = 3; % degree of monomial basis at most D
N = nchoosek(n+D,D); % Number of monomial basis functions
% A = diag([0:1:N-1]);
% A0 = diag([0:1:N-3],2);A = A - A0;
B = [0; 1; 1];

beta = 1;
alpha = 10;

Monom = repmat(x,1,D+1);
Monom(:, 1) = 1;
Monom = cumprod(Monom,2);
[X,Y,Z] = ndgrid(Monom(1,:),Monom(2,:),Monom(3,:));
Monom = X.*Y.*Z;

[I,J,K] = ndgrid(1:D+1,1:D+1,1:D+1);
Sum = I+J+K;
idx = find(Sum<=3+D);
[~,idx0] = sort(Sum(idx),'ascend');
Psi = Monom(idx(idx0)).'; % Monomial basis functions
Dpsi = jacobian(Psi,x);
Dpsidt = sum(Dpsi.*repmat(f.',N,1),2)
Dpsidt1 = sum(Dpsi.*repmat(B.',N,1),2)
for i = 1:N
    A(i,:) = coeffs_3D(Dpsidt(i),Psi);
    B1(i,:) = coeffs_3D(Dpsidt1(i),Psi);
%     m = min(length(Acoeff),N);
%     A0 = padarray(Acoeff,[N-size(Acoeff,1),N-size(Acoeff,2)],'pre');
%     Ai = [];
%     for i = D:-1:0
%         Ai = [Ai diag(A0,i).'];
%     end
    B0 = Dpsi*B*Psi(i);
    for j = 1:N
        Bi(j,:,i) = coeffs_3D(B0(j),Psi);
    end
end
pause

%% Test for stability conditions
% P = eye(N); % P should be symmetric, positive definite
cvx_begin sdp
variable P(N,N) symmetric
maximize(trace(B1'*P+P*B1))
subject to
A'*P + P*A == semidefinite(N)
P - 1e-5*eye(N) == semidefinite(N)
% alpha*eye(N) - P == semidefinite(N)
% trace(B'*P+P*B)<=100
cvx_end

% pause
% x1 = -20:2:20;
% x2 = -30:2:30;
% x3 = 0:2:50;
% for i = 1:length(x1)
%     for j = 1:length(x2)
%         for k = 1:length(x3)
%             psi_x = subs(Psi,x,[x1(i);x2(j);x3(k)])';
%             d = zeros(N,1);
%             diBi = 0;
%             diBiT= 0;
%             for m = 1:size(Bi,3)
% %                 tic
%                 d(m) = -beta*(Bi(:,:,m)*psi_x)'*P*psi_x;
%                 diBi = diBi + d(m)*Bi(:,:,m);
%                 diBiT = diBiT + d(m)*Bi(:,:,m)';
% %                 toc
%             end
%             dV(i,j,k) = psi_x'*(A'*P+P*A + diBiT*P + P*diBi)*psi_x;
%         end
%     end
% end
% 
% [xx,yy,zz] = meshgrid(x1,x2,x3);
% surf(xx,yy,eval(dV).')
% xlabel('x')
% ylabel('y')
% zlabel('$\dot{V}$','Interpreter','Latex')
% 
% max(eval(dV(:)))

%% Closed-loop simulation
noi = 1;
% z0 = diag([40 60 50])*rand(3,noi)-repmat([20,30,0]',1,noi)
z0 = [2;-2;10];
for j = 1:size(Bi,3)
    d_t(j) = -beta*(Bi(:,:,j)*Psi.').'*P*Psi.';
end
u = d_t*Psi.';
syms t;
f_c = matlabFunction(f+B*u,'Vars',{t,x});
figure
for i = 1:noi
    [~,z] = ode45(f_c,0:0.1:1,z0(:,i));
    plot3(z(:,1),z(:,2),z(:,3))
    hold on
end
xlabel('x')
ylabel('y')
zlabel('z')
title(['D=' num2str(D) ',\alpha=' num2str(beta)])
% title(['\alpha=' num2str(alpha)])