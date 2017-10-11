close all;clear all;clc
% Stablizing controller design for 2D nonlinear system(single control input)
% Assuming#1 monomial basis functions being used
% Assuming#2 A matrix corresponding to f(x)=x(1-x^2)
% tic
n = 2; % dimension of system
x = sym('x',[n,1]);
%% Dynamic system formulation
f = [x(2); -x(1)+x(2)*(1-x(1)^2)];% +B*u
D = 5; % degree of monomial basis at most D
N = nchoosek(n+D,D);	% Number of monomial basis functions
% A = diag([0:1:N-1]);
% A0 = diag([0:1:N-3],2); A = A - A0;
g = [0; 1];
alpha = 1;
beta = 10;

%% Generate basis function for EDMD/NSDMD
Monom = repmat(x,1,D+1);
Monom(:, 1) = 1;
Monom = cumprod(Monom,2);
[X,Y] = ndgrid(Monom(1,:),Monom(2,:));
Monom = rot90(X.*Y,3);
Psi = [];
for i = D:-1:0
    Psi = [Psi diag(Monom,i).'];
end

% Psi(1) = [];
% N = N-1;

Dpsi = jacobian(Psi,x);
Dpsif = Dpsi*f;%sum(Dpsi.*repmat(f.',N,1),2)
Dpsig = Dpsi*g;%sum(Dpsi.*repmat(g.',N,1),2)

for i = 1:N
    A(i,:) = coeffs_2D(Dpsif(i),Psi);
    B(i,:) = coeffs_2D(Dpsig(i),Psi);
%     m = min(length(Acoeff),N);
%     A0 = padarray(Acoeff,[N-size(Acoeff,1),N-size(Acoeff,2)],'pre');
%     Ai = [];
%     for i = D:-1:0
%         Ai = [Ai diag(A0,i).'];
%     end
%     B0 = Dpsi*g*Psi(i);
%     for j = 1:N
%         Bi(j,:,i) = coeffs_2D(B0(j),Psi);
%     end
end
% pause
%% Controller design problem
% P = eye(N); % P should be symmetric, positive definite
cvx_begin sdp
variable P(N,N) symmetric
variable t
% P = sdpvar(N,N)
% Lm = lambda_min(B1'*P+P*B1);
minimize(t - trace(B'*P+P*B))
subject to
t*eye(N) - (A'*P + P*A) == semidefinite(N);
P - 1e-3*eye(N) == semidefinite(N);
alpha*eye(N) - P == semidefinite(N);
% trace(B'*P+P*B)<=100
cvx_end
% pause
% %% Test for stability conditions
% 
% tic
% x1 = -1:0.025:1;
% x2 = -2:0.05:2;
% [xx,yy] = meshgrid(x1,x2);
% % [xx1,yy1] = meshgrid(x1,x2);
% 
% for i = 1:length(x1)
%     for j = 1:length(x2)
%         psi_x = subs(Psi,x,[x1(i);x2(j)])';
%         d = zeros(N,1);
%         diBi = 0;
%         diBiT= 0;
%         for k = 1:size(Bi,3)
%             d(k) = -beta*(Bi(:,:,k)*psi_x)'*P*psi_x;
%             diBi = diBi + d(k)*Bi(:,:,k);
%             diBiT = diBiT + d(k)*Bi(:,:,k)';
%         end
%            dV(i,j) = psi_x'*(A'*P+P*A + diBiT*P + P*diBi)*psi_x;
%         if dV(i,j) >= 0 %&& dV(i,j)<=1e-3
%            dV(i,j) = 1e9*dV(i,j);
%         end
%     end
% end
% figure
% surf(xx,yy,eval(dV).')
% xlabel('x')
% ylabel('y')
% zlabel('$\dot{V}$','Interpreter','Latex')
% toc
% 
% max(eval(dV(:)))
% % toc
%% Closed-loop simulation
close all
digits(3)
noi = 10;
% t = rand(1,noi)*2*pi;
% r = rand(1,noi)*0.1;
% x0 = r.*cos(t)-0.4;
% y0 = r.*sin(t)-0.3;

x0 = 4*rand(1,noi)-2;
y0 = 4*rand(1,noi)-2;
% x0 = 1;
% y0 = -1.5;
% for j = 1:size(Bi,3)
%     d_t(j) = -beta*(Bi(:,:,j)*Psi.').'*P*Psi.';
% end
u = simplify(vpa(-beta*(Psi*B'*P*Psi.'*(Psi*Psi.'))));
syms t;
f_c = matlabFunction(f+g*u,'Vars',{t,x});
% figure
for i = 1:noi
    [t,z] = ode15s(f_c,[0 20],[x0(i);y0(i)]);
%     plot(z(:,1),z(:,2))
%     hold on
    figure(3)
    plot(t,z(:,1))
    hold on
    figure(4)
    plot(t,z(:,2))
    hold on
%     pause
end




% syms t;
% z = sym('z',[N,1]);
% u = -beta*z'*B'*P*z*z'*z;
% f_z = A*z+B*z*u;
% f_c2 = matlabFunction(f_z,'Vars',{t,z});
% z0 = double(vpa(subs(V'*Psi.',{'x1','x2'},{x0,y0})));
% % z0 = 4*randn(N,noi);
% figure
% for i = 1:noi
% %     [t,z] = ode45(f_c,[0 10],[x0(i);y0(i)]);
%     [t,z_t] = ode15s(f_c2,[0,10],z0(:,i));
% %     plot(z(:,1),z(:,2))
% %     hold on
% %     figure(3)
% %     plot(t,z_t(:,1))
% %     hold on
% %     figure(4)
% %     plot(t,z_t(:,2))
% %     hold on
% 
% figure
% plot(t,z_t.')
% %     pause
% end
% 
% 
% % for i = 1:noi
%     [t,xy] = ode15s(f_c1,[0 10],[x0;y0]);
% 
% figure(3)
% plot(t,xy(:,1))
% hold on
% figure(4)
% plot(t,xy(:,2))
% hold on
% figure(5)
% plot(t,xy.')
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
% figure(3)
% xlabel('t')
% ylabel('x')
% title(['D=' num2str(D) ',\beta=' num2str(beta) ',\alpha=' num2str(alpha)])
% % title(['\alpha=' num2str(alpha)])
% figure(4)
% xlabel('t')
% ylabel('y')
% title(['D=' num2str(D) ',\beta=' num2str(beta) ',\alpha=' num2str(alpha)])
