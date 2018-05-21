close all;clear all;clc
% Stablizing controller design for 2D nonlinear system(single control input)
% Assuming#1 monomial basis functions being used
% Assuming#2 A matrix corresponding to f(x)
% tic
n = 8; % dimension of system
x = sym('x',[n,1]);
%% Dynamic system formulation
% x_dot = f(x) + g(x)u

% Invented pendulum on the cart with triple links
A21 = [0 -15.1861 2.7420 -0.3023;...
       0 147.1379 -82.4455 9.0896;...
       0 -182.6242 194.3229 -41.2087;...
       0 37.9468 -77.6656 45.8506];
A22 = [-5.2356 0.0136 -0.0107 0.0027;...
       17.1258 -0.1778 0.2052 -0.0822;...
       -6.8497 0.2967 -0.4352 0.2246;...
       1.4233 -0.0924 0.2138 -0.1355];
f = [zeros(4) eye(4);A21 A22]*x;
g = [0 0 0 0 3.7397 -12.2326 4.8926 -1.0166]';


% % Van der Pol Oscillator
% f =  [x(2); -x(1)+x(2)*(1-x(1)^2)];

% % 2D nonlinear
% delta = -0.01;
% syms t;
% f = [x(2);delta*x(2)-sin(x(1))];

% % Linear System
% f = [-1 2;0 0.1]*x; % Marginally stable case
% f_u =  [1 -2;0 0.9]*x; % Unstable case

% D = 3; % degree of monomial basis at most D
% N = nchoosek(n+D,D); % Number of monomial basis functions
% g = [0; 1];

alpha = 1;
beta = 10;

%% Generate basis function for EDMD/NSDMD
% Monom = repmat(x,1,D+1);
% Monom(:, 1) = 1;
% Monom = cumprod(Monom,2);
% [X,Y] = ndgrid(Monom(1,:),Monom(2,:));
% Monom = rot90(X.*Y,3);
% Psi = [];
% for i = D:-1:0
%     Psi = [Psi diag(Monom,i).'];
% end

global Deg;
Deg = 2;
Psi = Monomials(x);

% Psi = [x(1) x(2)];
% Psi(1) = [];
N = length(Psi);
%% Approximate the (A,B) bilinear system
Tf = 5;
dt = .1;

xy_limit = [-5*ones(n,1) 5*ones(n,1)];
syms t;
Kdmd = Kpm_comp_EDMD(matlabFunction(f,'Vars',{t,x}),xy_limit,dt,Tf,matlabFunction(Psi,'Vars',{x}));
[V, E] = eig(Kdmd);

lambda = log(diag(E))/dt;
A = diag(lambda);
i = 1;
while i <=length(lambda)
    lam = lambda(i);
    if ~isreal(lam)
        A([i i+1],[i i+1]) = [real(lam) imag(lam);...
                              -imag(lam) real(lam)];
        Re = real(V(:,i));Im = imag(V(:,i));
        V(:,i) = 2*Re;
        V(:,i+1) = -2*Im;
        i = i + 1;
    end
        i = i + 1;
end

% A = -A;
Dpsi = jacobian(Psi,x);
% Dpsif = sum(Dpsi.*repmat(f.',N,1),2)
% Dpsig = sum(Dpsi.*repmat(g.',N,1),2)
Dphi = V.'*Dpsi;

% Dphif = Dphi*f;
Dphig = Dphi*g;

for i = 1:N
%     A(i,:) = coeffs_nD(Dpsif(i),Psi);
    BV(i,:) = coeffs_nD(Dphig(i),Psi);
%     m = min(length(Acoeff),N);
%     A0 = padarray(Acoeff,[N-size(Acoeff,1),N-size(Acoeff,2)],'pre');
%     Ai = [];
%     for i = D:-1:0
%         Ai = [Ai diag(A0,i).'];
%     end

%     B0 = Dphig*Psi*V(:,i);
%     for j = 1:N
%         Bk(j,:,i) = coeffs_nD(B0(j),Psi);
%     end
%     Bk(:,:,i) = Bk(:,:,i)*inv(V.');
end
% B = BV/(V.');
B = BV*pinv(V.');

% pause
% %% Controller design problem
% % P = eye(N); % P should be symmetric, positive definite
% cvx_begin sdp
% variable P(N,N) symmetric
% variable t
% % P = sdpvar(N,N)
% % Lm = lambda_min(B'*P+P*B);
% minimize(t - trace(B'*P+P*B))
% subject to
% t*eye(N) - (A'*P+P*A) == semidefinite(N);
% P - 1e-3*eye(N) == semidefinite(N);
% alpha*eye(N) - P == semidefinite(N);
% % trace(B'*P+P*B)<=100
% cvx_end
% % pause
% % %% Test for stability conditions
% % tic
% % x1 = -1:0.025:1;
% % x2 = -2:0.05:2;
% % [xx,yy] = meshgrid(x1,x2);
% % % [xx1,yy1] = meshgrid(x1,x2);
% % 
% % for i = 1:length(x1)
% %     for j = 1:length(x2)
% %         psi_x = subs(Psi,x,[x1(i);x2(j)])';
% %         d = zeros(N,1);
% %         diBi = 0;
% %         diBiT= 0;
% %         for k = 1:size(Bk,3)
% %             d(k) = -beta*(Bk(:,:,k)*psi_x)'*P*psi_x;
% %             diBi = diBi + d(k)*Bk(:,:,k);
% %             diBiT = diBiT + d(k)*Bk(:,:,k)';
% %         end
% %            dV(i,j) = psi_x'*(A'*P+P*A + diBiT*P + P*diBi)*psi_x;
% %         if dV(i,j) >= 0 %&& dV(i,j)<=1e-3
% %            dV(i,j) = 1e9*dV(i,j);
% %         end
% %     end
% % end
% % figure
% % surf(xx,yy,eval(dV).')
% % xlabel('x')
% % ylabel('y')
% % zlabel('$\dot{V}$','Interpreter','Latex')
% % toc
% % 
% % max(eval(dV(:)))
% % toc
% %% Closed-loop simulation
% close all
% % digits(3)
% noi = 1;
% idx = find(abs(diag(A))<=1e-4)
% 
% % t = rand(1,noi)*2*pi;
% % r = rand(1,noi)*0.1;
% % x0 = r.*cos(t)-0.4;
% % y0 = r.*sin(t)-0.3;
% beta = 1e5;
% x0 = 2*rand(1,noi)-1;
% y0 = 2*rand(1,noi)-1;
% % x0 = 1;
% % y0 = -1.5;
% Phi = V'*Psi.';
% u = simplify(-beta*(Phi.'*B'*P*Phi*(Phi.'*Phi)));
% u = u - vpa(subs(u,{'x1','x2'},{0,0}));
% syms t;
% f_c1 = matlabFunction(f+g*u,'Vars',{t,x});
% for i = 1:noi
%     [t,xy] = ode15s(f_c1,[0 100],[x0;y0]);
% 
% %     figure
% %     plot(t,xy(:,1))
% %     hold on
% %     xlabel('t')
% %     ylabel('x')
% % 
% %     figure
% %     plot(t,xy(:,2))
% %     hold on
% %     xlabel('t')
% %     ylabel('y')
% % 
% %     figure
% %     plot(xy(:,1),xy(:,2))
% %     xlabel('x')
% %     ylabel('y')
% 	figure
%     plot(t,xy.')
%     xlabel('t')
%     ylabel('x and y')
% %     pause
% end
% 
% % z0 = double(vpa(subs(V'*Psi.',{'x1','x2'},{x0,y0})));
% % z0(idx) = 0;
% % % z0 = 4*randn(N,noi);
% % z = sym('z',[N,1],'real');
% % u = -beta*z'*B'*P*z*z'*z;
% % f_z = A*z+B*z*u;
% % syms t;
% % f_c2 = matlabFunction(f_z,'Vars',{t,z});
% % for i = 1:noi
% % %     [t,z] = ode45(f_c,[0 10],[x0(i);y0(i)]);
% %     [t,z_t] = ode15s(f_c2,[0,10],z0(:,i));
% % %     plot(z(:,1),z(:,2))
% % %     hold on
% % %     figure(3)
% % %     plot(t,z_t(:,1))
% % %     hold on
% % %     figure(4)
% % %     plot(t,z_t(:,2))
% % %     hold on
% % 
% % figure
% % plot(t,z_t.')
% % xlabel('t')
% % ylabel('z_t')
% % %     pause
% % end
% 
% 
% 
% % syms t;
% % xz = sym('xz',[N+2 1]);
% % z = xz(3:end);
% % u = -beta*z'*B'*P*z*z'*z;
% % f_xz = [xz(2);-xz(1)+xz(2)*(1-xz(1)^2)+u; A*z + B*z*u];
% % syms t;
% % f_c3 = matlabFunction(f_xz,'Vars',{t,xz});
% % for i = 1:noi
% %     [t,xz_t] = ode15s(f_c3,[0,1000],[x0;y0;z0(:,i)]);
% % figure
% % plot(t,xz_t(:,[1 2]).')
% % xlabel('t')
% % ylabel('x and y')
% % figure
% % plot(t,xz_t(:,3:end).')
% % xlabel('t')
% % ylabel('z')
% % %     pause
% % end
% % 


