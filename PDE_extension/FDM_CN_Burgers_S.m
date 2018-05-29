clear all;close all;clc
%Scheme implimentation


Basic_CN_S



for i=1:1:nt-1

    % making matrix D in AU=D
    D(1) = 0;
    D(nx) = 0;
    for j=2:1:nx-1
            D(j) = 0.5*s*U(i,j-1) + (1-s)*U(i,j) + 0.5*s*U(i,j+1);
%             D(j) = D(j) + 0.01*randn*((dt*dx)^(-1/2))*dt;
    end

    %making matrix A in AU=D
    
    A(1,1) = 1;
    A(nx,nx) = 1;
    for j=2:1:nx-1
        A(j,j-1) = -1*dt/(4*dx)*U(i,j-1) - s/2;
        A(j,j) = 1+s;
        A(j,j+1) = dt/(4*dx)*U(i,j+1) - s/2;
    end
    
    %making Inverse matrix B
    B = inv(A);
    
    %Solution for time t = i+1
    
    Ut = B * D;
    
    %setting values of Ut to U(x,t+1)
    
    for j=2:1:nx-1
        U(i+1,j) = Ut(j);
    end

end

%Making the surface
Ur1 =U;
x = linspace(x0,x1,nx);
t = linspace(t0,t1,nt);
mesh(t,x,U.');
% title('Surface plot of Burgers Equation solution(r=0.01)');
ylabel('$x$','Interpreter','Latex');
xlabel('$t$','Interpreter','Latex');

%Making the graphs
%making plot at t=0.1
%Ux = zeros(nx);
%ux = zeros(nx);
%for i=1:1:nx
%    Ux(i)=U(20,i);
%    ux(i)=U(40,i);
%end
%plot(x,Ux,x,ux);