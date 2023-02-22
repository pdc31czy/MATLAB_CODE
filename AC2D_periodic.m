%Solving 2D Allen-Cahn Eq using exact solution
%u_t= epsilon^2*(u_{xx}+u_{yy}) - 5*(u^3 - u) + g(u), x:[-pi,pi], y:[-pi,pi] and t:[0,1]
%where g(u) = 5{e^(-6*t*epsilon^2) * (sin(x+y))^3 - e^(-2*t*epsilon^2) * sin(x+y)}
%BC = Periodic; u(-pi,y,t) = u(pi,y,t) and u(x,-pi,t) = u(x,pi,t)
%IC=u(x,y,0)=sin(x+y);
%Exact solution+: u(x,y,t) = e^(-2*t*epsilon^2)*sin(x+y)
%
%save('AC2D.mat','x','y','t','uuExact')

clear all; clc;
%discretization
a = -pi; b = pi;
c = 0; d = 1;
M = 200; %the number of mesh of space
N = 200;%the number of mesh of time 
dx = (b-a)/M; % x-direction space step
dy = (b-a)/M; % y-direction space step
dt = (d -c)/N; % time step size 
epsilon = 0.0001; 
x = a:dx:b;
y = a:dy:b;
t = c:dt:d;

%uExact = zeros(M+1,M+1);
uuExact = zeros(M+1,M+1,N+1);

%Exact Solution
for k = 1:N+1 
    uExact = zeros(M+1,M+1);
    for i = 1:M+1
        for j = 1:M+1
            uExact(i,j) = exp(-2*t(k)*epsilon^2)*sin(x(i)+ y(j));
        end
    end
    uuExact(:,:,k)=uExact;
end

