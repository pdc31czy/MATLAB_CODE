%2D AC
%save('AC_2D_data.mat','x','y','tt','uuExact')
clc;
clear all;
close all;
tic;
%discretization
a = -pi; b = pi;
t = 0;
M = 511; %the number of mesh of space of x-direction
N = 511; %the number of mesh of space of y-direction
T = 100;%the number of mesh of time %pinn200
dx =(b-a)/M; %space step
dy = (b-a)/N;
dt = 0.01; % time step size %pinn0.005
x = a:dx:b;
y = a:dy:b;
uExact = zeros(M+1,N+1);


tt=zeros(T+1,1);
tt(1)=0;

%Exact Solution for t = 0
for i = 1:M+1
    for j = 1:N+1
        uExact(i,j) = exp(-2*0.0001*t)*sin(x(i)+y(j));
    end
end
uuExact = zeros(M+1,N+1,T+1);
uuExact(:,:,1)=uExact;

for k = 1:T %time loop
    t = t+dt; %next time level
    tt(k+1) = t;
    
    %Exact Solution
    for i = 1:M+1
        for j = 1:N+1
        uExact(i,j) = exp(-2*0.0001*t)*sin(x(i)+y(j));
        end
    end
    uuExact(:,:,k+1)=uExact;

    mesh(x,y,uExact); 
    title('Exact Solution U');


    %pause(0.5);

end
toc;