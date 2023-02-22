%Solving 2D Allen-Cahn Eq using pseudo-spectral with Implicit/Explicit
%u_t= u_{xx}+u_{yy} + u - u^3
%where u-u^3 is treated explicitly and u_{xx} and u_{yy} is treated implicitly
%BC = Periodic
%IC=v=sin(2*pi*x)+0.001*cos(16*pi*x;
clear all; clc;

%Grid
N = 64; h = 1/N; x = h*(1:N); 
dt = .01;

%x and y meshgrid
y=x';
[xx,yy]=meshgrid(x,y);

%initial conditions
v=sin(2*pi*xx)+0.001*cos(16*pi*xx);
epsilon=.01;

%(ik) and (ik)^2 vectors in x and y direction
kx=(1i*[0:N/2-1 0 -N/2+1:-1]);
ky=(1i*[0:N/2-1 0 -N/2+1:-1]');
k2x=kx.^2;
k2y=ky.^2;

ii=1:N;
%(ik)^2 matricies in x and y direction
for m= 1:N
    uxx(m,:)= k2x(ii); %Second derivative in the x-direction
end
 
for j= 1:N 
    uyy(:,j)= k2y(ii); %Second derivative in the y-direction  
end  
        
for n = 1:500
    
    v_nl=v.^3;  %calculates nonlinear term in real space
    
    for m=1:N %FFT in x-direction on linear and nonlinear term
        v_nl(m,:) = fft(v_nl(m,:));
        v_hat(m,:)=fft(v(m,:));
    end
    
    for j=1:N %FFT in y-direction on linear and nonlinear term
        v_nl(:,j) = fft(v_nl(:,j));
        v_hat(:,j)=fft(v_hat(:,j));
    end
    
    vnew=(v_hat/dt-v_nl)./ ...
       (-(uxx+uyy)*epsilon+1/dt-1); %Implicit/Explicit timestepping
  
    for m=1:N %converts to real space in x-direction
        v(m,:)=ifft(vnew(m,:));
    end
    
    for j=1:N %converts to real space in y-direction
        v(:,j)=real(ifft(v(:,j)));
    end
    
   %Plots each timestep
   surf(v); 
   title(num2str(n)); 
   axis([0 N 0 N -1 1]); 
   view(43,22); 
   drawnow;
     
end
