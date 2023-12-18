close all;
clear;
L=10; Tmax=0.5; N=200; dt=0.0001; q=2; om=1;
dx=2*L/N;
x=(0:N)'*dx-L;
eps=10^-6;
% D0=0.1;
cc=D0(0)*dt/(dx)^2;
du=zeros(N+1,1);

% f=exp(-x.^2);
[~, idx] = min(abs(x - 0));
f=zeros(N+1,1);
f(idx)=1/dx;
urk=f;
f0=f;
% f=x;
v_new=f;
du=zeros(N+1,1);
Tri=(-2*diag(ones(N+1,1))+diag(ones(N,1),1)+diag(ones(N,1),-1));
Tri(1,:)=0;
Tri(N+1,:)=0;

v_new=f;
om=1;

for k=1:Tmax/dt
    cc=D0(k*dt)*dt/(dx)^2;
    TriD=Tri*cc;
    L=-1*diag(-q*cc*(f(1:N).^(q-1)),-1);
    D=diag(1+2*cc*q*(f.^(q-1)));
    U=-1*diag(-cc*q*(f(2:N+1).^(q-1)),1);
    A=L+D+U;
    F=v_new-TriD*(v_new.^q)-f;
    while 1
        % du=(1-om)*du+om*inv(D-L)*(U*du-F);
        du=-A\F;
        v_new=v_new+du;
        F=v_new-TriD*(v_new.^q)-f;
        if(norm(du)<eps) 
            break;
        end
    end
    f=v_new;
    % f(1)=-10;
    % f(N+1)=10;
%%%%%%%%%%%%%%RK4%%%%%%%%%%%%%%
    f1=TriD*(urk.^q)/dt;
    f2=TriD*((urk+dt*f1/2).^q)/dt;
    f3=TriD*((urk+dt*f2/2).^q)/dt;
    f4=TriD*((urk+dt*f3).^q)/dt;

    urk = urk + dt * (f1 + 2*f2 + 2*f3 + f4) / 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ((mod(k*dt,0.05)==0)&&(k*dt<0.26))
    figure(1)
    % clf    
    hold on
    xlabel('x');
    ylabel('u');
    % title('t=',k*dt);
    txt = ['t = ',num2str(k*dt)];
    plot(x,f,'b',x,f0,'r','DisplayName',txt);
    legend();
    % hold on
    % plot(x,x-k*dt,'r')
    hold off
end
end

% uc=exp(-x.^2/(4*D0*Tmax))/(4*pi*D0*Tmax)^0.5;
figure(2)
hold on
title('m=2 T_{end}=0.5');
xlabel('x');
ylabel('u');
plot(x,f,'*b',x,urk,'r')
legend('numerical','RK4');
% err=abs(f-uc);