function NR_Burgers_KS_ARK3_PS
% function <a href="matlab:NR_Burgers_KS_CNRKW3_PS">NR_Burgers_KS_CNRKW3_PS</a>
% Simulate the 1D Burgers or KS equation on 0<x<L with periodic BCs using CN/RKW3 in time
% (explicit on nonlinear terms, implicit on linear terms) & pseudospectral in space.
% See <a href="matlab:NRweb">Numerical Renaissance: simulation, optimization, & control</a>, Section 11.2.2.
% Part of <a href="matlab:help NRC">Numerical Renaissance Codebase 1.0</a>, <a href="matlab:help NRchap11">Chapter 11</a>; please read the <a href="matlab:help NRcopyleft">copyleft</a>.
% See also NR_Burgers_CNRKW3_FD and NR_Burgers_CNRKW3_FD_RS with FD implementations.

%%%%%%%%%%%%%%%%%%%% Initialize the simulation paramters (user input) %%%%%%%%%%%%%%%%%%%%
L=1000; Tmax=100; N=8192; dt=0.05; PlotInt=10; alpha=1; % alpha=0 for Burgers, alpha=1 for KS
dx=L/N; x=(0:N-1)'*dx; u=0.15*randn(N,1); uhat=NR_RFFT(u,N);
%%%%%%%%%%%% Precalculate the time-stepping coefficients used in the simulation %%%%%%%%%%
kx=(2*pi/L)*[0:N/2-1]'; if alpha==0; Aop=-kx.^2; else Aop=kx.^2-kx.^4; end;

Aim=dt*[0 0 0 0; 0 3375509829940/4525919076317 0 0; 0 -11712383888607531889907/32694570495602105556248 566138307881/912153721139 0; 0 673488652607/2334033219546 493801219040/853653026979 184814777513/1389668723319];
Aex=dt*[0 0 0 0; 3375509829940/4525919076317 0 0 0; 0 272778623835/1039454778728 0 0; 0 673488652607/2334033219546 1660544566939/2334033219546 0];
b=dt*[0 673488652607/2334033219546 493801219040/853653026979 184814777513/1389668723319];

tst=tic;
for k=1:Tmax/dt
  for rk=1:4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL 3 RK SUBSTEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uhat(fix(N/3)+1:end)=0;  % Dealias (see Section 5.7).
    if (rk==1)
        y=uhat;
    else
        y = uhat + sum(Aim(rk,1:rk-1) .* f(:,1:rk-1), 2) + sum(Aex(rk,1:rk-1) .* g(:,1:rk-1), 2);
    end
    f(:,rk)=(Aop.*y)./(1-Aim(rk,rk)*Aop);
    r=NR_RFFTinv(y+Aim(rk,rk)*f(:,rk),N); r=-0.5*r.*r; g(:,rk)=i*kx.*NR_RFFT(r,N);
  end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF RK LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uhat = uhat + sum(b(1:4).*f(:,1:4), 2) + sum(b(1:4).*g(:,1:4), 2);

  % rs(k,:)=NR_RFFTinv(uhat,N)'; ts(k)=k*dt; % These variables are just used for plotting...
  % if (mod(k,PlotInt)==0) 
  %   pause(0.001); NR_PlotXY(x,rs(k,:),k*dt,0,L,-1.5,1.5);
  %   % Comment out the lines below to make some additional interesting plots.
  %   % figure(2); semilogy(kx(1:fix(N/3)),abs(uhat(1:fix(N/3))).^2); axis([0 3 1e-8 1e-1])
  %   % figure(3); loglog(kx(1:fix(N/3)),abs(uhat(1:fix(N/3))).^2); axis([3e-2 4 1e-8 1e-1])
  % end
end
tend=toc(tst);
fprintf('time: %d\n',tend);
% fprintf('  %d',NR_RFFTinv(uhat,N)');
figure(4); rs(:,N+1)=rs(:,1); xs=[0:N]*L/N;
contour(xs,ts,rs,[.25 .75 1.25],'r-'); hold on; contour(xs,ts,rs,[-.25 -.75 -1.25],'b-.')
end % function NR_Burgers_KS_CNRKW3_PS
