function NR_Burgers_KS_RKCB2c_PS
% function <a href="matlab:NR_Burgers_KS_CNRKW3_PS">NR_Burgers_KS_CNRKW3_PS</a>
% Simulate the 1D Burgers or KS equation on 0<x<L with periodic BCs using CN/RKW3 in time
% (explicit on nonlinear terms, implicit on linear terms) & pseudospectral in space.
% See <a href="matlab:NRweb">Numerical Renaissance: simulation, optimization, & control</a>, Section 11.2.2.
% Part of <a href="matlab:help NRC">Numerical Renaissance Codebase 1.0</a>, <a href="matlab:help NRchap11">Chapter 11</a>; please read the <a href="matlab:help NRcopyleft">copyleft</a>.
% See also NR_Burgers_CNRKW3_FD and NR_Burgers_CNRKW3_FD_RS with FD implementations.

%%%%%%%%%%%%%%%%%%%% Initialize the simulation paramters (user input) %%%%%%%%%%%%%%%%%%%%
L=1000  ; Tmax=100; N=8192; dt=0.05; PlotInt=10; alpha=1; % alpha=0 for Burgers, alpha=1 for KS
dx=L/N; x=(0:N-1)'*dx; u=0.15*randn(N,1); uhat=NR_RFFT(u,N);
%%%%%%%%%%%% Precalculate the time-stepping coefficients used in the simulation %%%%%%%%%%
h_bar=dt*[8/15 2/15 1/3]; beta_bar=[1 25/8 9/4]; zeta_bar=[0 -17/8 -5/4];
kx=(2*pi/L)*[0:N/2-1]'; if alpha==0; Aop=-kx.^2; else Aop=kx.^2-kx.^4; end;
hb2=h_bar/2; hbbb=beta_bar.*h_bar; hbzb=zeta_bar.*h_bar; Imhb2=1-h_bar/2;

Aim=dt*[0 0 0 0; 0 3375509829940/4525919076317 0 0; 0 -11712383888607531889907/32694570495602105556248 566138307881/912153721139 0; 0 673488652607/2334033219546 493801219040/853653026979 184814777513/1389668723319];
Aex=dt*[0 0 0 0; 3375509829940/4525919076317 0 0 0; 0 272778623835/1039454778728 0 0; 0 673488652607/2334033219546 1660544566939/2334033219546 0];
b=dt*[0 673488652607/2334033219546 493801219040/853653026979 184814777513/1389668723319];
for k=1:Tmax/dt
  for rk=1:4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL 3 RK SUBSTEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uhat(fix(N/3)+1:end)=0;
    if (rk==1)                                                       % 2 FFTs per RK step
        y=uhat;
    else
        y=uhat+(Aim(rk,rk-1)-b(rk-1))*(Aop.*y) + (Aex(rk,rk-1)-b(rk-1))*g;
    end
      y=y./(1-Aim(rk,rk)*Aop);
      r=NR_RFFTinv(y,N); r=-0.5*r.*r; g=i*kx.*NR_RFFT(r,N);                % Leading-order cost:
      uhat=uhat+b(rk)*(Aop.*y)+b(rk)*g;

  end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF RK LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  rs(k,:)=NR_RFFTinv(uhat,N)'; ts(k)=k*dt; % These variables are just used for plotting...
  % if (mod(k,PlotInt)==0) 
  %   pause(0.001); NR_PlotXY(x,rs(k,:),k*dt,0,L/100,-1.5,1.5);
  %   % % Comment out the lines below to make some additional interesting plots.
  %   % figure(2); semilogy(kx(1:fix(N/3)),abs(uhat(1:fix(N/3))).^2); axis([0 3 1e-20 1e-1])
  %   % figure(3); loglog(kx(1:fix(N/3)),abs(uhat(1:fix(N/3))).^2); axis([3e-2 4 1e-20 1e-1])
  % end
end

figure(4); rs(:,N+1)=rs(:,1); xs=[0:N]*L/N;
contour(xs,ts,rs,[.25 .75 1.25],'r-'); hold on; contour(xs,ts,rs,[-.25 -.75 -1.25],'b-.')
end % function NR_Burgers_KS_CNRKW3_PS
