clear all

figure('COlor','w')
Ne=1;                 Ni=0;
re=rand(Ne,1);          ri=rand(Ni,1);
a=[0.02*ones(Ne,1);     0.02+0.08*ri];
b=[0.2*ones(Ne,1);      0.25-0.05*ri];
c=[-65+15*0.8.^2;        -65*ones(Ni,1)];
d=[8-6*1.^2;           2*ones(Ni,1)];

S=[0.5*rand(Ne+Ni,Ne),  -rand(Ne+Ni,Ni)];

v=-65*ones(Ne+Ni,1);    % Initial values of v
u=b.*v;                 % Initial values of u
firings=[];             % spike timings
time_sim=40000;  % simulation time in ms
lowmod= medfilt1(randn(1,time_sim),100).*3;

% Used as a dummy way to affect the spike rate
dummy_rate = 4; % Default is 4

for t=1:time_sim         % simulation of 1000 ms
  I=[dummy_rate*randn(Ne,1)+lowmod(t);2*randn(Ni,1)]; % thalamic input
  fired=find(v>=30);    % indices of spikes
  firings=[firings; t+0*fired,fired];
  v(fired)=c(fired);
  u(fired)=u(fired)+d(fired);
  I=I+sum(S(:,fired),2);
  v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
  v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
  u=u+a.*(b.*v-u);  
  vx=v(1);vx(vx>30)=30;
  
  volt(t,1)=vx;
end;

noiselevel= 50;% 8; default% noise level
volt=volt+ randn(length(volt),1).*noiselevel;

plot(volt,'k')
hold on, plot(firings(:,1), firings(:,2)+40,'or')

spike_rate = (size(firings, 1))/length(volt)
size(firings)