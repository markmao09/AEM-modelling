clear

%% Inputs
F=96485; %C/mole, faraday constant
R=8.314; %J/(mole*K), gas constant
T=60+273.15; %K, temperature
f=F/R/T; %C/J=1/V

L=1e-4; %m, width of each region (100 µm)
n=100; %--, number points in each region [if increase n, might have to decrease dt]
dx=L/n; %m, space discretisation
x=linspace(0,3*L,3*n); %m

U=-0.6; %V, applied potential difference
E=(U-0)/(3*L-0); %V/m, linear potential gradient in regions 1 and 3

C0=5; %M, K+ and OH- liquid electrolyte bulk concentration
zk=1; %--, K+ valence
% Dk=1.96e-9; %m2/s, K+ diffusion coefficient
% Dkm=1.96e-10; %m2/s, K+ diffusion coefficient inside membrane
zh=-1; %--, OH- valence
Dh=5.27e-9; %m2/s, OH- diffusion coefficient
Dhm=5.27e-10; %m2/s, OH- diffusion coefficient inside membrane
Cm=5; %M, concentration of stationary fixed charges in membrane

% i0her=1e-3; %A/m2, exchange current density HER
% i0oer=1e-3; %A/m2, exchange current density OER
% A=0.8e1; %1/m, electrode surface area to volume ratio
% alfa=0.5; %--, transfercoëfficiënt for both HER and OER
% Eeq=1.229-0.9e-3*(T-298); %V, equilibrium potential

%% Initial condition (j=1)
CM=(1/2)*(sqrt(Cm^2+4*C0^2)+Cm); %concentration in membrane from equilibrium and electroneutrality
CMl=CM; %left interface
CMr=CM; %right interface
EM=E; %begin value of membrane potential gradient (where no potential U is applied yet)
C1(:,1)=C0*ones(n,1); %concentration region 1
C2(:,1)=CM*ones(n,1); %concentration region 2
C3(:,1)=C0*ones(n,1); %concentration region 3

%% Make space discretisation matrices for regios 1 and 3
dt=1e-5; %s, time discretisation [MUST BE SMALL!]
a=Dh*dt/dx^2;
b=zh*f/2*dx*E;
% Region 1 discretisation matrix, dimensions n-1 x n:
R1=[diag(-2*a*ones(n-1,1))+diag(a*(1+b)*ones(n-2,1),1)+diag(a*(1-b)*ones(n-2,1),-1),...
    [zeros(n-2,1);a*(1+b)]];
% Region 3 discretisation matrix, dimensions n-1 x n:
R3=[[a*(1-b);zeros(n-2,1)],diag(-2*a*ones(n-1,1))+diag(a*(1+b)*ones(n-2,1),1)+...
    diag(a*(1-b)*ones(n-2,1),-1)];

%% Explicit finite difference method: forward euler time integration
d=1e4; %duration
for j=2:d %end time (total elapsed time) is d*dt seconds
Em=EM(j-1); %V/m, linear potential gradient in region 2
am=Dhm*dt/dx^2;
bm=zh*f/2*dx*Em;
% Region 2 discretisation matrix, dimensions n x n:
R2=diag(-2*am*ones(n,1))+diag(am*(1+bm)*ones(n-1,1),1)+diag(am*(1-bm)*ones(n-1,1),-1);

% Calculate concentration profile at next timestep:
C1(1:n-1,j)=C1(1:n-1,j-1)+R1*C1(:,j-1)+[a*(1-b)*C0;zeros(n-2,1)]; %region 1
% ioer=i0oer*(exp(-alfa*f*(E*x(n-1)-Eeq))-(C1(n-1,j-1)/C0)*exp((1-alfa)*f*(E*x(n-1)-Eeq))); %A/m2
% C1(n-1,j)=C1(n-1,j-1)+R1(end,:)*C1(:,j-1)-A*ioer/F*dt; %mole/m3, OER

C2(1:n,j)=C2(1:n,j-1)+R2*C2(:,j-1); %region 2
C2(1,j)=C2(1,j)+am*(1-bm)*CMl(j-1); %left interface BC
C2(n,j)=C2(n,j)+am*(1+bm)*CMr(j-1); %right interface BC

C3(2:n,j)=C3(2:n,j-1)+R3*C3(:,j-1)+[zeros(n-2,1);a*(1+b)*C0]; %region 3
% iher=i0her*(exp(-alfa*f*(E*x(2*n+2)-Eeq))-(C3(2,j-1)/C0)*exp((1-alfa)*f*(E*x(2*n+2)-Eeq))); %A/m2
% C3(2,j)=C3(2,j-1)+R3(1,:)*C3(:,j-1)+A*iher/F*dt; %mole/m3, HER

% Calculate interface concentrations CMl and C1n:
k0=-Dhm/dx*C2(1,j) + Dhm*Cm/2/dx - Dhm*0.5*Cm*zh*f*Em - Dh/dx*C1(n-1,j); %constant 0
k1=-Dh/dx - Dh*zh*f*E; %constant 1
k2=Dhm*zh*f*Em*0.5 - Dhm/dx/2; %constant 2
C1(n,j)=min(roots([(4*k2^2-k1^2),2*k1*k0,k2^2*Cm^2-k0^2])); %has to be min!
CMl(j)=(1/2)*(sqrt(Cm^2+4*C1(n,j)^2)+Cm);

% Calculate interface concentrations CMr and C31:
k0=Dhm*Cm/2/dx - Dhm*C2(n,j)/dx + Dhm*zh*f*Em*0.5*Cm - Dh*C3(2,j)/dx; %constant 0
k1=Dh*zh*f*E - Dh/dx; %constant 1
k2=-Dhm/2/dx - Dhm*zh*f*Em*0.5; %constant 2
C3(1,j)=min(roots([(4*k2^2-k1^2),2*k1*k0,k2^2*Cm^2-k0^2])); %has to be min!
CMr(j)=(1/2)*(sqrt(Cm^2+4*C3(1,j)^2)+Cm);

% Calculate new donnan potentials at each interface and new potential gradient in membrane:
DPl=R*T/zh/F*log(CMl(j)/C1(n,j)); %left interface
DPr=R*T/zh/F*log(C3(1,j)/CMr(j)); %right interface
EM(j)=E+(DPr-DPl)/L;
end

%% Make concentrations and potential vectors
Coh=cat(1,C1,C2,C3); %M
Ck1=-zh*C1/zk;
Ck2=(-zh*C2-Cm)/zk;
Ck3=-zh*C3/zk;
Ck=cat(1,Ck1,Ck2,Ck3); %M
E1=E*x(1:n);
E3=E*x(2*n+1:3*n);
index=[1,d];
tijd=[0,d*dt]; %s, times correspond to index vector above
for t=1:length(index)
    E2(t,:)=EM(index(t))*x(n+1:2*n); %#ok<SAGROW>
    Pot(t,:)=[E1,E2(t,:),E3]; %#ok<SAGROW> %V
end

%% Plot
close all
% Concentrations plot:
figure(1)
hold on
grid on
xlabel('Position (m)')
ylabel('Concentration (M)')
xticks(linspace(0,L*3,7))
xlim([0,L*3])
for t=1:length(index)
    plot(x,Coh(:,index(t)),'-','DisplayName',sprintf('OH^– %g s',tijd(t)));
    plot(x,Ck(:,index(t)),':','DisplayName',sprintf('K^+ %g s',tijd(t)));
end
legend('AutoUpdate','off')
rect=area([L,L*2],[max(get(gca,'YLim')),max(get(gca,'YLim'))],min(get(gca,'YLim')),...
    'EdgeColor','none','ShowBaseLine','off','FaceColor','#999999','FaceAlpha',0.1);
uistack(rect,'bottom')
% Potential plot:
figure(2)
hold on
grid on
xlabel('Position (m)')
ylabel('Potential (V)')
xticks(linspace(0,L*3,7))
xlim([0,L*3])
for t=1:length(index)
    plot(x,Pot(t,:),'-','DisplayName',sprintf('%g s',tijd(t)));
end
legend('AutoUpdate','off')
rect=area([L,L*2],[max(get(gca,'YLim')),max(get(gca,'YLim'))],min(get(gca,'YLim')),...
    'EdgeColor','none','ShowBaseLine','off','FaceColor','#999999','FaceAlpha',0.1);
uistack(rect,'bottom')