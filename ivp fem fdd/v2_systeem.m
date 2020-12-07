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
E=U/(3*L); %V/m, linear potential gradient in regions 1 and 3

C0=5; %M, K+ and OH- liquid electrolyte bulk concentration
zk=1; %--, K+ valence
Dk=1.96e-9; %m2/s, K+ diffusion coefficient
zoh=-1; %--, OH- valence
Doh=5.27e-9; %m2/s, OH- diffusion coefficient
D=(Doh*Dk*(zoh-zk))/(Doh*zoh-Dk*zk); %m2/s, KOH electrolyte averaged diffusion coefficient
Dm=D*0.5; %m2/s, averaged diffusion coefficient in membrane
Cm=5; %M, concentration of stationary fixed charges in membrane

j0her=1e-3; %A/m2, exchange current density HER
j0oer=1e-3; %A/m2, exchange current density OER
A=1e0; %1/m, electrode surface area to volume ratio
alfa_a=1; %--, transfer coefficient for OER
alfa_c=0.5; %--, transfer coefficient for HER
Eeq=1.229-0.9e-3*(T-298); %V, equilibrium potential

%% Initial condition (j=1)
CM=(1/2)*(sqrt(Cm^2+4*C0^2)+Cm); %concentration in membrane from equilibrium and electroneutrality
C1M=CM; %left interface
C3M=CM; %right interface
Em=E; %V/m, begin value of membrane potential gradient (where no potential U is applied yet)
C1(:,1)=C0*ones(n,1); %concentration region 1
C2(:,1)=CM*ones(n,1); %concentration region 2
C3(:,1)=C0*ones(n,1); %concentration region 3

%% Make space discretisation matrices for regions 1 and 3
dt=1e-5; %s, time discretisation [MUST BE SMALL!]
a=D*dt/dx^2;
b=zoh*f/2*dx*E;
S=diag(-2*a*ones(n-1,1)) + diag(a*(1+b)*ones(n-2,1),1) + diag(a*(1-b)*ones(n-2,1),-1);
% Region 1 discretisation matrix, dimensions n-1 x n:
R1=[S,[zeros(n-2,1);a*(1+b)]];
% Region 3 discretisation matrix, dimensions n-1 x n:
R3=[[a*(1-b);zeros(n-2,1)],S];

d=1e5; %duration
for j=2:d %end time (total elapsed time) is d*dt seconds
%% Make space discretisation matrix for region 2
am=Dm*dt/dx^2;
bm=zoh*f/2*dx*Em(j-1);
% Region 2 discretisation matrix, dimensions n x n:
R2=diag(-2*am*ones(n,1)) + diag(am*(1+bm)*ones(n-1,1),1) + diag(am*(1-bm)*ones(n-1,1),-1);

%% Calculate concentration profile at next timestep:
C1(1:n-1,j)=C1(1:n-1,j-1)+R1*C1(:,j-1)+[a*(1-b)*C0;zeros(n-2,1)]; %region1
joer=j0oer*(exp(-alfa_a*f*(E*x(n-1)-Eeq))-(C1(n-1,j-1)/C0)*exp((1-alfa_a)*f*(E*x(n-1)-Eeq))); %A/m2
C1(n-1,j)=C1(n-1,j-1)+R1(end,:)*C1(:,j-1)-A*joer/F*dt; %mole/m3, OER, comment out if no reaction
for e=1:length(C1(:,j))
    if C1(e,j)<0 %prevent complex numbers on anode side
        C1(e,j)=0;
    end
end

C2(1:n,j)=C2(1:n,j-1)+R2*C2(:,j-1)+[am*(1-bm)*C1M(j-1);zeros(n-2,1);am*(1+bm)*C3M(j-1)]; %region2

C3(2:n,j)=C3(2:n,j-1)+R3*C3(:,j-1)+[zeros(n-2,1);a*(1+b)*C0]; %region3
jher=j0her*(exp(-alfa_c*f*(E*x(2*n+2)-Eeq))-(C3(2,j-1)/C0)*exp((1-alfa_c)*f*(E*x(2*n+2)-Eeq))); %A/m2
C3(2,j)=C3(2,j-1)+R3(1,:)*C3(:,j-1)+A*jher/F*dt; %mole/m3, HER, comment out if no reaction

%% Calculate interface concentrations: C1(n), C1M(n), C3M(1) and C3(1)
% Calculate interface concentrations C1(n) and C1M(n):
k0=-Dm/dx*C2(1,j) + Dm*Cm/2/dx - Dm*0.5*Cm*zoh*f*Em(j-1) - D/dx*C1(n-1,j); %constant 0
k1=-D/dx - D*zoh*f*E; %constant 1
k2=Dm*zoh*f*Em(j-1)*0.5 - Dm/dx/2; %constant 2
aa=4*k2^2-k1^2;
bb=2*k1*k0;
cc=k2^2*Cm^2-k0^2;
C1(n,j)=min([(-bb+sqrt(bb^2-4*aa*cc))/(2*aa),(-bb-sqrt(bb^2-4*aa*cc))/(2*aa)]);
C1M(j)=(1/2)*(sqrt(Cm^2+4*C1(n,j)^2)+Cm);

% Calculate interface concentrations C3M(1) and C3(1):
k0=Dm*Cm/2/dx - Dm*C2(n,j)/dx + Dm*zoh*f*Em(j-1)*0.5*Cm - D*C3(2,j)/dx; %constant 0
k1=D*zoh*f*E - D/dx; %constant 1
k2=-Dm/2/dx - Dm*zoh*f*Em(j-1)*0.5; %constant 2
aa=4*k2^2-k1^2;
bb=2*k1*k0;
cc=k2^2*Cm^2-k0^2;
C3(1,j)=min([(-bb+sqrt(bb^2-4*aa*cc))/(2*aa),(-bb-sqrt(bb^2-4*aa*cc))/(2*aa)]);
C3M(j)=(1/2)*(sqrt(Cm^2+4*C3(1,j)^2)+Cm);

%% Calculate new donnan potentials at each interface and new potential gradient in membrane:
DPleft=1/f/zoh*log(C1M(j)/C1(n,j)); %left interface
DPright=1/f/zoh*log(C3(1,j)/C3M(j)); %right interface
Em(j)=E+(DPright-DPleft)/L;

if C1(n,j)<0 %prevent complex numbers on anode side
    C1(n,j)=0;
    break
end
end

%% Make concentrations and potential vectors
Coh=cat(1,C1,C2,C3); %M
Ck1=-zoh*C1/zk;
Ck2=(-zoh*C2-Cm)/zk;
Ck3=-zoh*C3/zk;
Ck=cat(1,Ck1,Ck2,Ck3); %M
E1=E*x(1:n);
E3=E*x(2*n+1:3*n);
% [~,mi]=max(Em(2:end));
index=[j];
tijd=[j*dt]; %s, times correspond to index vector above
for t=1:length(index)
    E2(t,:)=Em(index(t))*x(n+1:2*n); %#ok<*SAGROW>
    Pot(t,:)=[E1,E2(t,:),E3]; %V
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
%     plot(x,Ck(:,index(t)),'--','DisplayName',sprintf('K^+ %g s',tijd(t)));
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

%% Calculate and plot flux (OH-) at end time
% ddxC1(1)=(Coh(2,j)-Coh(1,j))/(dx); %forward fdd
% for v=2:n-1
%     ddxC1(v)=(Coh(v+1,j)-Coh(v-1,j))/(2*dx); %centred fdd
% end
% ddxC1(n)=(Coh(n,j)-Coh(n-1,j))/(dx); %backward fdd
% 
% ddxC2(1)=(Coh(n+2,j)-Coh(n+1,j))/(dx); %forward fdd
% for v=n+2:n+n-1
%     ddxC2(v-n)=(Coh(v+1,j)-Coh(v-1,j))/(2*dx); %centred fdd
% end
% ddxC2(n)=(Coh(n+n,j)-Coh(n+n-1,j))/(dx); %backward fdd
% 
% ddxC3(1)=(Coh(n+n+2,j)-Coh(n+n+1,j))/(dx); %forward fdd
% for v=n+n+2:n+n+n-1
%     ddxC3(v-n-n)=(Coh(v+1,j)-Coh(v-1,j))/(2*dx); %centred fdd
% end
% ddxC3(n)=(Coh(n+n+n,j)-Coh(n+n+n-1,j))/(dx); %backward fdd
% 
% Jdiff1=-D*ddxC1;
% Jdiff2=-Dm*ddxC2;
% Jdiff3=-D*ddxC3;
% Jdiff=[Jdiff1,Jdiff2,Jdiff3];
% 
% Jmig1=-zoh*f*D*Coh(1:n,j)'*E;
% Jmig2=-zoh*f*Dm*Coh(n+1:2*n,j)'*Em(end);
% Jmig3=-zoh*f*D*Coh(2*n+1:3*n,j)'*E;
% Jmig=[Jmig1,Jmig2,Jmig3];
% 
% Jtot=Jdiff+Jmig;
% 
% figure(3)
% hold on
% grid on
% xlabel('Position (m)')
% ylabel('OH^– flux (mole/m^2/s^2)')
% xticks(linspace(0,L*3,7))
% xlim([0,L*3])
% plot(x,Jtot,'-','DisplayName',sprintf('Total %g s',tijd(end)));
% plot(x,Jdiff,'--','DisplayName',sprintf('Diffusion %g s',tijd(end)));
% plot(x,Jmig,'--','DisplayName',sprintf('Migration %g s',tijd(end)));
% legend('AutoUpdate','off')
% rect=area([L,L*2],[max(get(gca,'YLim')),max(get(gca,'YLim'))],min(get(gca,'YLim')),...
%     'EdgeColor','none','ShowBaseLine','off','FaceColor','#999999','FaceAlpha',0.1);
% uistack(rect,'bottom')

%% Calculate and plot negative derivative of total flux (= accumulation of OH-)
% ddxJ1(1)=(Jtot(2)-Jtot(1))/(dx); %forward fdd
% for v=2:n-1
%     ddxJ1(v)=(Jtot(v+1)-Jtot(v-1))/(2*dx); %centred fdd
% end
% ddxJ1(n)=(Jtot(n)-Jtot(n-1))/(dx); %backward fdd
% 
% ddxJ2(1)=(Jtot(n+2)-Jtot(n+1))/(dx); %forward fdd
% for v=n+2:n+n-1
%     ddxJ2(v-n)=(Jtot(v+1)-Jtot(v-1))/(2*dx); %centred fdd
% end
% ddxJ2(n)=(Jtot(n+n)-Jtot(n+n-1))/(dx); %backward fdd
% 
% ddxJ3(1)=(Jtot(n+n+2)-Jtot(n+n+1))/(dx); %forward fdd
% for v=n+n+2:n+n+n-1
%     ddxJ3(v-n-n)=(Jtot(v+1)-Jtot(v-1))/(2*dx); %centred fdd
% end
% ddxJ3(n)=(Jtot(n+n+n)-Jtot(n+n+n-1))/(dx); %backward fdd
% 
% Acc=[-ddxJ1,-ddxJ2,-ddxJ3];
% 
% figure(4)
% hold on
% grid on
% xlabel('Position (m)')
% ylabel('OH^– accumulation (mole/m^3/s)')
% xticks(linspace(0,L*3,7))
% xlim([0,L*3])
% plot(x,Acc,'-','DisplayName',sprintf('%g s',tijd(end)));
% legend('AutoUpdate','off')
% rect=area([L,L*2],[max(get(gca,'YLim')),max(get(gca,'YLim'))],min(get(gca,'YLim')),...
%     'EdgeColor','none','ShowBaseLine','off','FaceColor','#999999','FaceAlpha',0.1);
% uistack(rect,'bottom')
% ylim([-50,200])