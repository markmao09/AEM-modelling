clear
inputs

%% Initial guess of membrane concentration and potential profiles
%estimated with equilibrium and electroneutrality equations:
    a=1;
    b=zm*Qfix/z1;
    c=C10*C20*z2/z1;
    r=roots([a b c]);
    C1m0=r(r>0);
    C2m0=(-zm*Qfix-z1*C1m0)/(z2);
    c=C1N*C2N*z2/z1;
    r=roots([a b c]);
    C1mN=r(r>0);
    C2mN=(zm*Qfix+z1*C1mN)/(-z2);
    C1M=(C1m0+C1mN)/2; %average of both sides
    C2M=(C2m0+C2mN)/2; %average of both sides
%estimated with donnan potential equations:
    PotiemL=Phi0-R*T/z2/F*log(C2M/C20);
    PotiemR=U+R*T/z2/F*log(C2N/C2M);
    PotiemAVG=(PotiemL+PotiemR)/2; %average of both sides

%% Solve with bvp5c
xmesh=[linspace(0,enadl,51) linspace(enadl,eniem,51) linspace(eniem,encdl,51)];
yinit=[PotiemAVG; PotiemAVG; C2M; C2M; C1M; C1M]; %mediocre guesses...
sol=bvpinit(xmesh,yinit);
sol=bvp5c(@f,@bc,sol);

%% Generate plots
close all
figure(1)
top=max(sol.y(1,:)*1000);
bot=min(sol.y(1,:)*1000);
area([1e-4,2e-4],[top,top],bot,'EdgeColor','none','FaceColor','#00B0F0','FaceAlpha',0.05)
hold on
plot(sol.x,sol.y(1,:)*1000,'-b')
grid on
xlabel('Position (m)')
ylabel('Potential (mV)')
ylim([bot,top])
xticks(linspace(0,encdl,7))
xlim([0,encdl])

figure(2)
top=max([sol.y(3,:)/1000,sol.y(5,:)/1000]);
bot=min([sol.y(3,:)/1000,sol.y(5,:)/1000]);
area([1e-4,2e-4],[top,top],bot,'EdgeColor','none','FaceColor','#00B0F0','FaceAlpha',0.05)
hold on
oh=plot(sol.x,sol.y(3,:)/1000,'-b');
k=plot(sol.x,sol.y(5,:)/1000,'--b');
legend([oh,k],'OH^–','K^+')
grid on
xlabel('Position (m)')
ylabel('Concentration (M)')
ylim([bot,top])
xticks(linspace(0,encdl,7))
xlim([0,encdl])