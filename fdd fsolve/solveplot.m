clear
inputs
initialguess
options=optimoptions('fsolve','Display','iter-detailed');
[solution,fval,exitflag,output]=fsolve('equations',x0,options);


%% Generate plots
close all
x=linspace(0,tnn,tnn)*dx; %m
potential=[solution(2:end-1,7)',solution(2:end-1,8)',solution(2:end-1,9)']*1000; %V
figure(1)
top=roundn(max(potential),1)+10;
bot=roundn(min(potential),1)-10;
area([1e-4,2e-4],[top,top],bot,'EdgeColor','none','FaceColor','#00B0F0','FaceAlpha',0.05)
hold on
plot(x,potential,'-k')
grid on
xlabel('Position (m)')
ylabel('Potential (mV)')
xticks(linspace(0,tnn,7)*dx)
ylim([bot,top])

Ck=[solution(2:end-1,1)',solution(2:end-1,2)',solution(2:end-1,3)']/1000; %M=mole/L
Coh=[solution(2:end-1,4)',solution(2:end-1,5)',solution(2:end-1,6)']/1000; %M=mole/L
figure(2)
area([1e-4,2e-4],[5,5],0,'EdgeColor','none','FaceColor','#00B0F0','FaceAlpha',0.05)
hold on
oh=plot(x,Coh,'k-');
k=plot(x,Ck,'k--');
legend([oh,k],'OH^–','K^+')
grid on
xlabel('Position (m)')
ylabel('Concentration (M)')