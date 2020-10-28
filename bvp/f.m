function dydx=f(~,y,~) %equations being solved
inputs
dydx=zeros(6,1);

dydx(1)=y(2); %variable: d/dx(Phi)
dydx(2)=0; %Poisson/electroneutrality, variable: d2/dx2(Phi)

dydx(3)=y(4); %variable: d/dx(Coh)
dydx(4)=-F/R/T*z2*y(2)*y(4); %OH- nernst-planck/mass-balance, variable: d2/dx2(Coh)

dydx(5)=y(6); %variable: d/dx(Ck)
dydx(6)=-F/R/T*z1*y(2)*y(6); %K+ nernst-planck/mass-balance, variable: d2/dx2(Ck)
end