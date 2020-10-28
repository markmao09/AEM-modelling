%% Initial guess of membrane concentration profile
%estimated with equilibrium and electroneutrality equations
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

%% Initial guess of membrane potential profile
%estimated with donnan potential equations
PotiemL=Phi0-R*T/z2/F*log(C2M/C20);
PotiemR=U+R*T/z2/F*log(C2N/C2M);
PotiemAVG=(PotiemL+PotiemR)/2; %average of both sides

%% Store initial guesses in matrix x0
x0(1:nnadl+2,1)=C10; %Conc. K+ in ADL
x0(1:nniem+2,2)=C1M; %Conc. K+ in IEM
x0(1:nncdl+2,3)=C1N; %Conc. K+ in CDL
x0(1:nnadl+2,4)=C20; %Conc. OH- in ADL
x0(1:nniem+2,5)=C2M; %Conc. OH- in IEM
x0(1:nncdl+2,6)=C2N; %Conc. OH- in CDL
x0(1:nnadl+2,7)=Phi0; %Pot. in ADL
x0(1:nniem+2,8)=PotiemAVG; %Pot. in IEM
x0(1:nncdl+2,9)=U; %Pot. in CDL