function output=equations(x0)
inputs

%Define variables
C1adl(:,1)=x0(:,1);
C1iem(:,1)=x0(:,2);
C1cdl(:,1)=x0(:,3);
C2adl(:,1)=x0(:,4);
C2iem(:,1)=x0(:,5);
C2cdl(:,1)=x0(:,6);
Phiadl(:,1)=x0(:,7);
Phiiem(:,1)=x0(:,8);
Phicdl(:,1)=x0(:,9);

%Poisson with electroneutrality everywhere
ddxPhiadl=(Phiadl(2:nnadl+2)-Phiadl(1:nnadl+1))/dx; %FFD
Poi1=(ddxPhiadl(2:nnadl+1)-ddxPhiadl(1:nnadl))/dx; %FFD
ddxPhiiem=(Phiiem(2:nniem+2)-Phiiem(1:nniem+1))/dx; %FFD
Poi2=(ddxPhiiem(2:nniem+1)-ddxPhiiem(1:nniem))/dx; %FFD
ddxPhicdl=(Phicdl(2:nncdl+2)-Phicdl(1:nncdl+1))/dx; %FFD
Poi3=(ddxPhicdl(2:nncdl+1)-ddxPhicdl(1:nncdl))/dx; %FFD
    
%Define average concentrations to get fsolve to converge, interesting
C1adlminus1=(C1adl(2:nnadl+1)+C1adl(1:nnadl))/2;
C1adlplus1=(C1adl(3:nnadl+2)+C1adl(2:nnadl+1))/2;
C1iemminus1=(C1iem(2:nniem+1)+C1iem(1:nniem))/2;
C1iemplus1=(C1iem(3:nniem+2)+C1iem(2:nniem+1))/2;
C1cdlminus1=(C1cdl(2:nncdl+1)+C1cdl(1:nncdl))/2;
C1cdlplus1=(C1cdl(3:nncdl+2)+C1cdl(2:nncdl+1))/2;
C2adlminus1=(C2adl(2:nnadl+1)+C2adl(1:nnadl))/2;
C2adlplus1=(C2adl(3:nnadl+2)+C2adl(2:nnadl+1))/2;
C2iemminus1=(C2iem(2:nniem+1)+C2iem(1:nniem))/2;
C2iemplus1=(C2iem(3:nniem+2)+C2iem(2:nniem+1))/2;
C2cdlminus1=(C2cdl(2:nncdl+1)+C2cdl(1:nncdl))/2;
C2cdlplus1=(C2cdl(3:nncdl+2)+C2cdl(2:nncdl+1))/2;

%Nernst-planck and mass balance for ion 1 (K+) everywhere
J1adlmin1=-D1*(C1adl(2:nnadl+1)-C1adl(1:nnadl))/dx-z1*F/R/T*D1.*C1adlminus1(1:nnadl).*...
    (Phiadl(2:nnadl+1)-Phiadl(1:nnadl))/dx;
J1adlplu1=-D1*(C1adl(3:nnadl+2)-C1adl(2:nnadl+1))/dx-z1*F/R/T*D1.*C1adlplus1(1:nnadl).*...
    (Phiadl(3:nnadl+2)-Phiadl(2:nnadl+1))/dx;
J1iemmin1=-D1*(C1iem(2:nniem+1)-C1iem(1:nniem))/dx-z1*F/R/T*D1.*C1iemminus1(1:nniem).*...
    (Phiiem(2:nniem+1)-Phiiem(1:nniem))/dx;
J1iemplu1=-D1*(C1iem(3:nniem+2)-C1iem(2:nniem+1))/dx-z1*F/R/T*D1.*C1iemplus1(1:nniem).*...
    (Phiiem(3:nniem+2)-Phiiem(2:nniem+1))/dx;
J1cdlmin1=-D1*(C1cdl(2:nncdl+1)-C1cdl(1:nncdl))/dx-z1*F/R/T*D1.*C1cdlminus1(1:nncdl).*...
    (Phicdl(2:nncdl+1)-Phicdl(1:nncdl))/dx;
J1cdlplu1=-D1*(C1cdl(3:nncdl+2)-C1cdl(2:nncdl+1))/dx-z1*F/R/T*D1.*C1cdlplus1(1:nncdl).*...
    (Phicdl(3:nncdl+2)-Phicdl(2:nncdl+1))/dx;
Mb1=J1adlplu1-J1adlmin1;
Mb2=J1iemplu1-J1iemmin1;
Mb3=J1cdlplu1-J1cdlmin1;

%Nernst-planck and mass balance for ion 2 (OH-) everywhere
J2adlmin1=-D2*(C2adl(2:nnadl+1)-C2adl(1:nnadl))/dx-z2*F/R/T*D2.*C2adlminus1(1:nnadl).*...
    (Phiadl(2:nnadl+1)-Phiadl(1:nnadl))/dx;
J2adlplu1=-D2*(C2adl(3:nnadl+2)-C2adl(2:nnadl+1))/dx-z2*F/R/T*D2.*C2adlplus1(1:nnadl).*...
    (Phiadl(3:nnadl+2)-Phiadl(2:nnadl+1))/dx;
J2iemmin1=-D2*(C2iem(2:nniem+1)-C2iem(1:nniem))/dx-z2*F/R/T*D2.*C2iemminus1(1:nniem).*...
    (Phiiem(2:nniem+1)-Phiiem(1:nniem))/dx;
J2iemplu1=-D2*(C2iem(3:nniem+2)-C2iem(2:nniem+1))/dx-z2*F/R/T*D2.*C2iemplus1(1:nniem).*...
    (Phiiem(3:nniem+2)-Phiiem(2:nniem+1))/dx;
J2cdlmin1=-D2*(C2cdl(2:nncdl+1)-C2cdl(1:nncdl))/dx-z2*F/R/T*D2.*C2cdlminus1(1:nncdl).*...
    (Phicdl(2:nncdl+1)-Phicdl(1:nncdl))/dx;
J2cdlplu1=-D2*(C2cdl(3:nncdl+2)-C2cdl(2:nncdl+1))/dx-z2*F/R/T*D2.*C2cdlplus1(1:nncdl).*...
    (Phicdl(3:nncdl+2)-Phicdl(2:nncdl+1))/dx;
Mb4=J2adlplu1-J2adlmin1;
Mb5=J2iemplu1-J2iemmin1;
Mb6=J2cdlplu1-J2cdlmin1;

%Boundary condition: continuity of molar flux across adl-iem and iem-cdl interfaces
Con1=J1iemmin1(1)-J1adlplu1(end);
Con2=J1cdlmin1(1)-J1iemplu1(end);
Con3=J2iemmin1(1)-J2adlplu1(end);
Con4=J2cdlmin1(1)-J2iemplu1(end);

%Boundary condition: electroneutrality at adl-iem and iem-cdl interfaces
En1=z2*C2adl(end)+z1*C1adl(end);
En2=z2*C2iem(1)+z1*C1iem(1)+zm*Qfix;
En3=z2*C2iem(end)+z1*C1iem(end)+zm*Qfix;
En4=z2*C2cdl(1)+z1*C1cdl(1);

%Boundary condition: donnan potential at adl-iem and iem-cdl interfaces
Don1=(Phiadl(nnadl+2)-Phiiem(1))-R*T/F/z2*log(C2iem(1)/C2adl(nnadl+2));
Don2=(Phiiem(nniem+2)-Phicdl(1))-R*T/F/z2*log(C2cdl(1)/C2iem(nniem+2));
Don3=R*T/F/z2*log(C2iem(1)/C2adl(nnadl+2))-R*T/F/z1*log(C1iem(1)/C1adl(nnadl+2)); %interesting
Don4=R*T/F/z2*log(C2cdl(1)/C2iem(nniem+2))-R*T/F/z1*log(C1cdl(1)/C1iem(nniem+2));

%Boundary conditions: fixed values at 0 and N edges
BC1=C1adl(1)-C10;
BC2=C2adl(1)-C20;
BC3=Phiadl(1)-Phi0;
BC4=C1cdl(nncdl+2)-C1N;
BC5=C2cdl(nncdl+2)-C2N;
BC6=Phicdl(nncdl+2)-U;

output=[Poi1;Poi2;Poi3;Mb1;Mb2;Mb3;Mb4;Mb5;Mb6;Con1;Con2;Con3;Con4;En1;En2;En3;En4;...
        Don1;Don2;Don3;Don4;BC1;BC2;BC3;BC4;BC5;BC6];
end