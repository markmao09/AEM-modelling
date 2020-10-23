%% Constants
F=96485; %C/mole
R=8.314; %J/(mole*K)

%% Discretisation
adl=1e-4; %m, length anode DL
iem=1e-4; %m, length ion exchange membrane
cdl=1e-4; %m, length cathode DL
N=adl+iem+cdl; %m, total length
dx=1e-6; %m, length scale discretisation
tnn=fix(N/dx); %--, total number nodes
enadl=fix(adl/dx); %--, end node ADL
eniem=enadl+fix(iem/dx); %--, end node IEM
encdl=eniem+fix(cdl/dx); %--, end node CDL
nnadl=enadl; %--, number nodes ADL
nniem=eniem-enadl; %--, number nodes IEM
nncdl=encdl-eniem; %--, number nodes CDL

%% Chemical and material properties
z1=1; %--
z2=-1; %--
zm=1; %--
D1=2e-9; %m2/s
D2=2e-9; %m2/s
Qfix=4800; %C/m3

%% Operating conditions
T=60+273.15; %K
Phi0=0; %V
U=-0.1; %V
C10=500; %mole/m3
C20=500; %mole/m3
C1N=200; %mole/m3
C2N=200; %mole/m3