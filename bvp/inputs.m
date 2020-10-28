%% Constants
F=96485; %C/mole
R=8.314; %J/(mole*K)

%% Region partitioning
adl=1e-4; %m, length anode DL
iem=1e-4; %m, length ion exchange membrane
cdl=1e-4; %m, length cathode DL
enadl=adl; %m, end position ADL
eniem=adl+iem; %m, end position IEM
encdl=adl+iem+cdl; %m, end position CDL

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
C1N=500; %mole/m3
C2N=500; %mole/m3