function res=bc(yl,yr) %boundary conditions
%{
solving 3 2° DE's in 3 regions requires 18 BC's.
there are 3 regions and 4 boundaries: b1 r1 b2 r2 b3 r3 b4
syntax is: yl(n,r) and yr(n,r) are left and right edges of regions r=1,2,3 for variables n=1...6
the n's are in the order of the differential equations listed in f.m
n=1 is Phi
n=2 is d/dx(Phi)
n=3 is Coh
n=4 is d/dx(Coh)
n=5 is Ck
n=6 is d/dx(Ck)
%}
inputs

res=[yl(1,1)-Phi0 %fixed potential @b1
     yl(3,1)-C20 %fixed OH- concentration @b1
     yl(5,1)-C10 %fixed K+ concentration @b1
     yr(1,3)-U %fixed potential @b4
     yr(3,3)-C2N %fixed OH- concentration @b4
     yr(5,3)-C1N %fixed K+ concentration @b4
     
     -D2*yl(4,2)-F/R/T*z2*D2*yl(3,2)*yl(2,2)-...
     (-D2*yr(4,1)-F/R/T*z2*D2*yr(3,1)*yr(2,1)) %continuity molar flux OH- @b2
     -D2*yl(4,3)-F/R/T*z2*D2*yl(3,3)*yl(2,3)-...
     (-D2*yr(4,2)-F/R/T*z2*D2*yr(3,2)*yr(2,2)) %continuity molar flux OH- @b3
     -D1*yl(6,2)-F/R/T*z1*D1*yl(5,2)*yl(2,2)-...
     (-D1*yr(6,1)-F/R/T*z1*D1*yr(5,1)*yr(2,1)) %continuity molar flux K+ @b2
     -D1*yl(6,3)-F/R/T*z1*D1*yl(5,3)*yl(2,3)-...
     (-D1*yr(6,2)-F/R/T*z1*D1*yr(5,2)*yr(2,2)) %continuity molar flux K+ @b3
     
     z1*yl(5,2)+z2*yl(3,2)+zm*Qfix %electroneutrality membrane @b2
     z1*yr(5,1)+z2*yr(3,1) %electroneutrality adl @b2
     z1*yr(5,2)+z2*yr(3,2)+zm*Qfix %electroneutrality membrane @b3
     z1*yl(5,3)+z2*yl(3,3) %electroneutrality cdl @b3

     (yr(1,1)-yl(1,2))-R*T/F/z2*log(yl(3,2)/yr(3,1)) %donnan potential @b2
     (yr(1,2)-yl(1,3))-R*T/F/z2*log(yl(3,3)/yr(3,2)) %donnan potential @b3
     R*T/F/z2*log(yl(3,2)/yr(3,1))-R*T/F/z1*log(yl(5,2)/yr(5,1)) %donnan potential @b2
     R*T/F/z2*log(yl(3,3)/yr(3,2))-R*T/F/z1*log(yl(5,3)/yr(5,2))]; %donnan potential @b3
end