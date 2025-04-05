function [fel1,fel2,fel3] = RK4functionDLC(vgiss,thetagiss) % borde vi ha ut t_ut o u_ut också?

% Function som kommer användas för sekantmetoden, vi tar in en vinkel
% gissning i radianer som input, och outputen är felet från r värdet vi
% söker, t värderna och u värderna

% Begynnelsevärden
t(1) = 0 ; x(1)=0; xprim(1)=vgiss*cos(thetagiss); y(1)=-1.496e11; yprim(1) = vgiss*sin(thetagiss);
u=[x(1); xprim(1); y(1); yprim(1)];

h= 1000;
tslut=20000000; % cirka 240 dagar

for i=1:(tslut/h);
    %Nästa steg
    t(i+1) = t(i)+h;
    
    %RK4 
    k1=dudt(t(i), u(:,i));
    k2=dudt(t(i)+0.5*h,u(:,i)+0.5*k1*h);
    k3=dudt(t(i)+0.5*h,u(:,i)+0.5*k2*h);
    k4=dudt(t(i)+h,u(:,i)+k3*h);

    % Hittar r 
    r(i) = sqrt((u(1,i))^2+(u(3,i))^2);
    % Separerar x och y för grafen
    %xx(i) = u(1,i);
    %yy(i) = u(3,i);

    u(:,i+1)=u(:,i)+(h/6)*(k1+2*k2+2*k3+k4);

end

r(end);
minr = min(r);

fel1 = min(r)-2.7e10

% med 1000s per intervall är dag 210 iteration 18144
xC = u(1,18144);
yC = u(3,18144);  %IS THIS SHIT EVEN CORRECT?
fel2 = 0-xC
fel3 = 2.7e10-yC

%fel2 =sign([0-xC,2.7e10-yC])*abs([0-xC,2.7e10-yC])

% ÄR DET HÄR HIGHKEY 3 FEL, DÅ VI HAR XFELET OCH YFELET, vi testar
%isåfall behöver vi en tredje parameter
