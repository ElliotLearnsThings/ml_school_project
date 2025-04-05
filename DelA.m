clear all, close all, clc
format long

% Begynnelsevärden
t(1) = 0 ; x(1)=0; xprim(1)=25720*cos(-0.86649495308); y(1)=-1.496e11; yprim(1) = 25720*sin(-0.86649495308); % og vinkel pi/6, og v =15000
u=[x(1); xprim(1); y(1); yprim(1)];

h= 100
tslut=100000000

for i=1:(tslut/h)
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
    xx(i) = u(1,i);
    yy(i) = u(3,i);

    u(:,i+1)=u(:,i)+(h/6)*(k1+2*k2+2*k3+k4);
end

r(end)
minr = min(r)

% Med h=100, tslut = 10000000
% 1.565502420305857e+10
% Vi räknar var hundrade sekund, med en total tid på 115.74 dagar

plot(xx,yy)
hold on 
sunandmoon1 = [0,0]
sunandmoon2 = [0,-1.496e11]
scatter(sunandmoon1,sunandmoon2)
