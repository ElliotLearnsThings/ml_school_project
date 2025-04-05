clear all, close all, clc
% SEKANTMETODEN med 2 dimensioner
giss1 = [26000, -0.7];
giss0 = [27000, -0.75];

v_theta1 = giss1; v_theta0 = giss0;
% Initialising variables
dx = 10000
iter = 1;
f_1 = []
f_0 = []
f_3 = []
f_4 = []
while norm(dx)>1e-4
    % Sekantmetoden
    [f_1(1),f_1(2),f_1(3)] = RK4functionDLC(v_theta1(1),v_theta1(2));
    [f_0(1),f_0(2),f_0(3)] = RK4functionDLC(v_theta0(1),v_theta0(2));
    f_1
    
    [f_3(1),f_3(2),f_3(3)] = RK4functionDLC(v_theta1(1),v_theta0(2)); % ny v
    [f_4(1),f_4(2),f_4(3)] = RK4functionDLC(v_theta0(1),v_theta1(2)); % ny theta



    %J = [(f_1(1)-f_0(1))/(v_theta1(1)-v_theta0(1)), (f_1(1)-f_0(1))/(v_theta1(2)-v_theta0(2))  %dfel1/dv, dfel1/dtheta
    %     (f_1(2)-f_0(2))/(v_theta1(1)-v_theta0(1)), (f_1(2)-f_0(2))/(v_theta1(2)-v_theta0(2))] %dfel2/dv, dfel2/dtheta

    J = [(f_3(1)-f_0(1))/(v_theta1(1)-v_theta0(1)), (f_4(1)-f_0(1))/(v_theta1(2)-v_theta0(2))  %dfel1/dv, dfel1/dtheta
        (f_3(2)-f_0(2))/(v_theta1(1)-v_theta0(1)), (f_4(2)-f_0(2))/(v_theta1(2)-v_theta0(2))   %dfel2/dv, dfel2/dtheta
        (f_3(3)-f_0(3))/(v_theta1(1)-v_theta0(1)), (f_4(3)-f_0(3))/(v_theta1(2)-v_theta0(2))]  %dfel3/dv, dfel2/dtheta
    %det(J)
    dx = J\f_1'

    
    % Uppdaterar värden
    v_theta0 = v_theta1;
    v_theta1 = v_theta1-dx';  % ETT problem var dx vert, v_theta1 horiz
    iter = iter+1;
end
iter
format long
v_theta1
%  Verify result with DEL A
