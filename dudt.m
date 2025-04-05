function [uprim] = dudt(t,u) % u är (x,x',y,y')
%dudt, hittar uprim från u
uprim=[u(2);-(6.67430e-11)*(1.989e30)*u(1)/(((u(1)^2)+(u(3)^2))^(3/2));u(4);-(6.67430e-11)*(1.989e30)*u(3)/(((u(1)^2)+(u(3)^2))^(3/2))];
end
