function ddt=nonlinear(t,X,id3,id1)
% Modelagem da dinâmica não linear do MM de 3 polos
% pág.7 - Optimal Design of a 3-pole AMB.
% Leis de Newton
% fx = m*x_ddot -> x_ddot = (1/m)*fx
% fy-mg = m*y_ddot -> y_ddot = (1/m)fy -g

%fx e fy estão definidos na p.47 da minha dissertação
mu = 4*pi*10^(-7); %permeabilidade mag. do vácuo
A = 4*10^(-4); %área superficial do polo
h = 0.5*10^(-3); %entreferro
n = 300; %número de espiras da bobina
i0 = 0.2325; %corrente de base, na dissertação é iB1
g = 9.81;
m = 0.5;
%x = X(1);
%y = X(2);

%N1 = 6*(i0)^2*h*x + 3*(i0)^2*x*y + 4*sqrt(3)*(i0)*h^2*(id3) + 12*i0*h*(id1)*x + 6*i0*(id1)*x*y + sqrt(3)*(i0)*(id3)*x^2 -sqrt(3)*(i0)*(id3)*y^2 + 4*sqrt(3)*h^2*(id1)*(id3) + 6*h*(id1)^2*x + 2*h*(id3)^2*x + 3*(id1)^2*x*y + sqrt(3)*(id1)*(id3)*x^2 - sqrt(3)*(id1)*(id3)*y^2 - (id3)^2*x*y;

%N2 = 12*(i0)^2*h^2 + 12*(i0)^2*h*y - 3*(i0)^2*x^2 + 3*(i0)^2*y^2 + 24*(i0)*h^2*(id1) + 24*(i0)*h*y*(id1) - 6*(i0)*(id1)*x^2 + 6*(i0)*(id1)*y^2 + 4*sqrt(3)*i0*(id3)*x*y + 12*h^2*(id1)^2 - 4*h^2*(id3)^2 + 12*h*(id1)^2*y + 4*h*(id3)^2*y - 3*(id1)^2*x^2 + 3*(id1)^2*y^2 + 4*sqrt(3)*(id1)*(id3)*x*y + (id3)^2*x^2 -(id3)^2*y^2;

%DELTA = (x^2 + y^2 - 4*h^2);

%fx = (4/3)*mu*A*n^2*N1/DELTA^2;
%fy = (2/3)*mu*A*n^2*N2/DELTA^2;

%x = [x; y; xdot; ydot] = [x(1); x(2); x(3); x(4)] (equação 4.57 da
%dissertação)

%xs = xnl(1);
%ys = xnl(2);
%xs_dot = xnl(3);
%ys_dot = xnl(4);

ddt = zeros(size(X));
ddt(1) = X(3); %xs_dot =ddt (xs)
ddt(2) = X(4); %ys_dot =ddt (ys)
%ddt(3) = (1/m)*fx; % x_ddot = (1/m)*fx
%ddt(4) = (1/m)*fy - g; % x_ddot = (1/m)*fy - g
ddt(3) = (1/m)* (4/3)*mu*A*n^2*(6*(i0)^2*h*X(1) + 3*(i0)^2*X(1)*X(2) + 4*sqrt(3)*(i0)*h^2*(id3) + 12*i0*h*(id1)*X(1) + 6*i0*(id1)*X(1)*X(2) + sqrt(3)*(i0)*(id3)*X(1)^2 -sqrt(3)*(i0)*(id3)*X(2)^2 + 4*sqrt(3)*h^2*(id1)*(id3) + 6*h*(id1)^2*X(1) + 2*h*(id3)^2*X(1) + 3*(id1)^2*X(1)*X(2) + sqrt(3)*(id1)*(id3)*X(1)^2 - sqrt(3)*(id1)*(id3)*X(2)^2 - (id3)^2*X(1)*X(2))/(X(1)^2 + X(2)^2 - 4*h^2)^2;
ddt(4) = - g + (1/m)*(2/3)*mu*A*n^2*(12*(i0)^2*h^2 + 12*(i0)^2*h*X(2) - 3*(i0)^2*X(1)^2 + 3*(i0)^2*X(2)^2 + 24*(i0)*h^2*(id1) + 24*(i0)*h*X(2)*(id1) - 6*(i0)*(id1)*X(1)^2 + 6*(i0)*(id1)*X(2)^2 + 4*sqrt(3)*i0*(id3)*X(1)*X(2) + 12*h^2*(id1)^2 - 4*h^2*(id3)^2 + 12*h*(id1)^2*X(2) + 4*h*(id3)^2*X(2) - 3*(id1)^2*X(1)^2 + 3*(id1)^2*X(2)^2 + 4*sqrt(3)*(id1)*(id3)*X(1)*X(2) + (id3)^2*X(1)^2 -(id3)^2*X(2)^2)/(X(1)^2 + X(2)^2 - 4*h^2)^2 ;
end