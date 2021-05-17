%  Arquivo para simular a resposta de perturbação no mancal magnético de 3
%  pólos.
%
%  As constantes utilizadas aqui são as mesmas do artigo do prof. Afonso
%  Simulação utilizada para gerar figuras com os erros de observação dos
%  estados com velocidade omega = 1050 rad/s
format SHORT
close all;
clc,clf; %limpar figuras 

tfinal = 0.2; %tempo de simulação
dt = 0.0001;
t = 0:dt:tfinal;

Iz = 0.0017;
J = 0.0592;
Ca = 0.0303;
gp = 19563.5; %I_base = 0.2325 A
gi = 42.07; %I_base = 0.2325 A
%kp15 = 23082; %I_base = 1 A
%ki15 = 9.23; %I_base = 1 A

b = 0.137; % valor do artigo = 0.137
d = 0.203; % valor do artigo = 0.203
r = 0.060;
m = 0.001; %massa de desbalanço
M = 0.5; %massa do rotor
q = 0.252;
omega = 9*356; % valor quase 3x maior = 1050 rad/s = 10026 rpm, 1rpm = pi/30 rad/s 
Ke = (1/2)*(1/J)*(b^2)*gp; % p/ mancal 3 pólos
Ge = (1/J)*[Ca omega*Iz;-omega*Iz Ca];

% matrizes do sistema dinamico x_dot = Ax + Bu p/ mancal 3 polos
A = [zeros(2) eye(2);Ke*eye(2) -Ge];
B = (1/J)*d*b*gi*[zeros(2);eye(2)];
C = [eye(2) zeros(2)];

%Matriz D de inserção de desbalanço no sistema
% com desbalanço: D = (1/J)*d*m*r*q*omega^2*[zeros(2); eye(2)];
% sem desbalanço : D = zeros(4,2);
%D = (1/J)*d*m*r*q*omega^2*[zeros(2); eye(2)];
D = zeros(4,2);

% Matrizes Q e R do regulador linear quadrático
Q = eye(4);
R = eye(2);

%ganhos de realimentação ótimos
[F,~,e]=lqr(A,B,Q,R); %centralizado ou convencional

% Verificar a resposta devido a alteração das condições iniciais
sys = ss(A-B*F,B,C,[]); %% ss(A-B*F,B,C,[])

x0 = [0.0002;-0.0002;0;0];

%% Observador
%https://www.mathworks.com/help/control/ref/place.html
%https://www.mathworks.com/help/control/ref/lqr.html
L = place(A',C',5*e).';

%% Simulação no simulink - modelo3
%resposta do mancal de 3 polos, correntes de controle, fx e fy
sim('modelo3',tfinal);
[y3,~,~] = initial(sys,x0,t);

%espessura da linha nos plots
line = 1.5;

%correntes de controle
id3 = (3/sqrt(3))*u(1,:);
id1 = u(2,:);

%Forças de relutância
fxl = 0.5*gp*y(1,:) + sqrt(3)/3*gi*id3;
fyl = 0.5*gp*y(2,:) + gi*id1;


%% Plots e informações com mudança de velocidade omega

%20 é pq divide por 2 e depois pelo número de degraus que se quer chegar até omegamax
n = (size(t,2)-1)/20;
nn = (size(t,2)-1)/2 + 1;
nnend = [2501:size(t,2)];
omegav = [0*ones(1,n) 100*ones(1,n) 200*ones(1,n) 300*ones(1,n) 400*ones(1,n) 500*ones(1,n) 600*ones(1,n) 700*ones(1,n) 800*ones(1,n) 900*ones(1,n) 1000*ones(1,size(nnend,2))]';

figure(6)
plot(t,y(1,:),t,xnl,'--','color',[0 0.4470 0.7410],'LineWidth',line)
hold on
plot(t,y(2,:),t,ynl,'--','color',[0.8500 0.3250 0.0980],'LineWidth',line)
grid on
tl = title('Output for LQR control input \boldmath${u}$=F\boldmath${x}$, for $\omega$=3204 rad/s','FontName','Times')
yl = ylabel('{$x_s,y_s$} (mm)','FontName','Times')
yt = get(gca,'YTick');
set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
% yticklabels({'-0.2','-0.15','-0.1','-0.05','0','0.05','0.1','0.15','0.2'})
xlabel('time (ms)','FontName','Times')
lg = legend('$x_s$ - linearized','$x_s$ - nonlinear','$y_s$ - linearized','$y_s$ - nonlinear','FontName','Times','Location','northeast')
xt = get(gca,'XTick');
set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)
xticklabels({'','50','100','150','200'})
set(tl,'Interpreter','latex');
set(yl,'Interpreter','latex');
set(lg,'Interpreter','latex');
print -depsc2 C:\Users\vinic\Documents\Mestrado\Pesquisa\artigo_ISMB17\ISMB17_artigo_Vinicius\lqr_lin_nonlin3204.eps

%% Simulação no simulink - modelo4
%observador oredm reduzida, omega = 1050 rad/s


A2 = mat2cell(A,[2 2],[2 2]);
Aaa = A2{1,1};
Aau = A2{1,2};
Aua = A2{2,1};
Auu = A2{2,2};

Ba = B(1:2,1:2);
Bu = B(3:4,1:2);

er = 5*[real(e(1)) real(e(1))];%autovalores referentes à parte não observavel do estado. 5x maior que o autovalor mais lento. 

Lr = place(Auu',Aau',er).';

A_hat = Auu - Lr*Aau;
B_hat = A_hat*Lr + Aua -Lr*Aaa;
C_hat = [zeros(2);eye(2)];
D_hat = [eye(2);Lr];
G = Bu - Lr*Ba;
sim('modelo4',tfinal);

figure(7)
subplot(2,1,1)
plot(t,soerror1(1,:),t,soerror1(2,:),'LineWidth',line)
hold on
%plot(t,y(2,:),t,ynl,'--','color',[0.8500 0.3250 0.0980],'LineWidth',line)
grid on
title('Observation error for the full order observer','FontName','Times')
yl1= ylabel('{$x_s,y_s$}(mm)','FontName','Times')
yt = get(gca,'YTick');
set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
% yticklabels({'-0.2','-0.15','-0.1','-0.05','0','0.05','0.1','0.15','0.2'})
xlabel('time (ms)','FontName','Times')
lg1 = legend('$x_s$','$y_s$','FontName','Times','Location','northeast')
xt = get(gca,'XTick');
set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)
set(yl1,'Interpreter','latex');
set(lg1,'Interpreter','latex');

subplot(2,1,2)
plot(t,soerror1(3,:),t,soerror1(4,:),'LineWidth',line)
hold on
%plot(t,y(2,:),t,ynl,'--','color',[0.8500 0.3250 0.0980],'LineWidth',line)
grid on
yl1= ylabel('{$\dot{x}_s,\dot{y}_s$}(mm/s)','FontName','Times')
yt = get(gca,'YTick');
set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
yticklabels({'-20','-10','0','10','20'})
xlabel('time (ms)','FontName','Times')
lg1 = legend('$\dot{x}_s$','$\dot{y}_s$','FontName','Times','Location','northeast')
ylim([-25 20]*10^(-3))
xt = get(gca,'XTick');
set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)
set(yl1,'Interpreter','latex');
set(lg1,'Interpreter','latex');
%print -depsc2 C:\Users\vinic\Documents\Mestrado\Pesquisa\artigo_ISMB17\ISMB17_artigo_Vinicius\obseoc2.eps

% Erro de estimação para Observador de Ordem Reduzida, omega = 1050 rad/s
figure(8)
subplot(2,1,1)
plot(t,soerror2(1,:),t,soerror2(2,:),'LineWidth',line)
hold on
grid on
title('Observation error for the reduced order observer','FontName','Times')
yl1= ylabel('{$x_s,y_s$}(mm)','FontName','Times')
yt = get(gca,'YTick');
set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
yticklabels({'-0.2','-0.1','0','0.1','0.2'})
xlabel('time (ms)','FontName','Times')
lg1 = legend('$x_s$','$y_s$','FontName','Times','Location','northeast')
xt = get(gca,'XTick');
set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)

set(yl1,'Interpreter','latex');
set(lg1,'Interpreter','latex');

subplot(2,1,2)
plot(t,soerror2(3,:),t,soerror2(4,:),'LineWidth',line)
hold on
%plot(t,y(2,:),t,ynl,'--','color',[0.8500 0.3250 0.0980],'LineWidth',line)
grid on
yl1= ylabel('{$\dot{x}_s,\dot{y}_s$}(mm/s)','FontName','Times')
yt = get(gca,'YTick');
set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
%yticklabels({'-20','-10','0','10','20'})
xlabel('time (ms)','FontName','Times')
lg1 = legend('$\dot{x}_s$','$\dot{y}_s$','FontName','Times','Location','northeast')
ylim([-50 50]*10^(-3))
xt = get(gca,'XTick');
set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)
set(yl1,'Interpreter','latex');
set(lg1,'Interpreter','latex');
%print -depsc2 C:\Users\vinic\Documents\Mestrado\Pesquisa\artigo_ISMB17\ISMB17_artigo_Vinicius\obseor2.eps

% Estados estimados no Observador de Ordem Reduzida
figure(9)
plot(t,x_hat_ro(3,:),t,x_hat_ro(4,:),'LineWidth',line)
hold on
%plot(t,y(2,:),t,ynl,'--','color',[0.8500 0.3250 0.0980],'LineWidth',line)
grid on
tl = title('Estimated states $\dot{x}_s$ and $\dot{y}_s$ for the reduced order observer','FontName','Times')
yl = ylabel('{$\dot{x}_s,\dot{y}_s$} (mm/s)','FontName','Times')

yt = get(gca,'YTick');
set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
% yticklabels({'-0.2','-0.15','-0.1','-0.05','0','0.05','0.1','0.15','0.2'})
xlabel('time (ms)','FontName','Times')
lg = legend('$\dot{x}_s$','$\dot{y}_s$','FontName','Times','Location','northeast')
xt = get(gca,'XTick');
set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)

set(yl,'Interpreter','latex');
set(lg,'Interpreter','latex');
set(tl,'Interpreter','latex');
%print -depsc2 C:\Users\vinic\Documents\Mestrado\Pesquisa\artigo_ISMB17\ISMB17_artigo_Vinicius\xsysroobs.eps

% Comparação dos estados estimados x_s_dot e y_s_dot para o Observador de ordem completa e de ordem reduzida
figure(10)
plot(t,soerror1(3,:),t,soerror1(4,:),'LineWidth',line)
hold on
plot(t,soerror2(3,:),t,soerror2(4,:),'LineWidth',line)
grid on
tl = title('Observation errors  - states $\dot{x}_s$ and $\dot{y}_s$','FontName','Times')
yl1= ylabel('{$\dot{x}_s,\dot{y}_s$}(mm/s)','FontName','Times')
yt = get(gca,'YTick');
set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
xticklabels({'','25','50','75','100'})
xlabel('time (ms)','FontName','Times')
lg1 = legend('$\dot{x}_s$ - full order','$\dot{y}_s$ - full order','$\dot{x}_s$ - reduced order','$\dot{y}_s$ - reduced order','FontName','Times','Location','northeast')
xt = get(gca,'XTick');
xlim([0 100]*10^(-3))
set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)

set(yl1,'Interpreter','latex');
set(lg1,'Interpreter','latex');
set(tl,'Interpreter','latex');
%print -depsc2 C:\Users\vinic\Documents\Mestrado\Pesquisa\artigo_ISMB17\ISMB17_artigo_Vinicius\oerrors.eps
