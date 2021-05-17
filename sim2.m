%  Arquivo para simular a resposta de perturbação no mancal magnético de 3
%  pólos.
%
%  As constantes utilizadas aqui são as mesmas do artigo do prof. Afonso
format shortG
close all;
clc,clf; %limpar figuras 

tfinal = 0.080; %tempo de simulação
dt = 0.0001;
t = 0:dt:tfinal;

Iz = 0.0017;
J = 0.0592;
Ca = 0.0303;
kp = 19563.5; %I_base = 0.2325 A
ki = 42.07; %I_base = 0.2325 A
%kp15 = 23082; %I_base = 1 A
%ki15 = 9.23; %I_base = 1 A

b = 0.137; %0.137
d = 0.203; %0.203
omega =356; % valor do artigo=356 rad/s = 3400 rpm, 1rpm = pi/30 rad/s 
Ke = (1/2)*(1/J)*(b^2)*kp; % p/ mancal 3 pólos
Ke8 = (1/J)*(b^2)*kp; % p/ mancal 8 pólos
Ke4 = (1/J)*(b^2)*2*kp;
Ge = (1/J)*[Ca omega*Iz;-omega*Iz Ca];


% matrizes do sistema dinamico x_dot = Ax + Bu p/ mancal 3 polos
A = [zeros(2) eye(2);Ke*eye(2) -Ge];
B = (1/J)*d*b*ki*[zeros(2);eye(2)];
C = [eye(2) zeros(2)];
D = zeros(2);


% matrizes do sistema dinamico x_dot = Ax + Bu p/ mancal 8 polos
A8 = [zeros(2) eye(2);Ke8*eye(2) -Ge];
%B8 = (1/J)*d*b*ki*[zeros(2);eye(2)];

% matrizes do sistema dinamico x_dot = Ax + Bu p/ mancal 4 polos
A4 = [zeros(2) eye(2);Ke4*eye(2) -Ge];
B4 = (1/J)*d*b*4*ki*[zeros(2);eye(2)];

% Matrizes Q e R do regulador linear quadrático
Q = eye(4);
R = eye(2);

%ganhos de realimentação ótimos
[F,P,e]=lqr(A,B,Q,R);
[F8,P8,e8]=lqr(A8,B,Q,R);
[F4,P4,e4]=lqr(A4,B4,Q,R);
%plot(real(e(1)), imag(e(1)),'g*',real(e(2)), imag(e(2)),'r+',real(e(3)), imag(e(3)),'bo',real(e(4)), imag(e(4)),'mx')

% Ganho do observador
%L = place(A',C',[3*real(e(1))+i*imag(e(1)) 3*real(e(2))+i*imag(e(2))...
 %                3*real(e(3))+i*imag(e(3)) 3*real(e(4))+i*imag(e(4))]).';

% Verificar a resposta devido a alteração das condições iniciais
sys = ss(A-B*F,B,C,[]);
sys8 = ss(A8-B*F8,B,C,[]);
sys4 = ss(A4-B4*F4,B4,C,[]);
x0 = [0.0002;-0.0002;0;0];

%% Figura 1
% Comparação resposta 3 polos com 8 e 4 polos.
%

%Simula no simulink
sim('modelo2',tfinal);

[y3,t3,x3] = initial(sys,x0,t);

[y8,t8,x8] = initial(sys8,x0,t);

[y4,t4,x4] = initial(sys4,x0,t);

figure(1);
subplot(2,1,1)
plot(t3,y3(:,1),t8,y8(:,1),'--',t4,y4(:,1),':','LineWidth',2)
grid on
title('Resposta à condição inicial para x_s ','FontName','Times')
ylabel('{x_s} (mm)','FontName','Times')
yt = get(gca,'YTick');
set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
xlabel('tempo (ms)','FontName','Times')
legend('3 Polos','8 Polos','4 Polos','FontName','Times')
xt = get(gca,'XTick');
set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)
%set(gca,'XTickLabel',[]) %remover números eixo x

subplot(2,1,2)
plot(t3,y3(:,2),t8,y8(:,2),'--',t4,y4(:,2),':','LineWidth',2)
grid on
title('Resposta à condição inicial para y_s','FontName','Times')
ylabel('{y_s} (mm)','FontName','Times')
yt = get(gca,'YTick');
set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
yticklabels({'-0.2','-0.15','-0.1','-0.05','0'})
xlabel('tempo (ms)','FontName','Times')
xt = get(gca,'XTick');
set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)
 
print -depsc2 C:\Users\vinic\OneDrive\Documentos\Mestrado\Pesquisa\Relatorios\initial.eps
%% Figura 2
% Estado real e estimado

tfinal2 = 0.04;
t2 = 0:dt:tfinal2;

sim('modelo2',tfinal2);


% figure(2);
% subplot(2,1,1)
% plot(t2,y(:,1),t2,y_hat(:,1),'--','LineWidth',2)
% grid on
% title('Estado real e estimado - x_s','FontName','Times')
% ylabel('{x_s} (mm)','FontName','Times')
% yt = get(gca,'YTick');
% set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
% xlabel('tempo (ms)','FontName','Times')
% legend('Estado real','Estado estimado','FontName','Times')
% xt = get(gca,'XTick');
% set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)
% 
% subplot(2,1,2)
% plot(t2,y(:,2),t2,y_hat(:,2),'--','LineWidth',2)
% grid on
% title('Estado real e estimado - y_s','FontName','Times')
% ylabel('{y_s} (mm)','FontName','Times')
% yt = get(gca,'YTick');
% set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
% xlabel('tempo (ms)','FontName','Times')
% xt = get(gca,'XTick');
% set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)
% 
% %print -depsc2 C:\Users\vinic\OneDrive\Documentos\Mestrado\Pesquisa\Relatorios\estados_obs.eps

%% Figura 3

[y8,t2,x8] = initial(sys8,x0,t2);

[y4,t2,x4] = initial(sys4,x0,t2);


% figure(3);
% subplot(2,1,1)
% plot(t2,y(:,1),t2,y8(:,1),'--',t2,y4(:,1),':','LineWidth',2)
% grid on
% title('Resposta à condição inicial para x_s  (Com Observador)','FontName','Times')
% ylabel('{x_s} (mm)','FontName','Times')
% yt = get(gca,'YTick');
% set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
% xlabel('tempo (ms)','FontName','Times')
% legend('3 Polos','8 Polos','4 Polos','FontName','Times')
% xt = get(gca,'XTick');
% set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)
% %set(gca,'XTickLabel',[])
% 
% subplot(2,1,2)
% plot(t2,y(:,2),t2,y8(:,2),'--',t2,y4(:,2),':','LineWidth',2)
% grid on
% title('Resposta à condição inicial para y_s (Com Observador)','FontName','Times')
% ylabel('{y_s} (mm)','FontName','Times')
% yt = get(gca,'YTick');
% set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
% %yticklabels({'-0.2','-0.15','-0.1','-0.05','0','0'})
% xlabel('tempo (ms)','FontName','Times')
% xt = get(gca,'XTick');
% set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)
%  
% %print -depsc2 C:\Users\vinic\OneDrive\Documentos\Mestrado\Pesquisa\Relatorios\initial_obs.eps

[Fdn, delta, Jc, Jd] = descentralizado(x0,A,B,Q,R);

%% Figura 4

% Controle Ótimo descentralizado
sim('modelo2',tfinal2);
y3c=y; % re3sposta do controle centralizado

%Fd=[5057.3 0 27.9 0;0 5057.3 0 27.9]; % impõe elementos nulos aos termos de acoplamento
F=Fdn;

sim('modelo2',tfinal2);
yd=y; % resposta do controle descentralizado

figure(4)
subplot(2,1,1)
plot(t2,y3c(:,1),t2,yd(:,1),'--','LineWidth',2)
grid on
title('Controle centralizado e descentralizado, {\omega} = 356 rad/s','FontName','Times')
ylabel('{x_s} (mm)','FontName','Times')
yt = get(gca,'YTick');
set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
xlabel('tempo (ms)','FontName','Times')
legend('Centralizado','Descentralizado','FontName','Times')
xt = get(gca,'XTick');
set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)
%set(gca,'XTickLabel',[])

subplot(2,1,2)
plot(t2,y3c(:,2),t2,yd(:,2),'--','LineWidth',2)
grid on
%title('Comparação entre controle centralizado e descentralizado','FontName','Times')
ylabel('{y_s} (mm)','FontName','Times')
yt = get(gca,'YTick');
set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
%yticklabels({'-0.2','-0.15','-0.1','-0.05','0','0'})
xlabel('tempo (ms)','FontName','Times')
xt = get(gca,'XTick');
set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)

%print -depsc2 C:\Users\vinic\OneDrive\Documentos\Mestrado\Pesquisa\Relatorios\cent_descent.eps
print -depsc2 C:\Users\vinic\OneDrive\Documentos\Mestrado\Pesquisa\artigo ISMB17\ISMB17 Abstract Vinicius\cent_descent.eps

% %Não vejo muita necessidade dessa figura 5. Renumerar depois. 
% figure (5)
% plot(yd(:,1),yd(:,2),y3c(:,1),y3c(:,2),'--','LineWidth',1)
% hold on
% plot(y(1,1),y(1,2),'+','Color','m','MarkerSize',6,'LineWidth',2) %estado inicial
% hold on
% plot(y(end,1),y(end,2),'o','MarkerSize',6,'LineWidth',2) % estado final
% viscircles([0 0],0.0005,'Color','k','LineStyle',':','LineWidth',1)
% xlabel('{x_s} (m)','FontName','Times')
% ylabel('{y_s} (m)','FontName','Times')
% title('Trajetórias do rotor','FontName','Times')
% legend('Descentralizado','Centralizado','Estado inicial','Estado em regime','FontName','Times')
% axis([-0.0006 0.0006 -0.0006 0.0006])

%print -depsc2 C:\Users\vinic\OneDrive\Documentos\Mestrado\Pesquisa\Relatorios\trajetorias.eps


%% Trajetória dos autovalores - Figura 5

for omega=0:40:4000
Ge = (1/J)*[0 omega*Iz;-omega*Iz 0];
A = [zeros(2) eye(2);Ke*eye(2) -Ge];
B = (1/J)*d*b*ki*[zeros(2);eye(2)];
C = [eye(2) zeros(2)];
D = zeros(2);

[~,~,e]=lqr(A,B,Q,R);

figure(5)
plot(real(e),imag(e),'.','Color','[0, 0.4470, 0.7410]')
grid on
hold on

title('Autovalores da malha fechada, Centralizado, \omega=356 a \omega=4000 rad/s','FontName','Times')
xlabel('Parte Real','FontName','Times','FontSize',12)
ylabel('Parte Imaginária','FontName','Times','FontSize',12)

end

figure(5)
omega = 356;
Ge = (1/J)*[0 omega*Iz;-omega*Iz 0];
A = [zeros(2) eye(2);Ke*eye(2) -Ge];
B = (1/J)*d*b*ki*[zeros(2);eye(2)];
C = [eye(2) zeros(2)];
D = zeros(2);

disp('Autovalores centralizado omega = 356')
[~,~,e]=lqr(A,B,Q,R)

plot(real(e),imag(e),'+','Color','[0.4660, 0.6740, 0.1880]','LineWidth',3,'MarkerSize',14)

omega = 4000;
Ge = (1/J)*[0 omega*Iz;-omega*Iz 0];
A = [zeros(2) eye(2);Ke*eye(2) -Ge];
B = (1/J)*d*b*ki*[zeros(2);eye(2)];
C = [eye(2) zeros(2)];
D = zeros(2);

[F,P,e]=lqr(A,B,Q,R);

disp('Autovalores centralizado omega = 3000')
e
plot(real(e),imag(e),'x','Color','[0.8500, 0.3250, 0.0980]','LineWidth',3,'MarkerSize',14)

h=zeros(2,1);
%h(1)=plot(NaN,NaN,'+','Color','r','LineWidth',3,'MarkerSize',14);
h(1)=plot(NaN,NaN,'+','Color','[0.4660, 0.6740, 0.1880]','LineWidth',3,'MarkerSize',14);
h(2)=plot(NaN,NaN,'x','Color','[0.8500, 0.3250, 0.0980]','LineWidth',3,'MarkerSize',14);
legend(h,'\omega = 356 rad/s','\omega = 4000 rad/s','Location','northwest','FontName','Times')
%axis([-190 -163 -80 80])

print -depsc2 C:\Users\vinic\OneDrive\Documentos\Mestrado\Pesquisa\Relatorios\autovalores.eps

%% Figura 6 - similar à fig. 5, descentralizado

for omega=0:20:4000
Ge = (1/J)*[0 omega*Iz;-omega*Iz 0];
A = [zeros(2) eye(2);Ke*eye(2) -Ge];
B = (1/J)*d*b*ki*[zeros(2);eye(2)];
C = [eye(2) zeros(2)];
D = zeros(2);

[Fdn,~, ~,~] = descentralizado(x0,A,B,Q,R);

e=eig(A-B*Fdn);

figure(6)
plot(real(e),imag(e),'.','Color','[0, 0.4470, 0.7410]')
grid on
hold on

title('Autovalores da malha fechada, Descentralizado, \omega=356 a \omega=4000 rad/s','FontName','Times')
xlabel('Parte Real','FontName','Times','FontSize',12)
ylabel('Parte Imaginária','FontName','Times','FontSize',12)

end

figure(6)
omega = 356;
Ge = (1/J)*[0 omega*Iz;-omega*Iz 0];
A = [zeros(2) eye(2);Ke*eye(2) -Ge];
B = (1/J)*d*b*ki*[zeros(2);eye(2)];
C = [eye(2) zeros(2)];
D = zeros(2);

[Fdn,~, ~,~] = descentralizado(x0,A,B,Q,R);

disp('Autovalores descentralizado omega = 356')
e=eig(A-B*Fdn)

plot(real(e),imag(e),'+','Color','[0.4660, 0.6740, 0.1880]','LineWidth',3,'MarkerSize',14)

omega = 4000;
Ge = (1/J)*[0 omega*Iz;-omega*Iz 0];
A = [zeros(2) eye(2);Ke*eye(2) -Ge];
B = (1/J)*d*b*ki*[zeros(2);eye(2)];
C = [eye(2) zeros(2)];
D = zeros(2);

[Fdn,~, ~,~] = descentralizado(x0,A,B,Q,R);

disp('Autovalores descentralizado omega = 4000')
e = eig(A-B*Fdn)

plot(real(e),imag(e),'x','Color','[0.8500, 0.3250, 0.0980]','LineWidth',3,'MarkerSize',14)

h=zeros(2,1);
%h(1)=plot(NaN,NaN,'+','Color','r','LineWidth',3,'MarkerSize',14);
h(1)=plot(NaN,NaN,'+','Color','[0.4660, 0.6740, 0.1880]','LineWidth',3,'MarkerSize',14);
h(2)=plot(NaN,NaN,'x','Color','[0.8500, 0.3250, 0.0980]','LineWidth',3,'MarkerSize',14);
legend(h,'\omega = 356 rad/s','\omega = 4000 rad/s','Location','northeast','FontName','Times')
axis([-80 50 -150 150])

print -depsc2 C:\Users\vinic\OneDrive\Documentos\Mestrado\Pesquisa\Relatorios\autovalores_descent.eps

%% Figura 6 - similar a fig.4, mas omega = 8000 rad/s
%Figura não utilizada na dissertação


% omega = 8000;
% Ge = (1/J)*[0 omega*Iz;-omega*Iz 0];
% A = [zeros(2) eye(2);Ke*eye(2) -Ge];
% B = (1/J)*d*b*ki*[zeros(2);eye(2)];
% C = [eye(2) zeros(2)];
% D = zeros(2);
% 
% [F,P,e]=lqr(A,B,Q,R);
% 
% % Controle Ótimo descentralizado, omega = 8000 rad/s
% sim('modelo1',tfinal2);
% y3c=y; % re3sposta do controle centralizado
% 
% %Fd=[5057.3 0 27.9 0;0 5057.3 0 27.9]; % impõe elementos nulos aos termos de acoplamento
% F=[F(1,1) 0 F(1,3) 0;0 F(2,2) 0 F(2,4)];
% 
% sim('modelo1',tfinal2);
% yd=y; % resposta do controle descentralizado
% 
% figure(6)
% subplot(2,1,1)
% plot(t2,y3c(:,1),t2,yd(:,1),'--','LineWidth',2)
% grid on
% title('Controle centralizado e descentralizado, {\omega} = 356 rad/s','FontName','Times')
% ylabel('{x_s} (mm)','FontName','Times')
% yt = get(gca,'YTick');
% set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
% xlabel('tempo (ms)','FontName','Times')
% legend('Centralizado','Descentralizado','FontName','Times')
% xt = get(gca,'XTick');
% set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)
% %set(gca,'XTickLabel',[])
% 
% subplot(2,1,2)
% plot(t2,y3c(:,2),t2,yd(:,2),'--','LineWidth',2)
% grid on
% %title('Comparação entre controle centralizado e descentralizado','FontName','Times')
% ylabel('{y_s} (mm)','FontName','Times')
% yt = get(gca,'YTick');
% set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
% %yticklabels({'-0.2','-0.15','-0.1','-0.05','0','0'})
% xlabel('tempo (ms)','FontName','Times')
% xt = get(gca,'XTick');
% set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)
% 
% %print -depsc2 C:\Users\vinic\OneDrive\Documentos\Mestrado\Pesquisa\Relatorios\cent_descent_8000.eps

%% Simulação desbalanço de massa
% % Figura 7
% 
% omega = 360;
% Ge = (1/J)*[0 omega*Iz;-omega*Iz 0];
% A = [zeros(2) eye(2);Ke*eye(2) -Ge];
% B = (1/J)*d*b*ki*[zeros(2);eye(2)];
% C = [eye(2) zeros(2)];
% D = zeros(2);
% 
% [F,~,e]=lqr(A,B,Q,R);
% %Fd=[5057.3 0 27.9 0;0 5057.3 0 27.9]; % impõe elementos nulos aos termos de acoplamento
% F=[F(1,1) 0 F(1,3) 0;0 F(2,2) 0 F(2,4)];
% 
% sim('modelo_pid',t2);
% figure(7)
% plot(y1(:,1),y1(:,2),'LineWidth',0.5)
% axis([-1e-4 3e-4 -2e-4 1e-4])
% hold on
% 
% Ke8 = (1/J)*(b^2)*kp; % p/ mancal 8 pólos
% Ke4 = (1/J)*(b^2)*4*kp;
% Ge = (1/J)*[0 omega*Iz;-omega*Iz 0];
% 
% % matrizes do sistema dinamico x_dot = Ax + Bu p/ mancal 8 polos
% A = [zeros(2) eye(2);Ke8*eye(2) -Ge];
% %B8 = (1/J)*d*b*ki*[zeros(2);eye(2)];
% sim('modelo_pid',t2);
% plot(y1(:,1),y1(:,2),'LineWidth',0.5)
% hold on
% 
% % matrizes do sistema dinamico x_dot = Ax + Bu p/ mancal 4 polos
% A = [zeros(2) eye(2);Ke4*eye(2) -Ge];
% B = (1/J)*d*b*4*ki*[zeros(2);eye(2)];
% 
% sim('modelo_pid',t2);
% plot(y1(:,1),y1(:,2),'LineWidth',0.5)

%% Figura 8
% %Igual a fig. 4, mas com velocidade 8000 
% 
% omega = 4000;
% Ge = (1/J)*[0 omega*Iz;-omega*Iz 0];
% A = [zeros(2) eye(2);Ke*eye(2) -Ge];
% B = (1/J)*d*b*ki*[zeros(2);eye(2)];
% C = [eye(2) zeros(2)];
% D = zeros(2);
% F = lqr(A,B,Q,R);
% 
% % Controle Ótimo descentralizado
% sim('modelo1',tfinal2);
% y3c=y; % resposta do controle centralizado
% 
% [Fdn,~, ~,~] = descentralizado(x0,A,B,Q,R);
% 
% F = Fdn;
% 
% sim('modelo1',tfinal2);
% yd=y; % resposta do controle descentralizado
% 
% figure(8)
% subplot(2,1,1)
% plot(t2,y3c(:,1),t2,yd(:,1),'--','LineWidth',2)
% grid on
% title('Controle centralizado e descentralizado, {\omega} = 4000 rad/s','FontName','Times')
% ylabel('{x_s} (mm)','FontName','Times')
% yt = get(gca,'YTick');
% set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
% xlabel('tempo (ms)','FontName','Times')
% legend('Centralizado','Descentralizado','FontName','Times')
% xt = get(gca,'XTick');
% set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)
% %set(gca,'XTickLabel',[])
% 
% subplot(2,1,2)
% plot(t2,y3c(:,2),t2,yd(:,2),'--','LineWidth',2)
% grid on
% %title('Comparação entre controle centralizado e descentralizado','FontName','Times')
% ylabel('{y_s} (mm)','FontName','Times')
% yt = get(gca,'YTick');
% set(gca,'fontsize',12,'YTick',yt,'YTickLabel',yt*1000)
% %yticklabels({'-0.2','-0.15','-0.1','-0.05','0','0'})
% xlabel('tempo (ms)','FontName','Times')
% xt = get(gca,'XTick');
% set(gca,'fontsize',12,'XTick',xt,'XTickLabel',xt*1000)%}
% %print -depsc2 C:\Users\vinic\OneDrive\Documentos\Mestrado\Pesquisa\Relatorios\cent_descent_4000.eps