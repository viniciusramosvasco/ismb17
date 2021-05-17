% Função que retorna o ganho de realimentação descentralizado
%
% Baseado na dissertação Leonardo Sodré Rodrigues-Controle Ótimo Descentralizado a 2 parâmetros para MMM 
% e também a tese de PhD de Bleuler 1984.
%
% Fdn - ganho descentralizado após última iteração
% Fd - ganho descentralizado na 1ª iteração
% F - ganho centralizado
%
% delta - critério de parada, diferença entre norma matricial de duas
% iterações da matriz de ganho descentralizado
%
function [Fdn, delta, Jc, Jd]=descentralizado(x0,A,B,Q,R)
X0 = x0 * x0';

F = lqr(A,B,Q,R);

A0 = A-B*F;
B1 = B(:,1);
B2 = B(:,2);

C1 = [1 0 0 0;0 0 1 0];
C2 = [0 1 0 0;0 0 0 1];

%1º passo - usando matriz F centralizada p/ o algoritmo
X  = lyap(A0,X0);
N = Q + F(1,:)'*R(1,1)*F(1,:) + F(2,:)'*R(2,2)*F(2,:) ;
Pd = lyap(A0',N);
F1 = inv(R(1,1))*B1'*Pd*X*C1' * inv(C1*X*C1');
F2 = inv(R(2,2))*B2'*Pd*X*C2' * inv(C2*X*C2');

Fd = [F1*C1;F2*C2];

x = 0; % variavel para fazer a contagem do loop de controle 

Fdi = Fd; % Fdi = Fd da iteração
while x < 10

   %Fdi = Fd; % Fdi = Fd da iteração
   A0i = A-B*Fdi;
   Xi  = lyap(A0i,X0);
   Ni = Q + Fdi(1,:)'*R(1,1)*Fdi(1,:) + Fdi(2,:)'*R(2,2)*Fdi(2,:) ;
   Pdn = lyap(A0i',Ni);
   F1i = R(1,1) \ B1'*Pdn*Xi*C1' /(C1*Xi*C1');
   F2i = R(2,2) \ B2'*Pdn*Xi*C2' /(C2*Xi*C2');
   Fdn = [F1i*C1;F2i*C2];
   delta = norm(Fdi-Fdn);
   Fdi = Fdn;
      if delta < 1e-5
       break
   end

   x = x+1;
   end
F;
Fd;
Fdn;

% Índice de desempenho
Jc = x0' * Pd * x0 %índice caso centralizado
Jd = x0' * Pdn * x0 %índice caso centralizado
end