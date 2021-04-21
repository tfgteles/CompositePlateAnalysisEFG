function [ni,cint,wi]=GaussA

% Integracao numerica de um triangulo
format long

% Fornece:
% numeros dos pontos de integracao
ni=12;
% Coordenadas parametricas dos pontos de integracao (cint)
a=0.063089014491502;
b=0.249286745170910;
c=0.310352451033785;
d=0.053145049844816;
% coordenadas x, linha 1 de cint
cint(1,:)=[a (1-2*a) a b (1-2*b) b c d (1-(c+d)) (1-(c+d)) c d];
% coordenadas y, linha 2 de cint
cint(2,:)=[a a (1-2*a) b b (1-2*b) d c c d (1-(c+d)) (1-(c+d))];

% Pesos
e=0.025422453185103;
f=0.058393137863189;
g=0.041425537809187;
wi=[e e e f f f g g g g g g];