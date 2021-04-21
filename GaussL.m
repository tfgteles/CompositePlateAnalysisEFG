function [ni,cint,wi]=GaussL

% Integracao numerica linear
format long

% Fornece:
% numeros dos pontos de integracao
ni=3;
% Coordenadas parametricas dos pontos de integracao
cint=[-sqrt(0.6) 0 sqrt(0.6)];
% Pesos
wi=[5/9 8/9 5/9];