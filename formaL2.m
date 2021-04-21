function [N,Ns]=formaL2(s)

% Calcula o valor das funcoes isoparametricas no ponto parametrizado (s)
format long

% Elemento unidimensional linear (4 nos por elemento)
N(1)=(1-s)/2;
N(2)=(1+s)/2;

Ns(1)=-1/2;
Ns(2)=1/2;
