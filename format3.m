function [N,Nr,Ns]=format3(r,s)

% Calcula o valor das funcoes isoparametricas no ponto parametrizado (r,s)
format long

% Elemento triangular linear (3 nos por elemento)
N(1)=1-r-s;
N(2)=r;
N(3)=s;

Nr(1)=-1;
Nr(2)=1;
Nr(3)=0;

Ns(1)=-1;
Ns(2)=0;
Ns(3)=1;
