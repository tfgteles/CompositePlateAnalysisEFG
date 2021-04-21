function [cno,cel,nno,nel]=geometria(ladoA,ladoB,nA,nB)

% Fornece as coordenadas nodais e a conectividade dos elementos

% Corpo retangular com [(nA+1)*(nB+1)] nos e (2*nA*nB) elementos t3
format long
dA=ladoA/nA;
dB=ladoB/nB;
nno=0;
for i=1:(nB+1)
    for j=1:(nA+1)
        x=(j-1)*dA;
        y=(i-1)*dB;
        nno=nno+1;
        cno(nno,:)=[x y];
    end
end
% Elemento t4
nel=0;
for i=1:nB
    for j=1:nA
        no1=j+(i-1)*(nA+1);
        no2=no1+1; no3=no1+nA+2; no4=no1+nA+1;
        nel=nel+1;
        cel(nel,:)=[no1 no2 no3];
        nel=nel+1;
        cel(nel,:)=[no3 no4 no1];
    end
end