function [C]=constitutiva

% Fornece a matriz constitutiva a partir das propriedades materiais

format long
% E1 - Modulo de Young (MPa) na direcao 1 (fibra)
% E2 - Modulo de Young (MPa) na direcao 2 (matriz)
% G12 - Modulo de cisalhamento (MPa) no plano 1-2
% G13 - Modulo de cisalhamento (MPa) no plano 1-3
% G23 - Modulo de cisalhamento (MPa) no plano 2-3
% v12 - Coeficiente de Poisson referente a deformacao na direcao 2 causada por uma solicitacao na direcao 1
% k=5/6 - Fator de correcao da tensao cisalhante
% h - Espessura da placa
% emp=[o1 o2 o3 ...] - Sequencia de empilhamento das laminas (graus)
% Material composto laminado e estado plano de tensao

load geometria h
load CteMat E1 E2 G12 G13 G23 v12 k emp
SLx=[1/E1 -v12/E1 0; -v12/E1 1/E2 0; 0 0 1/G12];
QLx=inv(SLx);
SLz=[1/G23 0; 0 1/G13];
QLz=inv(SLz);
A=zeros(3); B=zeros(3);
D=zeros(3); F=zeros(2);
[i,n]=size(emp);
hk=h/n; % espessura de cada lamina
for i=1:n
    teta=emp(i)*pi/180;
    s=sin(teta); c=cos(teta);
    Tx=[c*c s*s -2*s*c; s*s c*c 2*s*c; s*c -s*c (c*c-s*s)];
    Tz=[c s; -s c];
    Qx=Tx*QLx*Tx';
    Qz=Tz*QLz*Tz';
    zki=(i-1)*hk-h/2; zks=zki+hk; % coordenada z inferior e superior da lamina
    Ai=zeros(3); Ai=(zks-zki)*Qx; A=A+Ai;
    Bi=zeros(3); Bi=((zks^2-zki^2)/2)*Qx; B=B+Bi;
    Di=zeros(3); Di=((zks^3-zki^3)/3)*Qx; D=D+Di;
    Fi=zeros(2); Fi=(zks-zki)*Qz;
    F=F+Fi;
end
C(1:3,1:3)=A;
C(1:3,4:6)=B;
C(4:6,1:3)=B;
C(4:6,4:6)=D;
C(7:8,7:8)=k*F;