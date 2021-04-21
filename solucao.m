function [u,w,Sxx,Syy,Sxy,Syz,Sxz]=solucao(x,y,q0,d)

% Fornece a tensao e o deslocamento de um ponto (x,y,z) qualquer

format long
load geometria cno nno h ladoA ladoB
load CteMat E1 E2 G12 G13 G23 v12 emp
[n,Jx]=suporte(x,y);
[fi,fix,fiy,fixx,fixy,fiyy]=formaEFG(x,y,n,Jx);

% Parametros nodais do suporte de (x,y)
pno=zeros(1,5*n);
for j=1:n
    noj=Jx(j);
    pno((j-1)*5+1)=d(1+(noj-1)*5,1); % parametro nodal u
    pno((j-1)*5+2)=d(2+(noj-1)*5,1); % parametro nodal v
    pno((j-1)*5+3)=d(3+(noj-1)*5,1); % parametro nodal w
    pno((j-1)*5+4)=d(4+(noj-1)*5,1); % parametro nodal -Oy
    pno((j-1)*5+5)=d(5+(noj-1)*5,1); % parametro nodal Ox
end

% Deslocamento u0(x,y)=[u; v; w]
mfi=zeros(5,5*n);
for i=1:n
    mfi(1,(i-1)*5+1)=fi(i);
    mfi(2,(i-1)*5+2)=fi(i);
    mfi(3,(i-1)*5+3)=fi(i);
    mfi(4,(i-1)*5+4)=fix(i);
    mfi(5,(i-1)*5+5)=fiy(i);
end
u=mfi*pno';

% Deslocamento em z normalizado
w=u(3,1)*(E2*h^3)*100/(q0*ladoA^4);

% B = matriz de deformacao no ponto de interesse
B=zeros(8,n*5);
for l=1:n
    B(1,((l-1)*5+1):(l*5))=[fix(l) 0 0 0 0];
    B(2,((l-1)*5+1):(l*5))=[0 fiy(l) 0 0 0];
    B(3,((l-1)*5+1):(l*5))=[fiy(l) fix(l) 0 0 0];
    B(4,((l-1)*5+1):(l*5))=[0 0 0 fixx(l) 0];
    B(5,((l-1)*5+1):(l*5))=[0 0 0 0 fiyy(l)];
    B(6,((l-1)*5+1):(l*5))=[0 0 0 fixy(l) fixy(l)];
    B(7,((l-1)*5+1):(l*5))=[0 0 fiy(l) 0 fiy(l)];
    B(8,((l-1)*5+1):(l*5))=[0 0 fix(l) fix(l) 0];
end

% Deformacoes generalizadas
e_gen=B*pno';

% Fatores de normalizacao das tensoes
f1=(h^2)/(q0*ladoA^2);
f2=h/(q0*ladoA);

[i,nk]=size(emp); % numero de laminas
hk=h/nk; % espessura de cada lamina
SLx=[1/E1 -v12/E1 0; -v12/E1 1/E2 0; 0 0 1/G12];
QLx=inv(SLx);
SLz=[1/G23 0; 0 1/G13];
QLz=inv(SLz);
for i=1:nk
    zki=(i-1)*hk-h/2;
    zks=zki+hk;
    exxi=e_gen(1,1)+zki*e_gen(4,1);
    eyyi=e_gen(2,1)+zki*e_gen(5,1);
    gxyi=e_gen(3,1)+zki*e_gen(6,1);
    exxs=e_gen(1,1)+zks*e_gen(4,1);
    eyys=e_gen(2,1)+zks*e_gen(5,1);
    gxys=e_gen(3,1)+zks*e_gen(6,1);
    gyz=e_gen(7,1);
    gxz=e_gen(8,1);
    teta=emp(i)*pi/180;
    s=sin(teta); c=cos(teta);
    Qx=zeros(3); Qz=zeros(2);
    Tx=[c*c s*s -2*s*c; s*s c*c 2*s*c; s*c -s*c (c*c-s*s)];
    Tz=[c s; -s c];
    Qx=Tx*QLx*Tx';
    Qz=Tz*QLz*Tz';
    Sxi=Qx*[exxi; eyyi; gxyi]*f1;
    Sxs=Qx*[exxs; eyys; gxys]*f1;
    Sz=Qz*[gyz; gxz]*f2;
    Sxx(i,:)=[Sxi(1,1) Sxs(1,1)];
    Syy(i,:)=[Sxi(2,1) Sxs(2,1)];
    Sxy(i,:)=[Sxi(3,1) Sxs(3,1)];
    Syz(i,1)=Sz(1,1);
    Sxz(i,1)=Sz(2,1);
    stress_xx=Sxx(i,2)/f1
    stress_yy=Syy(i,2)/f1
end
