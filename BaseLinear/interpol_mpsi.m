function [mpsi,fix,fiy]=interpol_mpsi(x,y,n,Jx)

% Funçao de Forma EFG - Polinomio Linear
format long

% xi,yi: coordenadas do no i
% x,y: coordenada do ponto que se quer calcular os valores das funcoes de forma

load geometria cno

A=zeros(6); Ax=zeros(6); Ay=zeros(6);
B=zeros(6,n); Bx=zeros(6,n); By=zeros(6,n);
for i=1:n
    noi=Jx(i); xi=cno(noi,1); yi=cno(noi,2);
    [w,wx,wy,wxx,wxy,wyy]=splineEFG(x,y,noi);
    pi2=[1; xi; yi];
    pixpi(:,1)=[1; xi; yi];
    pixpi(:,2)=[xi; xi^2; yi*xi];
    pixpi(:,3)=[yi; xi*yi; yi*yi];
    B(:,i)=w*pi2;
    Bx(:,i)=wx*pi2;
    By(:,i)=wy*pi2;
    Ai=w*pixpi;
    Aix=wx*pixpi;
    Aiy=wy*pixpi;
    A=A+Ai;
    Ax=Ax+Aix;
    Ay=Ay+Aiy;
end
Ainv=inv(A);
Ainvx=-Ainv*Ax*Ainv;
Ainvy=-Ainv*Ay*Ainv;
p=[1; x; y];
px=[0; 1; 0];
py=[0; 0; 1];
fi=p'*Ainv*B;
fix=px'*Ainv*B+p'*Ainvx*B+p'*Ainv*Bx;
fiy=py'*Ainv*B+p'*Ainvy*B+p'*Ainv*By;

mpsi=zeros(2,5*n);
for i=1:n
    mpsi(1,(i-1)*5+4)=fix(i);
    mpsi(2,(i-1)*5+5)=fiy(i);
end
