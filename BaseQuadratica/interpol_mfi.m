function [mfi,fi]=interpol_mfi(x,y,n,Jx)

% Funçao de Forma EFG - Polinomio Quadratico
format long

% xi,yi: coordenadas do no i
% x,y: coordenada do ponto que se quer calcular os valores das funcoes de forma


load geometria cno

A=zeros(6);
B=zeros(6,n);
for i=1:n
    noi=Jx(i); xi=cno(noi,1); yi=cno(noi,2);
    [w,wx,wy,wxx,wxy,wyy]=splineEFG(x,y,noi);
    pi2=[1; xi; yi; xi^2; xi*yi; yi^2];
    pixpi(:,1)=[1; xi; yi; xi^2; xi*yi; yi^2];
    pixpi(:,2)=[xi; xi^2; yi*xi; xi^3; xi*xi*yi; xi*yi^2];
    pixpi(:,3)=[yi; xi*yi; yi*yi; yi*xi^2; xi*yi^2; yi^3];
    pixpi(:,4)=[xi^2; xi^3; yi*xi^2; xi^4; yi*xi^3; (xi^2)*(yi^2)];
    pixpi(:,5)=[xi*yi; yi*xi^2; xi*yi^2; yi*xi^3; (xi^2)*(yi^2); xi*yi^3];
    pixpi(:,6)=[yi^2; xi*yi^2; yi^3; (xi^2)*(yi^2); xi*yi^3; yi^4];
    B(:,i)=w*pi2;
    Ai=w*pixpi;
    A=A+Ai;
end
Ainv=inv(A);
p=[1; x; y; x^2; x*y; y^2];
fi=p'*Ainv*B;

mfi=zeros(3,5*n);
for i=1:n
    mfi(1,(i-1)*5+1)=fi(i);
    mfi(2,(i-1)*5+2)=fi(i);
    mfi(3,(i-1)*5+3)=fi(i);
end