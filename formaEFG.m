function [fi,fix,fiy,fixx,fixy,fiyy]=formaEFG(x,y,n,Jx)

% Funçao de Forma EFG - Polinomio Quadratico
format long

% xi,yi: coordenadas do no i
% x,y: coordenada do ponto que se quer calcular os valores das funcoes de forma

load geometria cno

A=zeros(6); Ax=zeros(6); Ay=zeros(6); Axx=zeros(6); Axy=zeros(6); Ayy=zeros(6);
B=zeros(6,n); Bx=zeros(6,n); By=zeros(6,n); Bxx=zeros(6,n); Bxy=zeros(6,n); Byy=zeros(6,n);
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
    Bx(:,i)=wx*pi2;
    By(:,i)=wy*pi2;
    Bxx(:,i)=wxx*pi2;
    Bxy(:,i)=wxy*pi2;
    Byy(:,i)=wyy*pi2;
    Ai=w*pixpi;
    Aix=wx*pixpi;
    Aiy=wy*pixpi;
    Aixx=wxx*pixpi;
    Aixy=wxy*pixpi;
    Aiyy=wyy*pixpi;
    A=A+Ai;
    Ax=Ax+Aix;
    Ay=Ay+Aiy;
    Axx=Axx+Aixx;
    Axy=Axy+Aixy;
    Ayy=Ayy+Aiyy;
end
Ainv=inv(A);
Ainvx=-Ainv*Ax*Ainv;
Ainvy=-Ainv*Ay*Ainv;
Ainvxx=-Ainvx*Ax*Ainv-Ainv*Ax*Ainvx-Ainv*Axx*Ainv;
Ainvxy=-Ainvx*Ay*Ainv-Ainv*Ay*Ainvx-Ainv*Axy*Ainv;
Ainvyy=-Ainvy*Ay*Ainv-Ainv*Ay*Ainvy-Ainv*Ayy*Ainv;
p=[1; x; y; x^2; x*y; y^2];
px=[0; 1; 0; 2*x; y; 0];
py=[0; 0; 1; 0; x; 2*y];
pxx=[0; 0; 0; 2; 0; 0];
pxy=[0; 0; 0; 0; 1; 0];
pyy=[0; 0; 0; 0; 0; 2];
fi=p'*Ainv*B;
fix=px'*Ainv*B+p'*Ainvx*B+p'*Ainv*Bx;
fiy=py'*Ainv*B+p'*Ainvy*B+p'*Ainv*By;
fixx=pxx'*Ainv*B+px'*Ainvx*B+px'*Ainv*Bx+px'*Ainvx*B+p'*Ainvxx*B+p'*Ainvx*Bx+px'*Ainv*Bx+p'*Ainvx*Bx+p'*Ainv*Bxx;
fixy=pxy'*Ainv*B+px'*Ainvy*B+px'*Ainv*By+py'*Ainvx*B+p'*Ainvxy*B+p'*Ainvx*By+py'*Ainv*Bx+p'*Ainvy*Bx+p'*Ainv*Bxy;
fiyy=pyy'*Ainv*B+py'*Ainvy*B+py'*Ainv*By+py'*Ainvy*B+p'*Ainvyy*B+p'*Ainvy*By+py'*Ainv*By+p'*Ainvy*By+p'*Ainv*Byy;
