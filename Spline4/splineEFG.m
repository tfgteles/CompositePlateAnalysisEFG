function [w,wx,wy,wxx,wxy,wyy]=splineEFG(x,y,noi)

% Fun?ao Peso - Spline Quartica

format long

% xi,yi: coordenadas do no i
% x,y: coordenada do ponto que se quer calcular os valores das funcoes de forma
% ri: dominio de influencia da funcao de forma do no "i"

load geometria cno
load nuvens dinf
xi=cno(noi,1); yi=cno(noi,2);
ri=dinf(noi);
r=sqrt((x-xi)^2+(y-yi)^2)/ri;


if r<1e-16
    w=1; wx=0; wy=0; wxx=0; wxy=0; wyy=0;
elseif r>1
    w=0; wx=0; wy=0; wxx=0; wxy=0; wyy=0;
else
    rx=(x-xi)/(r*ri);
    ry=(y-yi)/(r*ri);
    rxx=((1/r)-(x-xi)*(x-xi)/(r^3))/ri;
    rxy=-(x-xi)*(y-yi)/(ri*r^3);
    ryy=((1/r)-(y-yi)*(y-yi)/(r^3))/ri;
    w=1-6*r^2+8*r^3-3*r^4;
    wr=-12*r+24*r^2-12*r^3;
    wrr=-12+48*r-36*r^2;
    wx=wr*rx;
    wy=wr*ry;
    wxx=wrr*rx*rx+wr*rxx;
    wxy=wrr*rx*ry+wr*rxy;
    wyy=wrr*ry*ry+wr*ryy;
end
