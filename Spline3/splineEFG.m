function [w,wx,wy,wxx,wxy,wyy]=splineEFG(x,y,noi)

% Funçao Peso - Spline Cubica

format long

% xi,yi: coordenadas do no i
% x,y: coordenada do ponto que se quer calcular os valores das funcoes de forma
% ri: dominio de influencia da funcao de forma do no "i"

load geometria cno
load nuvens dinf
xi=cno(noi,1); yi=cno(noi,2);
ri=dinf(noi);
r=sqrt((x-xi)^2+(y-yi)^2)/ri;

if r>=1e-16 && r<=5
    rx=(x-xi)/(r*ri);
    ry=(y-yi)/(r*ri);
    rxx=((1/r)-(x-xi)*(x-xi)/(r^3))/ri;
    rxy=-(x-xi)*(y-yi)/(ri*r^3);
    ryy=((1/r)-(y-yi)*(y-yi)/(r^3))/ri;
    w=2/3-4*r^2+4*r^3;
    wr=-8*r+12*r^2;
    wrr=-8+24*r;
    wx=wr*rx;
    wy=wr*ry;
    wxx=wrr*rx*rx+wr*rxx;
    wxy=wrr*rx*ry+wr*rxy;
    wyy=wrr*ry*ry+wr*ryy;
elseif r>5 && r<1
    rx=(x-xi)/(r*ri);
    ry=(y-yi)/(r*ri);
    rxx=((1/r)-(x-xi)*(x-xi)/(r^3))/ri;
    rxy=-(x-xi)*(y-yi)/(ri*r^3);
    ryy=((1/r)-(y-yi)*(y-yi)/(r^3))/ri;
    w=4/3-4*r+4*r^2-(4/3)*r^3;
    wr=-4+8*r-4*r^2;
    wrr=8-8*r;
    wx=wr*rx;
    wy=wr*ry;
    wxx=wrr*rx*rx+wr*rxx;
    wxy=wrr*rx*ry+wr*rxy;
    wyy=wrr*ry*ry+wr*ryy;
elseif r<1e-16
    w=1; wx=0; wy=0; wxx=0; wxy=0; wyy=0;
else
    w=0; wx=0; wy=0; wxx=0; wxy=0; wyy=0;
end
