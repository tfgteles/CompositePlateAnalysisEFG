function [n,Jx]=suporte(x,y)

% Nos cujos suportes contem o ponto (x,y)

load geometria cno nno
load nuvens dinf
n=0;
for i=1:nno
    noi=i; xi=cno(noi,1); yi=cno(noi,2);
    di=sqrt((xi-x)^2+(yi-y)^2);
    if dinf(noi)>di
        n=n+1;
        Jx(n)=noi;
    end
end