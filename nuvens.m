function [dinf,sno,spis,spif]=nuvens(a)

% Calcula snu, sno, spis, spif

% snu - suporte de cada nuvem (abrangencia)
% sno - suporte de cada no (n,Jx)
% spis - suporte de cada ponto de integracao da superficie  (n,Jx)
% spif - suporte de cada ponto de integracao das faces (arestas) (n,Jx)

load geometria nno cno cel nel

% Dominio de influencia de cada no (funcao peso)
snu=zeros(1,nno);
for i=1:nel
    no1=cel(i,1); x1=cno(no1,1); y1=cno(no1,2);
    no2=cel(i,2); x2=cno(no2,1); y2=cno(no2,2);
    no3=cel(i,3); x3=cno(no3,1); y3=cno(no3,2);
    d1=a*sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    d2=a*sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2));
    d3=a*sqrt((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3));
    if d1>snu(no1); snu(no1)=d1; end
    if d1>snu(no2); snu(no2)=d1; end
    if d2>snu(no2); snu(no2)=d2; end
    if d2>snu(no3); snu(no3)=d2; end
    if d3>snu(no1); snu(no1)=d3; end
    if d3>snu(no3); snu(no3)=d3; end
end
dinf=snu;
for i=1:nel
    no1=cel(i,1); no2=cel(i,2); no3=cel(i,3);
    if snu(no2)>snu(no1); dinf(no1)=snu(no2); end
    if snu(no3)>snu(no1); dinf(no1)=snu(no3); end
    if snu(no1)>snu(no2); dinf(no2)=snu(no1); end
    if snu(no3)>snu(no2); dinf(no2)=snu(no3); end
    if snu(no1)>snu(no3); dinf(no3)=snu(no1); end
    if snu(no2)>snu(no3); dinf(no3)=snu(no2); end
end

% Suporte de cada no - sno
% sno=[Jx ... n; ...]
sno=zeros(nno,(nno+1));
for j=1:nno
    x=cno(j,1); y=cno(j,2);
    n=0;
    for l=1:nno
        noi=l; xi=cno(noi,1); yi=cno(noi,2);
        di=sqrt((xi-x)^2+(yi-y)^2);
        if dinf(noi)>di
            n=n+1;
            Jx(n)=noi;
        end
    end
    sno(j,1:n)=Jx(1:n);
    sno(j,(nno+1))=n;
end
% Suporte de cada ponto de integracao da superficie (Gauss 12p) - spis
[ni,cint,wi]=GaussA;
spis=zeros(ni*nel,(nno+1));
for j=1:nel
    el=j; x_el=CoordElem(el);
    for i=1:ni
        r=cint(1,i); s=cint(2,i);
        [N,Nr,Ns]=format3(r,s);
        x=N*x_el(:,1); y=N*x_el(:,2);
        n=0;
        for l=1:nno
            noi=l; xi=cno(noi,1); yi=cno(noi,2);
            di=sqrt((xi-x)^2+(yi-y)^2);
            if dinf(noi)>di
                n=n+1;
                Jx(n)=noi;
            end
        end
        linha=(el-1)*ni+i;
        spis(linha,1:n)=Jx(1:n);
        spis(linha,(nno+1))=n;
    end
end
% Suporte de cada ponto de integracao da face (Gauss 3p) - spif
[ni,cint,wi]=GaussL;
spif=zeros(ni*3*nel,(nno+1));
for j=1:nel
    el=j;
    % Face 1
    x_face=cface(el,1);
    for i=1:ni
        s=cint(i);
        [N,Ns]=formaL2(s);
        x=N*x_face(:,1); y=N*x_face(:,2);
        n=0;
        for m=1:nno
            noi=m; xi=cno(noi,1); yi=cno(noi,2);
            di=sqrt((xi-x)^2+(yi-y)^2);
            if dinf(noi)>di
                n=n+1;
                Jx(n)=noi;
            end
        end
        linha=(el-1)*3*ni+i;
        spif(linha,1:n)=Jx(1:n);
        spif(linha,(1+nno))=n;
    end
    % Face 2
    x_face=cface(el,2);
    for i=1:ni
        s=cint(i); ws=wi(i);
        [N,Ns]=formaL2(s);
        x=N*x_face(:,1); y=N*x_face(:,2);
        n=0;
        for m=1:nno
            noi=m; xi=cno(noi,1); yi=cno(noi,2);
            di=sqrt((xi-x)^2+(yi-y)^2);
            if dinf(noi)>di
                n=n+1;
                Jx(n)=noi;
            end
        end
        linha=(el-1)*3*ni+ni+i;
        spif(linha,1:n)=Jx(1:n);
        spif(linha,(1+nno))=n;
    end
    % Face 3
    x_face=cface(el,3);
    for i=1:ni
        s=cint(i); ws=wi(i);
        [N,Ns]=formaL2(s);
        x=N*x_face(:,1); y=N*x_face(:,2);
        n=0;
        for m=1:nno
            noi=m; xi=cno(noi,1); yi=cno(noi,2);
            di=sqrt((xi-x)^2+(yi-y)^2);
            if dinf(noi)>di
                n=n+1;
                Jx(n)=noi;
            end
        end
        linha=(el-1)*3*ni+2*ni+i;
        spif(linha,1:n)=Jx(1:n);
        spif(linha,(1+nno))=n;
    end
end