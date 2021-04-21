function [Fq]=Fsuperficie(q)

% Calcula o vetor de carregamento superficial
format long

load geometria cno cel nno
load nuvens spis
% ni - numero de pontos de integracao
% cint=[r1 r2 ...: s1 s2 ...] coordenadas dos pontos de integracao (parametrizada)
% w=[w1 w2 ...] sao os pesos relativos aos pontos de integracao
[ni,cint,wi]=GaussA;

Fq=zeros(5*nno,1);
[nq,i]=size(q); % Numero de elementos com tracoes prescritas na superficie
for j=1:nq
    el=q(j,1); x_el=CoordElem(el);
    for i=1:ni
        r=cint(1,i); s=cint(2,i);
        [N,Nr,Ns]=format3(r,s);
        [jac]=jacobiano(el,r,s);
        x=N*x_el(:,1); y=N*x_el(:,2);
        linha=(el-1)*ni+i;
        n=spis(linha,(1+nno)); Jx=spis(linha,1:n);
        mfi=zeros(3,5*n); [mfi,fi]=interpol_mfi(x,y,n,Jx);
        Fp=zeros(5*n,1); Fp=mfi'*q(j,2:4)'*jac*wi(i);
        Fa=zeros(5*nno,1); [Fa]=union(Fp,Jx,nno);
        Fq=Fq+Fa;
    end
end
