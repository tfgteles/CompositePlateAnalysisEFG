function FN=FbordaN(Nt)

% Calcula o vetor de carregamento das forcas que atuam na borda

format long

[ni,cint,wi]=GaussL;
% ni - numero de pontos de integracao
% cint=[s_1 ... s_ni] - coordenadas dos pontos de integracao (parametrizada)
% w=[w_1 w_2 ... w_ni] - sao os pesos relativos aos pontos de integracao

load geometria nno
load nuvens spif

FN=zeros(5*nno,1);
[nN,i]=size(Nt); % Numero de faces com tracoes prescritas (borda)
for i=1:nN
    el=Nt(i,1); face=Nt(i,2);
    [x_face]=cface(el,face);
    for j=1:ni
        s=cint(j);
        linha=(el-1)*ni*3+(face-1)*3+j;
        n=spif(linha,(1+nno)); Jx=spif(linha,1:n);
        [N,Ns]=formaL2(s);
        x=N*x_face(:,1); y=N*x_face(:,2);
        mfi=zeros(3,5*n); [mfi,fi]=interpol_mfi(x,y,n,Jx);
        [jac]=jacobianoL2(Ns,x_face);
        Fp=mfi'*Nt(i,3:5)'*jac*wi(j);
        Fa=zeros(5*nno,1); [Fa]=union(Fp,Jx,nno);
        FN=FN+Fa;
    end
end
