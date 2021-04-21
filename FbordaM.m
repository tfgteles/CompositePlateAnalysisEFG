function FM=FbordaM(Mt)

% Calcula o vetor de carregamento dos momentos que atuam na borda
format long

[ni,cint,wi]=GaussL;
% ni - numero de pontos de integracao
% cint=[s_1 ... s_ni] - coordenadas dos pontos de integracao (parametrizada)
% w=[w_1 w_2 ... w_ni] - sao os pesos relativos aos pontos de integracao

load geometria nno
load nuvens spif
FM=zeros(5*nno,1);
[nM,i]=size(Mt); % Numero de faces/dof com momentos prescritas (borda)
for i=1:nM
    el=Mt(i,1); face=Mt(i,2);
    [x_face]=cface(el,face);
    for j=1:ni
        s=cint(j);
        linha=(el-1)*ni*3+(face-1)*3+j;
        n=spif(linha,(1+nno)); Jx=spif(linha,1:n);
        [N,Ns]=formaL2(s);
        x=N*x_face(:,1); y=N*x_face(:,2);
        mpsi=zeros(2,5*n); [mpsi,fix,fiy]=interpol_mpsi(x,y,n,Jx);
        [jac]=jacobianoL2(Ns,x_face);
        Fp=mpsi'*Mt(i,3:4)'*jac*wi(j);
        Fa=zeros(5*nno,1); [Fa]=union(Fp,Jx,nno);
        FM=FM+Fa;
    end
end
