function [Ku,Fu]=penalidade_u(CCu)

% Calcula a contribuicao do termo de penalidade no vetor forca e matriz
%  - devido a restricao de deslocamento
format long

[ni,cint,wi]=GaussL;
% ni - numero de pontos de integracao
% cint=[s_1 ... s_ni] - coordenadas dos pontos de integracao (parametrizada)
% w=[w_1 w_2 ... w_ni] - sao os pesos relativos aos pontos de integracao

load geometria cel cno nno
load nuvens spif

[nu,i]=size(CCu); % Numero de faces/dof com restricoes de deslocamento
Ku=zeros(5*nno);
Fu=zeros(5*nno,1);
for l=1:nu
    el=CCu(l,1); face=CCu(l,2); dofu=CCu(l,4);
    [x_face]=cface(el,face);
    mpen=zeros(3); mpen(dofu,dofu)=1;
    Fpen=zeros(3,1); Fpen(dofu,1)=CCu(l,3);
    for j=1:ni
        s=cint(j);
        [N,Ns]=formaL2(s);
        x=N*x_face(:,1); y=N*x_face(:,2);
        linha=(el-1)*ni*3+(face-1)*3+j;
        n=spif(linha,(1+nno)); Jx=spif(linha,1:n);
        fi=zeros(n); [mfi,fi]=interpol_mfi(x,y,n,Jx);
        [jac]=jacobianoL2(Ns,x_face);
        Fp=zeros(5*n,1);
        Fp=mfi'*mpen*Fpen*jac*wi(j);
        Fa=zeros(5*nno,1); [Fa]=union(Fp,Jx,nno);
        Fu=Fu+Fa;
        Kpen=zeros(n,n); Kp=zeros(5*n,5*n);
        Kp=mfi'*mpen*mfi*jac*wi(j);
        Ka=zeros(5*nno); [Ka]=assemb(Kp,Jx,nno);
        Ku=Ku+Ka;
    end
end
