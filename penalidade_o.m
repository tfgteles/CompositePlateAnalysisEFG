function [Ko,Fo]=penalidade_o(CCo)

% Calcula a contribuicao do termo de penalidade no vetor forca e matriz
%  - devido a restricao de deslocamento
format long

[ni,cint,wi]=GaussL;
% ni - numero de pontos de integracao
% cint=[s_1 ... s_ni] - coordenadas dos pontos de integracao (parametrizada)
% w=[w_1 w_2 ... w_ni] - sao os pesos relativos aos pontos de integracao

load geometria cel cno nno
load nuvens spif

[no,i]=size(CCo); % Numero de faces/dof com restricoes de deslocamento
Ko=zeros(5*nno);
Fo=zeros(5*nno,1);
if CCo(1,1)>0
    for l=1:no
        el=CCo(l,1); face=CCo(l,2); dofo=CCo(l,4);
        [x_face]=cface(el,face);
        mpen=zeros(2); mpen(dofo,dofo)=1;
        Fpen=zeros(2,1); Fpen(dofo,1)=CCo(l,3);
        for j=1:ni
            s=cint(j);
            [N,Ns]=formaL2(s);
            x=N*x_face(:,1); y=N*x_face(:,2);
            linha=(el-1)*ni*3+(face-1)*3+j;
            n=spif(linha,(1+nno)); Jx=spif(linha,1:n);
            fix=zeros(n); fiy=zeros(n); [mpsi,fix,fiy]=interpol_mpsi(x,y,n,Jx);
            [jac]=jacobianoL2(Ns,x_face);
            Fp=zeros(n*5,1);
            Fp=mpsi'*mpen*Fpen*jac*wi(j);
            Fa=zeros(5*nno,1); [Fa]=union(Fp,Jx,nno);
            Fo=Fo+Fa;
            Kp=zeros(5*n);
            Kp=mpsi'*mpen*mpsi*jac*wi(j);
            Ka=zeros(5*nno,5*nno); [Ka]=assemb(Kp,Jx,nno);
            Ko=Ko+Ka;
        end
    end
end