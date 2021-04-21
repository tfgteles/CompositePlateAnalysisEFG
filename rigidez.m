function [K]=rigidez(C)

% Calcula a matriz de rigidez
format long

load geometria cel cno nno nel
load nuvens spis

% ni - numero de pontos de integracao
% cint=[r1 r2 ...: s1 s2 ...] coordenadas dos pontos de integracao (parametrizada)
% w=[w1 w2 ...] sao os pesos relativos aos pontos de integracao
[ni,cint,wi]=GaussA;

% Matriz de Rigidez
K=zeros(5*nno,5*nno);
for i=1:nel
    el=i;
    Kel=zeros(5*nno,5*nno);
    x_el=CoordElem(el);
    for j=1:ni
        r=cint(1,j); s=cint(2,j);
        [N,Nr,Ns]=format3(r,s); % funcoes de forma para elemento t3
        x=N*x_el(:,1); y=N*x_el(:,2);
        linha=(el-1)*ni+j;
        n=spis(linha,(1+nno)); Jx=spis(linha,1:n);
        [fi,fix,fiy,fixx,fixy,fiyy]=formaEFG(x,y,n,Jx);
        [jac]=jacobiano(el,r,s);
        % B = matriz de deformacao no ponto de integracao
        B=zeros(8,n*5);
        for l=1:n
            B(1,((l-1)*5+1):(l*5))=[fix(l) 0 0 0 0];
            B(2,((l-1)*5+1):(l*5))=[0 fiy(l) 0 0 0];
            B(3,((l-1)*5+1):(l*5))=[fiy(l) fix(l) 0 0 0];
            B(4,((l-1)*5+1):(l*5))=[0 0 0 fixx(l) 0];
            B(5,((l-1)*5+1):(l*5))=[0 0 0 0 fiyy(l)];
            B(6,((l-1)*5+1):(l*5))=[0 0 0 fixy(l) fixy(l)];
            B(7,((l-1)*5+1):(l*5))=[0 0 fiy(l) 0 fiy(l)];
            B(8,((l-1)*5+1):(l*5))=[0 0 fix(l) fix(l) 0];
        end
        Kp=zeros(5*n,5*n); Kp=B'*C*B*jac*wi(j);
        Ka=zeros(5*nno,5*nno); [Ka]=assemb(Kp,Jx,nno);
        K=K+Ka;
    end
end