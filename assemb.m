function [Ka]=assemb(Kp,Jx,nno);

% Montagem da matriz de rigidez elementar na matriz global
format long

% Kp - matriz de rigidez do ponto de integracao
% Jx - sequencia de nos que compoe o suporte do ponto de integracao
% ndof - mapeamento dos graus de liberdade do sistema por no


[i,n]=size(Jx);
Ka=zeros(5*nno,5*nno);
for i=1:n
    no=Jx(i);
    for j=1:5
        dof((i-1)*5+j)=(no-1)*5+j;
    end
end
for i=1:n*5
    for j=1:n*5
        Ka(dof(i),dof(j))=Kp(i,j);
    end
end