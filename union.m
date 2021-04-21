function [Fa]=union(Fp,Jx,nno)

% Monta o vetor forca elementar no global

format long

[i,n]=size(Jx);
Fa=zeros(5*nno,1);
for i=1:n
    no=Jx(i);
    u=Fp(((i-1)*5+1),1);
    v=Fp(((i-1)*5+2),1);
    w=Fp(((i-1)*5+3),1);
    alfa=Fp(((i-1)*5+4),1);
    beta=Fp(((i-1)*5+5),1);
    dof1=(no-1)*5+1;
    dof2=(no-1)*5+2;
    dof3=(no-1)*5+3;
    dof4=(no-1)*5+4;
    dof5=(no-1)*5+5;
    Fa(dof1,1)=u;
    Fa(dof2,1)=v;
    Fa(dof3,1)=w;
    Fa(dof4,1)=alfa;
    Fa(dof5,1)=beta;
end