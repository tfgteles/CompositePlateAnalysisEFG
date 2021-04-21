function [u0]=deslocamento(x,y,d)

% Calcula o deslocamento q do ponto (x,y) da superficie neutra

format long

load geometria cno nno
[n,Jx]=suporte(x,y);
[mfi,fi]=interpol_mfi(x,y,n,Jx);
pno=zeros(1,5*n);
for j=1:n
    noj=Jx(j);
    pno((j-1)*5+1)=d(1+(noj-1)*5,1); % parametro nodal u
    pno((j-1)*5+2)=d(2+(noj-1)*5,1); % parametro nodal v
    pno((j-1)*5+3)=d(3+(noj-1)*5,1); % parametro nodal w
end
u0=mfi*pno'; % u0=[u; v; w]

