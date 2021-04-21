function [jac]=jacobiano(el,r,s)

% Calcula o determinante do jacobiano e sua funcao inversa no ponto (r,s)
format long
[N,Nr,Ns]=format3(r,s);
x_el=CoordElem(el);
[i,n]=size(Nr);
for i=1:n
    Nrs(1,i)=Nr(i);
    Nrs(2,i)=Ns(i);
end
jacm=Nrs*x_el; % matriz jacobiana
jac=jacm(1,1)*jacm(2,2)-jacm(1,2)*jacm(2,1); % determinante da matriz do jacobiano
