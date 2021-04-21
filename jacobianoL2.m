function [jac]=jacobianoL2(Ns,x_face)

% Calcula o determinante do jacobiano para elemento unidimensional linear
format long

xi=x_face(:,1);
yi=x_face(:,2);
dxds=Ns*xi;
dyds=Ns*yi;

jac=sqrt(dxds*dxds+dyds*dyds);