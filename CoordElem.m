function x_el=CoordElem(el)

% Fornece as coordenadas x,y dos nos do elemento

% x_el=[x1 y1; x2 y2; ...]

format long

load geometria cel cno
no1=cel(el,1); x1=cno(no1,1); y1=cno(no1,2);
no2=cel(el,2); x2=cno(no2,1); y2=cno(no2,2);
no3=cel(el,3); x3=cno(no3,1); y3=cno(no3,2);
x_el=[x1 y1; x2 y2; x3 y3];