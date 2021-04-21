function [x_face]=cface(el,face);

% Dado o elemento e a face calcula-se as coordenadas dos nos da face

load geometria cel cno
if face==1
    no1=cel(el,1);
    no2=cel(el,2);
    x1=cno(no1,1); y1=cno(no1,2);
    x2=cno(no2,1); y2=cno(no2,2);
    x_face=[x1 y1; x2 y2];
end
if face==2
    no1=cel(el,2);
    no2=cel(el,3);
    x1=cno(no1,1); y1=cno(no1,2);
    x2=cno(no2,1); y2=cno(no2,2);
    x_face=[x1 y1; x2 y2];
end
if face==3
    no1=cel(el,3);
    no2=cel(el,1);
    x1=cno(no1,1); y1=cno(no1,2);
    x2=cno(no2,1); y2=cno(no2,2);
    x_face=[x1 y1; x2 y2];
end
