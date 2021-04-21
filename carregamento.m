function [Nt,Mt,q]=carregamento(tpLoad,q0,V0)

% Fornece carregamentos pre-definidos

format long

% tpLoad=1: superficie uniformemente distribuida q0, borda nula
% tpLoad=2: superficie senoidalmente distrubuida q0, borda nula
% tpLoad=3: carregamento V0 direcao x na borda x=ladoA, superficie nula

% Carregamento no contorno dos elementos (borda)
% Nt=[el face Nx Ny Nz;...]
% Mt=[el face Mx My]
% Face 1: entre os nos locais 1 e 2
% Face 2: entre os nos locais 2 e 3
% Face 3: entre os nos locais 3 e 1

% Carregamento na superficie da placa - pt
% pt=[el Nx Ny Nz;...]


if tpLoad==1
    Nt=[1 1 0 0 0]; % carregamento nulo na borda
    Mt=[1 1 0 0];
    load geometria nel
    for i=1:nel
        q(i,:)=[i 0 0 q0];
    end
end
if tpLoad==2
    Nt=[1 1 0 0 0]; % carregamento nulo na borda
    Mt=[1 1 0 0];
    load geometria nel cel cno ladoA ladoB
    for i=1:nel
        no1=cel(i,1); x1=cno(no1,1); y1=cno(no1,2);
        no2=cel(i,2); x2=cno(no2,1); y2=cno(no2,2);
        no3=cel(i,3); x3=cno(no3,1); y3=cno(no3,2);
        x=(x1+x2+x3)/3; y=(y1+y2+y3)/3;
        V=q0*sin(pi*x/ladoA)*sin(pi*y/ladoB);
        q(i,:)=[i 0 0 V];
    end
end
if tpLoad==3
    Mt=[1 1 0 0];
    load geometria nA nB
    for i=1:nB
        el=i*(2*nA)-1;
        Nt(i,:)=[el 2 V0 0 0];
    q(1,:)=[1 0 0 0];
    end
end
