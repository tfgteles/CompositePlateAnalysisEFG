function [E1,E2,G12,G13,G23,v12,k]=prop_mat(tpMat)

% Fornece as constantes materiais a partir do tipo pre-definido

% E1:   Modulo de Young (Pa) na direcao 1 (fibra)
% E2:   Modulo de Young (Pa) na direcao 2 (matriz)
% G12:  Modulo de cisalhamento (Pa) no plano 1-2
% G13:  Modulo de cisalhamento (Pa) no plano 1-3
% G23:  Modulo de cisalhamento (Pa) no plano 2-3
% v12:  Coeficiente de Poisson referente a deformacao na direcao 2 causada por uma solicitacao na direcao 1
% k:    fator de correcao da tensao cisalhante
% emp:  Sequencia de empilhamento das laminas (graus)

format long

if tpMat==1
    E1=210e9;
    E2=210e9;
    v12=0.3;
    G12=E1/(2*(1+v12));
    G13=G12;
    G23=G12;
    k=5/6;
end
if tpMat==2
    E2=7e9;
    E1=25*E2;
    G12=E2/2;
    G13=G12;
    G23=E2/5;
    v12=0.25;
    k=5/6;
end
if tpMat==3
    E1=2e9;
    E2=2e9;
    v12=0.3;
    G12=E1/(2*(1+v12));
    G13=G12;
    G23=G12;
    k=5/6;
end
if tpMat==4
    E1=200e9;
    E2=5e9;
    G12=3e9;
    G13=3e9;
    G23=2.5e9;
    v12=0.25;
    k=5/6;
end
if tpMat==5
    E1=45e9;
    E2=12e9;
    G12=4.5e9;
    G13=1e9;
    G23=1e9;
    v12=0.3;
    k=1;
end
