% Element Free Galerkin Method
% Placa de Mindlin

clear all
format long

% -------------------------------------------
% Dados de Entrada
% -------------------------------------------

% --------------
% Geometria da placa

ladoA=1; % comprimento do lado A da placa (m)
ladoB=1; % comprimento do lado B da placa (m)
h=0.1; % espessura da placa (m)
nA=9; % numero de elementos (divisoes) do lado A da placa
nB=9; % numero de elementos (divisoes) do lado B da placa
malha=1; % malha 1 se regular, 0 se distorcida (5 x 5)
% --------------

% --------------
% Propriedades do material composto laminado
tpMat=1; emp=[0];
% emp:  Sequencia de empilhamento das laminas (graus), ex.: emp=[0 90 90 0]
% tpMat=1:  Material isotropico, E=210e9, v=0.3, k=5/6.
% tpMat=2:  Material composto, E1==25*7, E2=E1/25, G12=E2/2,
%           G13=G12, G23=E2/5, v12=0.25, k=5/6.
% tpMat=3:  Material isotropico, E=200e9, v=0.3, k=5/6.
% tpMat=4:  Material composto, E1==200e9, E2=5e9, G12=3e9,
%           G13=3e9, G23=2.5e9, v12=0.25, k=5/6.
% tpMat=5:  Material composto, teste.
% --------------

% --------------
% Carregamento
% tpload: carregamentos pre-definidos
q0=-10; % (N/m2) em z
V0=10; % (N/m) em z
tpLoad=3;
% tpLoad=1: superficie uniformemente distribuida q0, borda nula
% tpLoad=2: superficie senoidalmente distrubuida q0, borda nula
% tpLoad=3: carregamento V0 direcao x na borda x=ladoA, superficie nula
% --------------

% --------------
% Condicoes de Contorno Essenciais - Deslocamentos Prescritos
% tpCC: tipo de condicao de contorno predefinida
tpCC=6;
% tpCC=1: SS-1 (Reddy), apoio simples com restricao lateral
% tpCC=2: SS-2 (Reddy), apoio simples
% tpCC=3: placa suportada em y=0 e y=ladoB conforme SS-1
% tpCC=4: placa suportada em x=0 e x=ladoA conforme SS-2
% tpCC=5: placa suportada em x=0 (SS-1) e y=0 (SS-2)
% tpCC=6: placa engastada em x=0 (1/4 da placa)
% tpCC=7: placa simplesmente apoiada (1/4 da placa)
% tpCC=8: placa simplesmente apoiada (toda a placa)
% tpCC=9: placa totalmente engastada (toda a placa)
% --------------

% --------------
% Coeficiente de penalidade
%ep=10*E1*h;
% --------------

% --------------
% Abrangencia do suporte
a=1.5;
% --------------

% -------------------------------------------
% Pre-Processamento
% -------------------------------------------

% Geometria
% Coordenadas dos nos (cno) e conectividade dos elementos (cel)
% Numero de nos (nno) e de elementos (nel) do sistema
[cno,cel,nno,nel]=geometria(ladoA,ladoB,nA,nB);
if malha==1
    disp('malha regular');
elseif nA==4 && nB==4
    cno=PatchTeste(ladoA,ladoB);
    disp('malha distorcida');
end
save geometria ladoA ladoB h nA nB cno cel nno nel

% Carregamento
[Nt,Mt,q]=carregamento(tpLoad,q0,V0);

% Condicoes de Contorno Essenciais - Deslocamentos Prescritos
[CCu,CCo]=CondContorno(tpCC);

% Matriz constitutiva do laminado
[E1,E2,G12,G13,G23,v12,k]=prop_mat(tpMat);
save CteMat E1 E2 G12 G13 G23 v12 k emp
[C]=constitutiva;

% Nuvens: suporte dos nos, pontos de integracao e abrangencia das funcoes
[dinf,sno,spis,spif]=nuvens(a);
save nuvens dinf sno spis spif


% -------------------------------------------
% Solver
% -------------------------------------------

% Matriz de Rigidez
K=zeros(5*nno);
[K]=rigidez(C);

% Vetor de carregamento
FN=zeros(5*nno,1);
FN=FbordaN(Nt);
FM=zeros(5*nno,1);
FM=FbordaM(Mt);
Fq=zeros(5*nno,1);
Fq=Fsuperficie(q);

% Condicoes de Contorno Essenciais
Fu=zeros(5*nno,1); Ku=zeros(5*nno);
[Ku,Fu]=penalidade_u(CCu);
Fo=zeros(5*nno,1); Ko=zeros(5*nno);
[Ko,Fo]=penalidade_o(CCo);

% Solucao Nodal (parametros nodais)
ep=0.1*E1*h;
Kt=K+ep*1*Ku+ep*Ko;
Ft=FN+FM+Fq+ep*1*Fu+ep*Fo;
d=Kt\Ft;

% -------------------------------------------
% Pos-Processamento
% -------------------------------------------

% Grafico do deslocamento vertical (em z)
for i=1:(nB+1)
    for j=1:(nA+1)
        no=j+(i-1)*(nB+1); x=cno(no,1); y=cno(no,2);
        u0=zeros(3,1); [u0]=deslocamento(x,y,d);
        W(i,j)=u0(3,1); X(i,j)=x; Y(i,j)=y; Z(i,j)=0;
    end
end
surf(X,Y,W)
xlabel('x');
ylabel('y');
zlabel('W(x,y)');
title('Placa de Mindlin');


% Resultados normalizados
% u=[u; v; w; -Oy; Ox]
% w= deslocamento vertical normalizado
% Sxx (Syy e Sxy tambem)=[Sxx_inf_lam1 Sxx_sup_lam1; ...]
% Syz (Sxz tambem)=[Syz_lam1; Syz_lam2: ...]
[u,w,Sxx,Syy,Sxy,Syz,Sxz]=solucao(ladoA/2,ladoB/2,q0,d);
%w_centro=w
u
%Sxx_centro=Sxx
%Syy_centro=Syy
%Syy_centro=Syy(3,1)+(1/4-1/6)*(Syy(3,2)-Syy(3,1))/(1/3)
%[u,w,Sxx,Syy,Sxy,Syz,Sxz]=solucao(0,0,q0,d);
%Sxy_vertice=Sxy
%[u,w,Sxx,Syy,Sxy,Syz,Sxz]=solucao(0,ladoB/2,q0,d);
%Sxz_ladoB=Sxz
%[u,w,Sxx,Syy,Sxy,Syz,Sxz]=solucao(ladoA/2,0,q0,d);
%Syz_ladoB=Syz

% Condicao de contorno
%[u0]=deslocamento(ladoA/2,ladoB/2,d);
%w1=u0(3,1);
%[u0]=deslocamento(0,ladoB/2,d);
%w2=u0(3,1);
%erroCC=100*w2/w1;

% geo=[ladoA ladoB h nA nB];
% parametros=[tpMat emp tpLoad tpCC a];
% save exemplo1 geo parametros ep d w_centro Sxx_centro Syy_centro Sxz_ladoA Syz_ladoB