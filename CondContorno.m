function [CCu,CCo]=CondContorno(tpCC);

% Condicoes de Contorno Essenciais - Deslocamentos Prescritos

format long

% CCu=[el face valor dof_u;...]
% CCo=[el face valor dof_o;...]

% dof_u: grau de liberdade de deslocamento prescrito (u=1; v=2; w=3)
% dof_o: grau de liberdade de deslocamento prescrito (Oy=alfa=1; -Ox=beta=2)
load geometria nA nB
nu=0; no=0;
if tpCC==1 % Simplesmente suportado: SS-1 (Reddy)
    % suporte em x=0 (v=w=beta=0)
    for i=1:nB
        el=2+(i-1)*nA*2;
        nu=nu+1; CCu(nu,:)=[el 2 0 2];
        nu=nu+1; CCu(nu,:)=[el 2 0 3];
        no=no+1; CCo(no,:)=[el 2 0 2];
    end
    % suporte em y=0 (u=w=alfa=0)
    for i=1:nA
        el=2*i-1;
        nu=nu+1; CCu(nu,:)=[el 1 0 1];
        nu=nu+1; CCu(nu,:)=[el 1 0 3];
        no=no+1; CCo(no,:)=[el 1 0 1];
    end
    % suporte em x=ladoA (v=w=beta=0)
    for i=1:nB
        el=i*(2*nA)-1;
        nu=nu+1; CCu(nu,:)=[el 2 0 2];
        nu=nu+1; CCu(nu,:)=[el 2 0 3];
        no=no+1; CCo(no,:)=[el 2 0 2];

    end
    % suporte em y=ladoB (u=w=alfa=0)
    for i=1:nA
        el=(nB-1)*2*nA+2*i;
        nu=nu+1; CCu(nu,:)=[el 1 0 1];
        nu=nu+1; CCu(nu,:)=[el 1 0 3];
        no=no+1; CCo(no,:)=[el 1 0 1];
    end
end
if tpCC==2 % Simplesmente suportado: SS-2 (Reddy)
    % suporte em x=0 (u=w=beta=0)
    for i=1:nB
        el=2+(i-1)*nA*2;
        nu=nu+1; CCu(nu,:)=[el 2 0 1];
        nu=nu+1; CCu(nu,:)=[el 2 0 3];
        no=no+1; CCo(no,:)=[el 2 0 2];
    end
    % suporte em y=0 (v=w=alfa=0)
    for i=1:nA
        el=2*i-1;
        nu=nu+1; CCu(nu,:)=[el 1 0 2];
        nu=nu+1; CCu(nu,:)=[el 1 0 3];
        no=no+1; CCo(no,:)=[el 1 0 1];

    end
    % suporte em x=ladoA (u=w=beta=0)
    for i=1:nB
        el=i*(2*nA)-1;
        nu=nu+1; CCu(nu,:)=[el 2 0 1];
        nu=nu+1; CCu(nu,:)=[el 2 0 3];
        no=no+1; CCo(no,:)=[el 2 0 2];
    end
    % suporte em y=ladoB (v=w=alfa=0)
    for i=1:nA
        el=(nB-1)*2*nA+2*i;
        nu=nu+1; CCu(nu,:)=[el 1 0 2];
        nu=nu+1; CCu(nu,:)=[el 1 0 3];
        no=no+1; CCo(no,:)=[el 1 0 1];
    end
end
if tpCC==3 % Suportada em y=0 e y=ladoB conforme SS-1
    % suporte em y=0 (u=w=alfa=0)
    for i=1:nA
        el=2*i-1;
        nu=nu+1; CCu(nu,:)=[el 1 0 1];
        nu=nu+1; CCu(nu,:)=[el 1 0 3];
        no=no+1; CCo(no,:)=[el 1 0 1];
    end
    % suporte em y=ladoB (u=w=alfa=0)
    for i=1:nA
        el=(nB-1)*2*nA+2*i;
        nu=nu+1; CCu(nu,:)=[el 1 0 1];
        nu=nu+1; CCu(nu,:)=[el 1 0 3];
        no=no+1; CCo(no,:)=[el 1 0 1];
    end
end
if tpCC==4 % Suportada em x=0 e x=ladoA conforme SS-2
    % suporte em x=0 (u=w=beta=0)
    for i=1:nB
        el=2+(i-1)*nA*2;
        nu=nu+1; CCu(nu,:)=[el 2 0 1];
        nu=nu+1; CCu(nu,:)=[el 2 0 3];
        no=no+1; CCo(no,:)=[el 2 0 2];
    end
    % suporte em x=ladoA (u=w=beta=0)
    for i=1:nB
        el=i*(2*nA)-1;
        nu=nu+1; CCu(nu,:)=[el 2 0 1];
        nu=nu+1; CCu(nu,:)=[el 2 0 3];
        no=no+1; CCo(no,:)=[el 2 0 2];
    end
end
if tpCC==5 % suportada em x=0 (SS-1) e y=0 (SS-2)
    % suporte em x=0 (v=w=beta=0)
    for i=1:nB
        el=2+(i-1)*nA*2;
        nu=nu+1; CCu(nu,:)=[el 2 0 2];
        nu=nu+1; CCu(nu,:)=[el 2 0 3];
        no=no+1; CCo(no,:)=[el 2 0 2];
    end
        % suporte em y=0 (v=w=alfa=0)
    for i=1:nA
        el=2*i-1;
        nu=nu+1; CCu(nu,:)=[el 1 0 2];
        nu=nu+1; CCu(nu,:)=[el 1 0 3];
        no=no+1; CCo(no,:)=[el 1 0 1];
    end
end
if tpCC==6 % Engastado em x=0
    % suporte em x=0 (u=v=w=alfa=beta=0)
    for i=1:nB
        el=2+(i-1)*nA*2;
        nu=nu+1; CCu(nu,:)=[el 2 0 1];
        %nu=nu+1; CCu(nu,:)=[el 2 0 2];
        nu=nu+1; CCu(nu,:)=[el 2 0 3];
        %no=no+1; CCo(no,:)=[el 2 0 1];
        no=no+1; CCo(no,:)=[el 2 0 2];
    end
        % suporte em y=0 (v=w=alfa=0)
    for i=1:nA
        el=2*i-1;
        nu=nu+1; CCu(nu,:)=[el 1 0 2];
        nu=nu+1; CCu(nu,:)=[el 1 0 3];
        no=no+1; CCo(no,:)=[el 1 0 1];
    end
end
if tpCC==7 % simplesmente apoiada (1/4 da placa)
    % suporte em x=0 (w=Ox=0)
    for i=1:nB
        el=2+(i-1)*nA*2;
        nu=nu+1; CCu(nu,:)=[el 2 0 3];
        no=no+1; CCo(no,:)=[el 2 0 2];
    end
    % suporte em y=0 (w=Oy=0)
    for i=1:nA
        el=2*i-1;
        nu=nu+1; CCu(nu,:)=[el 1 0 3];
        no=no+1; CCo(no,:)=[el 1 0 1];
    end
    % suporte em x=ladoA (u=0)
    for i=1:nB
        el=i*(2*nA)-1;
        nu=nu+1; CCu(nu,:)=[el 2 0 1];
    end
    % suporte em y=ladoB (v=0)
    for i=1:nA
        el=(nB-1)*2*nA+2*i;
        nu=nu+1; CCu(nu,:)=[el 1 0 2];
    end
end
if tpCC==8 % simplesmente apoiada (toda a placa)
    CCo=-1;
    % suporte em x=0 (w=Ox=0)
    for i=1:nB
        el=2+(i-1)*nA*2;
        nu=nu+1; CCu(nu,:)=[el 2 0 1];
        nu=nu+1; CCu(nu,:)=[el 2 0 3];
    end
    % suporte em y=0 (w=Oy=0)
    for i=1:nA
        el=2*i-1;
        nu=nu+1; CCu(nu,:)=[el 1 0 2];
        nu=nu+1; CCu(nu,:)=[el 1 0 3];
    end
    % suporte em x=ladoA (u=0)
    for i=1:nB
        el=i*(2*nA)-1;
        nu=nu+1; CCu(nu,:)=[el 2 0 3];
    end
    % suporte em y=ladoB (v=0)
    for i=1:nA
        el=(nB-1)*2*nA+2*i;
        nu=nu+1; CCu(nu,:)=[el 1 0 3];
    end
end
if tpCC==9 % Totalmente engastada
    % suporte em x=0
    for i=1:nB
        el=2+(i-1)*nA*2;
        nu=nu+1; CCu(nu,:)=[el 2 0 1];
        nu=nu+1; CCu(nu,:)=[el 2 0 2];
        nu=nu+1; CCu(nu,:)=[el 2 0 3];
        no=no+1; CCo(no,:)=[el 2 0 1];
        no=no+1; CCo(no,:)=[el 2 0 2];
    end
    % suporte em y=0
    for i=1:nA
        el=2*i-1;
        %nu=nu+1; CCu(nu,:)=[el 1 0 1];
        %nu=nu+1; CCu(nu,:)=[el 1 0 2];
        %nu=nu+1; CCu(nu,:)=[el 1 0 3];
        %no=no+1; CCo(no,:)=[el 1 0 1];
        %no=no+1; CCo(no,:)=[el 1 0 2];
    end
    % suporte em x=ladoA
    for i=1:nB
        el=i*(2*nA)-1;
        nu=nu+1; CCu(nu,:)=[el 2 0 1];
        nu=nu+1; CCu(nu,:)=[el 2 0 2];
        nu=nu+1; CCu(nu,:)=[el 2 0 3];
        no=no+1; CCo(no,:)=[el 2 0 1];
        no=no+1; CCo(no,:)=[el 2 0 2];

    end
    % suporte em y=ladoB
    for i=1:nA
        el=(nB-1)*2*nA+2*i;
        nu=nu+1; CCu(nu,:)=[el 1 0 1];
        nu=nu+1; CCu(nu,:)=[el 1 0 2];
        nu=nu+1; CCu(nu,:)=[el 1 0 3];
        no=no+1; CCo(no,:)=[el 1 0 1];
        no=no+1; CCo(no,:)=[el 1 0 2];
    end
end
