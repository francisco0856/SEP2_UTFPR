% O programa llfault é projetado para a análise de faltas bifásicas
% em uma rede de sistema elétrico de potência. O programa requer
% as matrizes de impedância de barras de sequência positiva e negativa,
% Zbus1 e Zbus2. As matrizes de impedância de barras podem ser definidas
% pelo usuário, obtidas pela inversão de Ybus ou determinadas pelas
% funções Zbus = zbuild(zdata) ou Zbus = zbuildpi(linedata, gendata, yload).
% O programa solicita ao usuário que insira o número da barra com falta
% e a impedância de falta Zf. As tensões pré-falta nas barras são definidas
% pelo vetor reservado V. O vetor V pode ser definido ou retornado pelos
% programas de fluxo de potência lfgauss, lfnewton, decouple ou perturb.
% Se V não existir, as tensões pré-falta nas barras são automaticamente
% definidas como 1,0 por unidade. O programa calcula a corrente total de falta,
% as tensões nas barras e as correntes nas linhas durante a falta.
%
% Copyright (C) 1998 Hadi Saadat

function llfault(zdata1, Zbus1, zdata2, Zbus2, V)
if exist('zdata2') ~= 1
    zdata2 = zdata1;
else, end
if exist('Zbus2') ~= 1
    Zbus2 = Zbus1;
else, end

nl = zdata1(:,1); nr = zdata1(:,2);
R1 = zdata1(:,3); X1 = zdata1(:,4);
R2 = zdata2(:,3); X2 = zdata2(:,4);
ZB1 = R1 + j*X1;  
ZB2 = R2 + j*X2;
nbr = length(zdata1(:,1)); 
nbus = max(max(nl), max(nr));

if exist('V', 'var') == 1
    if length(V) == nbus
        V0 = V;
    else, end
else
    V0 = ones(nbus, 1) + j*zeros(nbus, 1);
end

fprintf('\nAnálise de falta bifásica \n')
ff = 999;
while ff > 0
    nf = input('Digite o número da barra com falta -> ');
    while nf <= 0 || nf > nbus
        fprintf('O número da barra com falta deve estar entre 1 e %g \n', nbus)
        nf = input('Digite o número da barra com falta -> ');
    end
    fprintf('\nInsira a impedância de falta Zf = R + j*X em ')
    Zf = input('forma complexa (para falta sólida, insira 0). Zf = ');
    fprintf('  \n')
    fprintf('Falta bifásica na barra número %g\n', nf)

    a = cos(2*pi/3) + j*sin(2*pi/3);
    sctm = [1   1   1; 1 a^2  a; 1 a  a^2];
    Ia0 = 0;
    Ia1 = V0(nf) / (Zbus1(nf,nf) + Zbus2(nf,nf) + Zf); 
    Ia2 = -Ia1;
    I012 = [Ia0; Ia1; Ia2];
    Ifabc = sctm * I012;
    Ifabcm = abs(Ifabc);

    fprintf('Corrente total de falta = %9.4f por unidade\n\n', Ifabcm(2))
    fprintf('Tensões nas barras durante a falta em por unidade \n\n')
    fprintf('     Barra    -------Magnitude da Tensão-------  \n')
    fprintf('     Nº       Fase a     Fase b     Fase c  \n')

    for n = 1:nbus
        Vf0(n) = 0;
        Vf1(n) = V0(n) - Zbus1(n, nf) * Ia1;
        Vf2(n) = 0 - Zbus2(n, nf) * Ia2;
        Vabc = sctm * [Vf0(n); Vf1(n); Vf2(n)];
        Va(n) = Vabc(1); 
        Vb(n) = Vabc(2); 
        Vc(n) = Vabc(3);
        fprintf(' %5g', n)
        fprintf(' %11.4f', abs(Va(n))), fprintf(' %11.4f', abs(Vb(n)))
        fprintf(' %11.4f\n', abs(Vc(n)))
    end

    fprintf('\n')
    fprintf('Correntes nas linhas para falta na barra número %g\n\n', nf)
    fprintf('     De       Para      -----Magnitude da Corrente na Linha----  \n')
    fprintf('     Barra    Barra     Fase a     Fase b     Fase c  \n')

    for n = 1:nbus
        for I = 1:nbr
            if nl(I) == n || nr(I) == n
                if nl(I) == n
                    k = nr(I);
                else
                    k = nl(I);
                end
                if k ~= 0
                    Ink0(n, k) = 0;
                    Ink1(n, k) = (Vf1(n) - Vf1(k)) / ZB1(I);
                    Ink2(n, k) = (Vf2(n) - Vf2(k)) / ZB2(I);

                    Inkabc = sctm * [Ink0(n, k); Ink1(n, k); Ink2(n, k)];
                    Inkabcm = abs(Inkabc);
                    th = angle(Inkabc);

                    fprintf('%7g %10g %11.4f %11.4f %11.4f\n', n, k, abs(Inkabc(1)), abs(Inkabc(2)), abs(Inkabc(3)));
                end
            end
        end
    end

    resp = '';
    while ~strcmp(resp, 'n') && ~strcmp(resp, 'N') && ~strcmp(resp, 'y') && ~strcmp(resp, 'Y')
        resp = input('Outra localização de falta? Digite ''y'' ou ''n'' entre aspas simples -> ', 's');
        if ~ismember(resp, {'y', 'Y', 'n', 'N'})
            fprintf('\nResposta incorreta, tente novamente.\n\n');
        end
    end
    if strcmp(resp, 'y') || strcmp(resp, 'Y')
        nf = 999;
    else
        ff = 0;
    end
end
