% O programa dlgfault é projetado para a análise de faltas bifásicas para 
% terra em uma rede de sistema elétrico de potência. O programa requer 
% as matrizes de impedância de barras de sequência positiva, negativa ou zero, 
% Zbus1, Zbus2 e Zbus0. As matrizes de impedância de barras podem ser definidas 
% pelo usuário, obtidas pela inversão de Ybus ou determinadas pelas funções 
% Zbus = zbuild(zdata) ou Zbus = zbuildpi(linedata, gendata, yload).
% O programa solicita ao usuário que insira o número da barra em falta 
% e a impedância de falta Zf. As tensões pré-falta nas barras são definidas 
% pelo vetor reservado V. O vetor V pode ser definido ou retornado pelos 
% programas de fluxo de potência lfgauss, lfnewton, decouple ou perturb.
% Se V não existir, as tensões pré-falta nas barras são automaticamente 
% definidas como 1,0 por unidade. O programa calcula a corrente total de falta, 
% as tensões nas barras e as correntes nas linhas durante a falta.
%
% Copyright (C) 1998 Hadi Saadat

function dlgfault(zdata0, Zbus0, zdata1, Zbus1, zdata2, Zbus2, V)

if exist('zdata2') ~= 1
    zdata2 = zdata1;
else, end
if exist('Zbus2') ~= 1
    Zbus2 = Zbus1;
else, end

nl = zdata1(:,1); nr = zdata1(:,2);
nl0 = zdata0(:,1); nr0 = zdata0(:,2);
nbr = length(zdata1(:,1)); nbus = max(max(nl), max(nr));
nbr0 = length(zdata0(:,1));
R0 = zdata0(:,3); X0 = zdata0(:,4);
R1 = zdata1(:,3); X1 = zdata1(:,4);
R2 = zdata2(:,3); X2 = zdata2(:,4);

for k = 1:nbr0
    if R0(k) == inf || X0(k) == inf
        R0(k) = 99999999; X0(k) = 999999999;
    else, end
end
ZB1 = R1 + j*X1;  ZB0 = R0 + j*X0;
ZB2 = R2 + j*X2;

if exist('V') == 1
    if length(V) == nbus
        V0 = V;
    else, end
else
    V0 = ones(nbus, 1) + j*zeros(nbus, 1);
end

fprintf('\nAnálise de falta bifásica para terra \n')
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
    fprintf('Falta bifásica para terra na barra número %g\n', nf)
    a = cos(2*pi/3) + j*sin(2*pi/3);
    sctm = [1   1   1; 1 a^2  a; 1 a  a^2];

    Z11 = Zbus2(nf, nf)*(Zbus0(nf, nf) + 3*Zf)/(Zbus2(nf, nf) + Zbus0(nf, nf) + 3*Zf);
    Ia1 = V0(nf)/(Zbus1(nf,nf) + Z11);
    Ia2 = -(V0(nf) - Zbus1(nf, nf)*Ia1)/Zbus2(nf,nf);
    Ia0 = -(V0(nf) - Zbus1(nf, nf)*Ia1)/(Zbus0(nf,nf) + 3*Zf);
    I012 = [Ia0; Ia1; Ia2];
    Ifabc = sctm*I012; Ifabcm = abs(Ifabc);
    Ift = Ifabc(2) + Ifabc(3);
    Iftm = abs(Ift);

    fprintf('Corrente total de falta = %9.4f por unidade\n\n', Iftm)
    fprintf('Tensões nas barras durante a falta em por unidade \n\n')
    fprintf('     Barra    -------Magnitude da Tensão-------  \n')
    fprintf('     Nº       Fase a     Fase b     Fase c  \n')

    for n = 1:nbus
        Vf0(n) = 0 - Zbus0(n, nf)*Ia0;
        Vf1(n) = V0(n) - Zbus1(n, nf)*Ia1;
        Vf2(n) = 0 - Zbus2(n, nf)*Ia2;
        Vabc = sctm*[Vf0(n); Vf1(n); Vf2(n)];
        Va(n) = Vabc(1); Vb(n) = Vabc(2); Vc(n) = Vabc(3);
        fprintf(' %5g', n)
        fprintf(' %11.4f', abs(Va(n))), fprintf(' %11.4f', abs(Vb(n)))
        fprintf(' %11.4f\n', abs(Vc(n)))
    end
    fprintf('  \n')
    fprintf('Correntes nas linhas para falta na barra número  %g\n\n', nf)
    fprintf('     De       Para      -----Magnitude da Corrente na Linha----  \n')
    fprintf('     Barra    Barra     Fase a     Fase b     Fase c  \n')

    % Continuação do código...
end
