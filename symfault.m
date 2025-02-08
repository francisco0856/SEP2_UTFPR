function symfault(zdata, Zbus, V)
    % O programa symfault é projetado para a análise de faltas trifásicas
    % equilibradas em uma rede de sistema elétrico de potência. O programa
    % requer a matriz de impedância de barras Zbus. Zbus pode ser definida
    % pelo usuário, obtida pela inversão de Ybus ou determinada pelas
    % funções Zbus = zbuild(zdata) ou Zbus = zbuildpi(linedata, gendata, yload).
    % O programa solicita ao usuário que insira o número da barra em falta
    % e a impedância de falta Zf. As tensões pré-falta nas barras são definidas
    % pelo vetor reservado V. O vetor V pode ser definido ou retornado pelos
    % programas de fluxo de potência lfgauss, lfnewton, decouple ou perturb.
    % Se V não existir, as tensões pré-falta nas barras são automaticamente
    % definidas como 1,0 por unidade. O programa calcula a corrente total de falta,
    % as tensões pós-falta nas barras e as correntes nas linhas.
    %
    % Copyright (C) 1998 H. Saadat

    nl = zdata(:, 1);
    nr = zdata(:, 2);
    R = zdata(:, 3);
    X = zdata(:, 4);
    nc = size(zdata, 2);
    
    % Define BC dependendo do número de colunas em zdata
    if nc > 4
        BC = zdata(:, 5);
    else
        BC = zeros(length(zdata(:, 1)), 1);
    end
    
    ZB = R + 1j * X;
    nbr = length(zdata(:, 1));
    nbus = max(max(nl), max(nr));
    
    % Verifica se V foi definido
    if exist('V', 'var') == 1
        if length(V) == nbus
            V0 = V;
        else
            error('O vetor V precisa ter o mesmo comprimento que o número de barramentos.');
        end
    else
        V0 = ones(nbus, 1) + 1j * zeros(nbus, 1);
    end
    
    fprintf('\nAnálise de falta trifásica equilibrada \n');
    
    ff = 999;
    while ff > 0
        % Entrada do barramento com falta
        nf = input('Digite o número da barra com falta -> ');
        while nf <= 0 || nf > nbus
            fprintf('O número da barra com falta deve estar entre 1 e %g \n', nbus);
            nf = input('Digite o número da barra com falta -> ');
        end

        % Entrada da impedância de falta
        fprintf('\nInsira a impedância de falta Zf = R + j*X em ');
        Zf = input('forma complexa (para falta sólida, insira 0). Zf = ');
        fprintf('  \n');
        fprintf('Falta trifásica equilibrada na barra número %g\n', nf);

        % Cálculo da corrente de falta
        If = V0(nf) / (Zf + Zbus(nf, nf));
        Ifm = abs(If);
        Ifmang = angle(If) * 180 / pi;
        fprintf('Corrente total de falta = %8.4f por unidade \n\n', Ifm);
        
        % Exibe tensões durante a falta
        fprintf('Tensões nas barras durante a falta em por unidade \n\n');
        fprintf('     Barra     Tensão       Ângulo\n');
        fprintf('     Nº        Magnitude     Graus\n');

        for n = 1:nbus
            if n == nf
                Vf(nf) = V0(nf) * Zf / (Zf + Zbus(nf, nf));
                Vfm = abs(Vf(nf));
                angv = angle(Vf(nf)) * 180 / pi;
            else
                Vf(n) = V0(n) - V0(n) * Zbus(n, nf) / (Zf + Zbus(nf, nf));
                Vfm = abs(Vf(n));
                angv = angle(Vf(n)) * 180 / pi;
            end
            fprintf('   %4g %13.4f %13.4f\n', n, Vfm, angv);
        end

        fprintf('\n');

        % Cálculo das correntes de linha durante a falta
        fprintf('Correntes nas linhas para falta na barra número  %g\n\n', nf);
        fprintf('     De       Para     Corrente     Ângulo\n');
        fprintf('     Barra    Barra    Magnitude    Graus\n');

        for n = 1:nbus
            for I = 1:nbr
                if nl(I) == n || nr(I) == n
                    if nl(I) == n
                        k = nr(I);
                    else
                        k = nl(I);
                    end
                    if k == 0
                        Ink = (V0(n) - Vf(n)) / ZB(I);
                    else
                        Ink = (Vf(n) - Vf(k)) / ZB(I) + BC(I) * Vf(n);
                    end
                    Inkm = abs(Ink);
                    th = angle(Ink) * 180 / pi;
                    if real(Ink) > 0 || (real(Ink) == 0 && imag(Ink) < 0)
                        fprintf('%7g %10g %12.4f %12.4f\n', n, k, Inkm, th);
                    end
                end
            end
            if n == nf
                fprintf('%7g         F %12.4f %12.4f\n', n, Ifm, Ifmang);
            end
        end

        % Pergunta ao usuário se deseja uma nova localização de falta
        resp = '';
        while ~strcmp(resp, 'n') && ~strcmp(resp, 'N') && ~strcmp(resp, 'y') && ~strcmp(resp, 'Y')
            resp = input('Outra localização de falta? Digite ''y'' ou ''n'' entre aspas simples -> ', 's');
            if ~ismember(resp, {'y', 'Y', 'n', 'N'})
                fprintf('\n Resposta incorreta, tente novamente \n\n');
            end
        end
        if strcmp(resp, 'y') || strcmp(resp, 'Y')
            nf = 999;
        else
            ff = 0;
        end
    end
end
