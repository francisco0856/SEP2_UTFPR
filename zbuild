function [Zbus] = zbuild(linedata)
    % This program forms the complex bus impedance matrix by the method
    % of building algorithm. Bus zero is taken as reference.
    % Copyright (c) 1998 by H. Saadat

    nl = linedata(:, 1);
    nr = linedata(:, 2);
    R = linedata(:, 3);
    X = linedata(:, 4);
    
    nbr = length(linedata(:, 1));
    nbus = max(max(nl), max(nr));
    
    % Define valores altos para R e X se forem inf
    for k = 1:nbr
        if R(k) == inf || X(k) == inf
            R(k) = 99999999;
            X(k) = 99999999;
        end
    end

    ZB = R + 1j * X;
    Zbus = zeros(nbus, nbus);
    tree = 0;  % Variável de controle para rastrear árvores construídas
    
    % Adicionando um ramo de um novo barramento ao barramento de referência (0)
    ntree = ones(nbr, 1);  % Inicializa o vetor de controle de árvore
    for I = 1:nbr
        if nl(I) == 0 || nr(I) == 0
            if nl(I) == 0
                n = nr(I);
            else
                n = nl(I);
            end
            if abs(Zbus(n, n)) == 0
                Zbus(n, n) = ZB(I);
                tree = tree + 1;
            else
                Zbus(n, n) = Zbus(n, n) * ZB(I) / (Zbus(n, n) + ZB(I));
            end
            ntree(I) = 2;
        end
    end

    % Adicionando um ramo de um novo barramento a um barramento existente
    while tree < nbus
        for n = 1:nbus
            nadd = 1;
            if abs(Zbus(n, n)) == 0
                for I = 1:nbr
                    if nadd == 1 && (nl(I) == n || nr(I) == n)
                        if nl(I) == n
                            k = nr(I);
                        else
                            k = nl(I);
                        end
                        if abs(Zbus(k, k)) ~= 0
                            for m = 1:nbus
                                if m ~= n
                                    Zbus(m, n) = Zbus(m, k);
                                    Zbus(n, m) = Zbus(m, k);
                                end
                            end
                            Zbus(n, n) = Zbus(k, k) + ZB(I);
                            tree = tree + 1;
                            nadd = 2;
                            ntree(I) = 2;
                        end
                    end
                end
            end
        end
    end

    % Adicionando uma ligação entre dois barramentos antigos
    for n = 1:nbus
        for I = 1:nbr
            if ntree(I) == 1 && (nl(I) == n || nr(I) == n)
                if nl(I) == n
                    k = nr(I);
                else
                    k = nl(I);
                end
                DM = Zbus(n, n) + Zbus(k, k) + ZB(I) - 2 * Zbus(n, k);
                DELZ = zeros(nbus, nbus);
                for jj = 1:nbus
                    AP = Zbus(jj, n) - Zbus(jj, k);
                    for kk = 1:nbus
                        AT = Zbus(n, kk) - Zbus(k, kk);
                        DELZ(jj, kk) = AP * AT / DM;
                    end
                end
                Zbus = Zbus - DELZ;
                ntree(I) = 2;
            end
        end
    end
end
