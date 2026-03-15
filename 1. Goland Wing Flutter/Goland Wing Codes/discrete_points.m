function close_points = discrete_points(x, x0, xf)

    Distancia0 = abs(x-x0);
    Distanciaf = abs(x-xf);
    [~, Index0] = sort(Distancia0);
    [~, Indexf] = sort(Distanciaf);
    close_points(1) = Index0(1);
    close_points(2) = Indexf(1);
