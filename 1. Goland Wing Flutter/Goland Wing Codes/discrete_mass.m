function m_vec = discrete_mass(x, x_c, m)



    Diferencia = abs(x - x_c);
    [~, IndiceOrdenado] = sort(Diferencia);
    IndicesCercanos = IndiceOrdenado(1:2);
    Deltax = x(max(IndicesCercanos)) - x(min(IndicesCercanos));
    m_vec = zeros(size(x));
    m_vec(min(IndicesCercanos)) = (x_c - x(min(IndicesCercanos)))/Deltax*m/Deltax;
    m_vec(max(IndicesCercanos)) = (1- (x_c - x(min(IndicesCercanos)))/Deltax)*m/Deltax;
    % m_vec(min(IndicesCercanos)) = (x_c - x(min(IndicesCercanos)))*m;
    % m_vec(max(IndicesCercanos)) = (1- (x_c - x(min(IndicesCercanos))))*m;

