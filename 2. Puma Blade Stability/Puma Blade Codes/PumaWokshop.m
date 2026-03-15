%% Workshop 2 - Main Execution Script 
clc; clearvars;

total_time = tic;

%% --- PRECOMPUTING GEOMETRY AND BASE MATRICES (Done only once) ---
disp('--- Precomputing blade geometry and structural matrices... ---');

D = puma_data();
orden = 35;
n_node = 301;
x = linspace(D.flap_hinge_axis, D.blade_radius, n_node);

% Interpolations and mesh update
[EI, x_int] = w_interpolate(D.blade_flap_chord_stiffness(:, 1), D.blade_flap_chord_stiffness(:, 2), x);
[m, ~]      = w_interpolate(D.blade_mass(:, 1), D.blade_mass(:, 2), x);
[J, ~]      = w_interpolate(D.blade_torsional_inertia(:, 1), D.blade_torsional_inertia(:, 2), x);
[GJ, ~]     = w_interpolate(D.blade_torsional_stiffness(:, 1), D.blade_torsional_stiffness(:, 2), x);

x = x_int(:); EI = EI(:); m = m(:); J = J(:); GJ = GJ(:);

% Centrifugal force N(x)
N = zeros(size(x));
N(end-1) = m(end)/2 * (x(end)^2 - x(end-1)^2);
for i = length(N)-2:-1:1
    N(i) = N(i+1) + m(i)/2 * (x(i+1)^2 - x(i)^2);
end

% Pitch bearing node location
[~, Ix_pl] = min(abs(x - D.pitch_bearing));

% Isoparametric variables for Legendre polynomials
u_flap = 2 * (x - D.flap_hinge_axis) / (D.blade_radius - D.flap_hinge_axis) - 1;
u_tors = 2 * (x - D.pitch_bearing) / (D.blade_radius - D.pitch_bearing) - 1;
jacob_flap = 2 / (D.blade_radius - D.flap_hinge_axis);
jacob_tors = 2 / (D.blade_radius - D.pitch_bearing);

% Initialize shape functions
phi = zeros(length(x), orden); phi_p = zeros(length(x), orden);
psi = zeros(length(x), orden); psi_p = zeros(length(x), orden); psi_pp = zeros(length(x), orden);

% Generate shape functions using Legendre polynomials
for n = 1:orden
    P_flap = get_legendre(n-1, u_flap);     
    P_flap_p = get_legendre_p(n-1, u_flap); 
    P_flap_pp = get_legendre_pp(n-1, u_flap); 
    
    P_tors = get_legendre(n-1, u_tors);
    P_tors_p = get_legendre_p(n-1, u_tors);
    
    phi(:, n) = P_tors;
    phi_p(:, n) = P_tors_p * jacob_tors; 
    phi(1:Ix_pl-1, n) = 0; phi_p(1:Ix_pl-1, n) = 0;
    
    term_fix = (1 + u_flap); 
    psi(:, n) = term_fix .* P_flap;
    psi_p(:, n) = (term_fix .* P_flap_p + P_flap) * jacob_flap;
    psi_pp(:, n) = (term_fix .* P_flap_pp + 2 * P_flap_p) * jacob_flap^2;
end

% Base structural matrices
Kww = zeros(orden); Mww = zeros(orden);
Ktt_base = zeros(orden); Mtt = zeros(orden);
Mwt_base = zeros(orden); Mtw_base = zeros(orden); Kwt_base = zeros(orden);

for i=1:orden
    for j=1:orden
        Kww(i, j) = trapz(x, psi_pp(:,i).*psi_pp(:, j).*EI) + D.omega^2*trapz(x, psi_p(:, i).*psi_p(:, j).*N);
        Mww(i, j) = trapz(x, psi(:, i).*psi(:,j).*m);
        % Calculate Ktt WITHOUT the pitch link stiffness
        Ktt_base(i, j) = trapz(x, phi_p(:, i).*phi_p(:, j).*GJ) + D.omega^2*trapz(x, phi(:, i).*J.*phi(:, j));
        Mtt(i, j) = trapz(x, phi(:, i).*phi(:, j).*J);
        % Coupling matrices evaluated with unit offset (xt_n = 1)
        Mwt_base(i, j) = -trapz(x, psi(:, i).*phi(:, j).*m);
        Mtw_base(i, j) = -trapz(x, phi(:, i).*psi(:, j).*m);
        Kwt_base(i, j) = -D.omega^2*trapz(x, psi_p(:, i).*x.*m.*phi(:, j));
    end
end

% Package everything into structure S
S.D = D; S.orden = orden; S.x = x; 
S.Kww = Kww; S.Mww = Mww; S.Ktt_base = Ktt_base; S.Mtt = Mtt;
S.Mwt_base = Mwt_base; S.Mtw_base = Mtw_base; S.Kwt_base = Kwt_base;
S.phi_pl = phi(Ix_pl, :); % Pitch link row for fast matrix updates

% Trim shape functions for aerodynamic calculations
xaero_idx = find(x >= 2.08); 
S.xaero = x(xaero_idx);
S.psi_aero = psi(xaero_idx, :);
S.phi_aero = phi(xaero_idx, :);
S.psi_p_aero = psi_p(xaero_idx, :);

disp('--- Precomputations complete. Starting parameter sweep loops ---');

% Color map: Red = Unstable (0), Green = Stable (1)
color_map = [1 0 0; 0 1 0];

%% SECTION 3 (2 Modes)
disp('Running Section 3...');
tic;
link_stiff = linspace(0, 33032, 100);
cg_off = linspace(24, 41, 64);

stability_map = zeros(length(cg_off), length(link_stiff));
damping_map = zeros(length(cg_off), length(link_stiff), 2);

parfor i = 1:length(link_stiff)
    col_stab = zeros(length(cg_off), 1);
    col_damp = zeros(length(cg_off), 2);
    for j = 1:length(cg_off)
        [stab, damp] = PumaBlade_3_fast(link_stiff(i), cg_off(j), S);
        col_stab(j) = stab;
        col_damp(j, :) = damp;
    end
    stability_map(:, i) = col_stab;
    damping_map(:, i, :) = col_damp;
end
toc;

figure(1); clf;
imagesc(cg_off, link_stiff, stability_map'); 
colormap(color_map); set(gca, 'YDir', 'normal');  
xlabel('CG Offset (%)'); ylabel('Link Stiffness (N/m)'); 
title('Section 3: Stability Map'); colorbar('Ticks', [0, 1], 'TickLabels', {'Unstable', 'Stable'});


%% SECTION 4 (10 Modes)
disp('Running Section 4...');
tic;
link_stiff2 = linspace(0, 33032, 50);
cg_off2 = linspace(24, 40, 50);

stability_map2 = zeros(length(cg_off2), length(link_stiff2));
damping_map2 = zeros(length(cg_off2), length(link_stiff2), 10);

parfor i = 1:length(link_stiff2)
    col_stab2 = zeros(length(cg_off2), 1);
    col_damp2 = zeros(length(cg_off2), 10);
    for j = 1:length(cg_off2)
        [stab, damp] = PumaBlade_4_fast(link_stiff2(i), cg_off2(j), S);
        col_stab2(j) = stab;
        col_damp2(j, :) = damp;
    end
    stability_map2(:, i) = col_stab2;
    damping_map2(:, i, :) = col_damp2;
end
toc;

figure(2); clf;
imagesc(cg_off2, link_stiff2, stability_map2'); 
colormap(color_map); set(gca, 'YDir', 'normal');  
xlabel('CG Offset (%)'); ylabel('Link Stiffness (N/m)'); 
title('Section 4: Stability Map (10 Modes)'); colorbar('Ticks', [0, 1], 'TickLabels', {'Unstable', 'Stable'});


%% SECTION 5 (Vectorial Theodorsen)
disp('Running Section 5...');
tic;
link_stiff3 = linspace(0, 33032, 50);
cg_off3 = linspace(24, 40, 50);

stability_map3 = zeros(length(cg_off3), length(link_stiff3));
damping_map3 = zeros(length(cg_off3), length(link_stiff3), 2);

parfor i = 1:length(link_stiff3)
    col_stab3 = zeros(length(cg_off3), 1);
    col_damp3 = zeros(length(cg_off3), 2);
    for j = 1:length(cg_off3)
        [stab, damp] = PumaBlade_5_fast(link_stiff3(i), cg_off3(j), S);
        col_stab3(j) = stab;
        col_damp3(j, :) = damp;
    end
    stability_map3(:, i) = col_stab3;
    damping_map3(:, i, :) = col_damp3;
end
toc;

figure(3); clf;
imagesc(cg_off3, link_stiff3, stability_map3'); 
colormap(color_map); set(gca, 'YDir', 'normal');  
xlabel('CG Offset (%)'); ylabel('Link Stiffness (N/m)'); 
title('Section 5: Stability Map (Vectorial C(k))'); colorbar('Ticks', [0, 1], 'TickLabels', {'Unstable', 'Stable'});

%% Finish
fprintf('\n--- Total execution time: %.2f seconds ---\n', toc(total_time));


%% --- AUXILIARY FUNCTIONS ---
function P = get_legendre(n, u)
    if n==0, P = ones(size(u)); return; end
    if n==1, P = u; return; end
    % Recursion formula: (n)*P_n = (2n-1)u*P_{n-1} - (n-1)P_{n-2}
    P_prev2 = ones(size(u)); % P0
    P_prev1 = u;             % P1
    for k = 2:n
        P = ((2*k-1).*u.*P_prev1 - (k-1).*P_prev2) / k;
        P_prev2 = P_prev1;
        P_prev1 = P;
    end
end

function dP = get_legendre_p(n, u)
    if n==0, dP = zeros(size(u)); return; end
    % Derivative formula: (u^2-1)P'_n = n(u P_n - P_{n-1})
    P_n = get_legendre(n, u);
    P_nm1 = get_legendre(n-1, u);
    dP = (n * (u.*P_n - P_nm1)) ./ (u.^2 - 1);
    % Singularity correction at edges (-1 and 1)
    dP(abs(u)==1) = n*(n+1)/2 * sign(u(abs(u)==1)).^(n+1); 
end

function ddP = get_legendre_pp(n, u)
    if n==0 || n==1, ddP = zeros(size(u)); return; end
    % Legendre ODE: (1-u^2)P'' - 2uP' + n(n+1)P = 0
    P_n = get_legendre(n, u);
    dP_n = get_legendre_p(n, u);
    ddP = (2*u.*dP_n - n*(n+1)*P_n) ./ (1 - u.^2);
    % Numerical bound for edge singularities
    idx = abs(1-u.^2) < 1e-10;
    ddP(idx) = 0; 
end