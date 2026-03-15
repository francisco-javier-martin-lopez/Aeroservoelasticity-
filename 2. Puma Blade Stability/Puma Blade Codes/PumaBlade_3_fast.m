function [stability, damping_box] = PumaBlade_3_fast(link_stiff, cg_off, S)
    No = 2;
    
    %% 1.  STRUCTURAL RECONSTRUCTION
    % The CG offset is uniform, so we can treat it as a scalar
    xt_n_scalar = (cg_off - 25)/100 * S.D.blade_chord;
    
    % Update the link stiffness
    Ktt = S.Ktt_base + (S.phi_pl' * S.phi_pl) * link_stiff;
    
    % Scale the precalculated coupling matrices
    Mwt = S.Mwt_base * xt_n_scalar;
    Mtw = S.Mtw_base * xt_n_scalar;
    Kwt = S.Kwt_base * xt_n_scalar;
    Ktw = Kwt';
    
    M = [S.Mww, Mwt; Mtw, S.Mtt];
    K = [S.Kww, Kwt; Ktw,  Ktt];
    
    %% 2. Structural problem
    [Autovectores, Autovalores] = eig(K, M);
    aux = Autovectores;
    [Autovalores, sortng] = sort(sqrt(diag(Autovalores)), 'ascend');
    for j=1:length(Autovalores)
        Autovectores(:, j) = aux(:, sortng(j));
    end
    
    Autovectores = Autovectores(:, [1, 4]);
    Autovalores = Autovalores([1, 4]);
    
    % Modal projection of the aero shape functions
    psi = S.psi_aero * Autovectores(1:S.orden, :);
    psi_p = S.psi_p_aero * Autovectores(1:S.orden, :);
    phi = S.phi_aero * Autovectores(S.orden+1:end, :);
    
    M_gen = Autovectores' * M * Autovectores;
    M_gen(1,2) = 0; M_gen(2,1) = 0;
    K_gen = Autovectores' * K * Autovectores;
    [Autovectores_gen, ~] = eig(K_gen, M_gen);
    
    %% 3. Modal Aero Functions 
    ep = -0.5; b = S.D.blade_chord/2; rho = 1.225;
    Mach = 3/4*S.D.blade_radius*S.D.omega/340;
    Cla = 2*pi/sqrt(1 - Mach^2);
    
    Uvec = S.D.omega * S.xaero;
    
    WW_K = zeros(No); WW_M = zeros(No); WT_K = zeros(No); WT_M = zeros(No);
    TW_K = zeros(No); TW_M = zeros(No); TT_K = zeros(No); TT_M = zeros(No);
    WW_C = zeros(No); WT_C_Theta = zeros(No); WT_C_Omega = zeros(No);
    TW_C = zeros(No); TT_C_Theta = zeros(No); TT_C_Omega = zeros(No);
    
    for i=1:No
        for j=1:No
            WW_K(i, j) = trapz(S.xaero, psi(:, i).*psi(:, j).*Uvec.^2);
            WW_M(i, j) = trapz(S.xaero, psi(:, i).*psi(:, j));
            WT_K(i, j) = trapz(S.xaero, psi(:, i).*phi(:, j).*Uvec.^2);
            WT_M(i, j) = trapz(S.xaero, psi(:, i).*phi(:, j));
            TW_K(i, j) = trapz(S.xaero, phi(:, i).*psi(:, j).*Uvec.^2);
            TW_M(i, j) = trapz(S.xaero, phi(:, i).*psi(:, j));
            TT_K(i, j) = trapz(S.xaero, phi(:, i).*phi(:, j).*Uvec.^2); 
            TT_M(i, j) = trapz(S.xaero, phi(:, i).*phi(:, j)); 
            WW_C(i, j) = trapz(S.xaero, psi(:, i).*psi(:, j).*Uvec);
            WT_C_Theta(i, j) = trapz(S.xaero, psi(:, i).*phi(:, j).*Uvec);
            WT_C_Omega(i, j) = trapz(S.xaero, S.D.omega*Uvec.*psi(:, i).*psi_p(:, j)); 
            TW_C(i, j) = trapz(S.xaero, phi(:, i).*psi(:, j).*Uvec); 
            TT_C_Theta(i, j) = trapz(S.xaero, phi(:, i).*phi(:, j).*Uvec);
            TT_C_Omega(i, j) = trapz(S.xaero, S.D.omega*Uvec.*phi(:, i).*psi_p(:, j)); 
        end
    end
    
    %% 4. Flutter Resolution
    tolerancia = 10e-6;
    Amatrix = zeros(2*No,  2*No); Amatrix(1:No, (No+1):end) = eye(No);
    tol = 100; counter1 = 0; k_box = zeros(1, No); damping_box = zeros(1, No);
    
    for modes=1:No
        k = 4/3*Autovalores(modes)*b/S.D.omega/S.D.blade_radius;
        while tol>tolerancia
            Ck = besselh(1, 2, k)/(besselh(1, 2, k) + 1i*besselh(0, 2, k));
            
            Ma = rho*b^2*[pi,-pi*b*ep; ep*pi*b,-pi*b^2*(1/8+ep^2)];
            Ca = rho*b^2*[0,pi; 0,-pi*b*(1/2-ep)];
            Ka = zeros(2,2);
            Bw = rho*b*Cla*[1;b*(ep+1/2)]; Cw = [0, 1]; Chatw = [1,b*(1/2-ep)];
            
            Kaero = (Ka + Bw*Ck*Cw); Caero = (Ca+Bw*Ck*Chatw); Maero = Ma;
        
            KK = -Kaero(1,1)*WW_K + Caero(1, 2)*WT_C_Omega + Kaero(1,2)*WT_K -Kaero(2,1)*TW_K + Caero(2,2)*TT_C_Omega+ Kaero(2,2)*TT_K;
            CC = -Caero(1,1)*WW_C+ Caero(1,2)*WT_C_Theta -Caero(2,1)*TW_C+ Caero(2,2)*TT_C_Theta;
            MM = -Maero(1,1)*WW_M+ Maero(1,2)*WT_M -Maero(2,1)*TW_M + Maero(2,2)*TT_M;
    
            TM2 = M_gen - MM; TM1 = -CC; TM0 = K_gen - KK;
            Amatrix((No+1):end, 1:No) = -TM2\TM0; Amatrix((No+1):end, (No+1):end) = -TM2\TM1;
            
            [AVEC, z] = eig(Amatrix);
            if counter1 == 0
                Autovectores_track = AVEC((No+1):end, :); counter1 = 1;
            else 
                Autovectores_track = AVEC;
            end
            
            ScalarProducts = zeros(1, length(Autovectores_track));
            for tracking=1:length(Autovectores_track)
               ScalarProducts(tracking) = ps_complex(Autovectores_track(:, tracking), Autovectores_gen(:, modes));
            end
            [~, Index2] = max(ScalarProducts); 
            z_target = z(Index2, Index2); 
            
            omega_new = sqrt(imag(z_target)^2 + real(z_target)^2);
            diff1 = abs(Autovalores(modes) - omega_new); 
            Autovalores(modes) = omega_new;
            k = 4/3*Autovalores(modes)*b/S.D.omega/S.D.blade_radius;
            damping = real(z_target)/Autovalores(modes);
            tol = diff1;
            
            Autovectores_gen(1:No, modes) = AVEC(1:No, Index2);
            Autovectores_gen((No+1):(2*No), modes) = AVEC((No+1):end, Index2);
        end
        tol = 100; k_box(modes) = k; damping_box(modes) = damping; counter1 = 0;
    end
    
    stability = true;
    for kkk =1:No
        if damping_box(kkk)>1e-3 && stability 
            stability = false;
        end
    end
end