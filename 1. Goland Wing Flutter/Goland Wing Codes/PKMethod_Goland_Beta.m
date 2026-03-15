clc
clear all


%% Parameters
orden = 5;
EJ = 9.77e6;
m = 35.75;
L = 6.096;
GJ = 9.876e5;
Jp = 8.65;
chord = 1.829;
xt_n = 0.1*chord;
ep = -0.34;
b = chord/2;
Cla = 2*pi;
rho = 1.225;
mf = 8.92/1.8;
cf = 0.45;
K_beta = 6.48e3;
x0 = 3.35;
xf = L-0.95;
sbeta = 0.67*chord-0.45;
Ibeta = 0.09/1.8 + mf*(cf/2)^2;
cp = (chord*0.5-0.45)/(chord/2);
x_c = (xf+x0)/2;
n_element = 18;
n_node = n_element+1;
x = linspace(0, L, n_node);
%% Rayleigh-Ritz
for n = 1:orden
    phi(:, n) = sin(n*pi/(2*L)*x);
    psi(:, n) = 1-cos(n*pi/(2*L)*x);
    psi_pp(:, n) = (n*pi/2/L)^2*cos(n*pi/2/L*x);
    phi_p(:,n) = (n*pi/2/L)*cos(n*pi/(2*L)*x);
end

%% K M matrices
m_vec = discrete_mass(x, x_c, mf*1.8)';
Ibeta_vec = discrete_mass(x, x_c, 0.09)';

for i=1:orden
    for j=1:orden
    Kww(i, j) = trapz(x, psi_pp(:,i).*psi_pp(:, j)*EJ);
    Mww(i, j) = trapz(x,psi(:, i).*psi(:,j)*m) + trapz(x,psi(:, i).*psi(:,j).*m_vec);
    Ktt(i, j) = trapz(x,phi_p(:, i).*phi_p(:, j)*GJ);
    Mtt(i, j) = trapz(x,phi(:, i).*phi(:, j)*Jp) + trapz(x,phi(:, i).*phi(:, j).*(Ibeta_vec + m_vec*(sbeta+cf/2)^2));
    Mwt(i, j) = -trapz(x,psi(:, i).*phi(:, j)*m*xt_n) - trapz(x,psi(:, i).*phi(:, j).*m_vec*(sbeta + cf/2));
    Mtw(i, j) = -trapz(x,phi(:, i).*psi(:, j)*m*xt_n) - trapz(x,phi(:, i).*psi(:, j).*m_vec*(sbeta + cf/2));
    end
   Mwb(i, 1) = -trapz(x,psi(:, i).*m_vec*cf/2);
   Mbw(1, i) = -trapz(x,psi(:, i).*m_vec*cf/2);
   Mtb(i, 1) = trapz(x,phi(:, i).*(Ibeta_vec + (sbeta + cf/2)*cf/2*m_vec));
   Mbt(1, i) = trapz(x,phi(:, i).*(Ibeta_vec + (sbeta + cf/2)*cf/2*m_vec));
end
Kbb = K_beta;
Mbb = Ibeta*1.8;
M = [Mww, Mwt, Mwb;Mtw, Mtt, Mtb;Mbw, Mbt, Mbb];
K = [Kww, zeros(length(Kww)), zeros(length(Kww),1); zeros(length(Kww)), Ktt, zeros(length(Kww),1); zeros(length(Kww),1)', zeros(length(Kww),1)', Kbb ];
%% Structural
[Autovectores, Autovalores] = eig(M\K);

aux = Autovectores;
[Autovalores, sortng] = sort(sqrt(diag(Autovalores)), 'ascend');
for j=1:length(Autovalores)
    Autovectores(:, j) = aux(:, sortng(j));
end
aux = Autovectores;
aux2 = Autovalores;
%% Theodorsen parameters
T = TheodorsenCoefficients(cp,ep);

%% Aerodynamic shape functions
Ix = discrete_points(x, x0, xf);
for i=1:orden
    for j=1:orden
        WW(i, j) = trapz(x, psi(:, i).*psi(:, j));
        TT(i, j) = trapz(x, phi(:, i).*phi(:, j));
        WT(i, j) = trapz(x, psi(:, i).*phi(:, j));
        TW(i, j) = trapz(x, phi(:, i).*psi(:, j));
    end
    WB(i, 1) = trapz(x(Ix(1):Ix(2)), psi((Ix(1):Ix(2)), i));
    BW(1, i) = trapz(x(Ix(1):Ix(2)), psi((Ix(1):Ix(2)), i));
    TB(i, 1) = trapz(x(Ix(1):Ix(2)), phi((Ix(1):Ix(2)), i));
    BT(1, i) = trapz(x(Ix(1):Ix(2)), phi((Ix(1):Ix(2)), i));
end
BB = xf-x0;
phi = zeros(orden, 1);
psi = zeros(orden, 1);
%% Solving parameters
tolerancia = 10e-4;
Uvec = linspace(0.5, 300, 300);
Amatrix = zeros(4*length(phi)+2,  4*length(phi)+2);
Amatrix (1:(2*length(phi)+1), (2*length(phi)+2):end) = eye(2*length(psi)+1);
tol = 100;
counter1 = 0;
%% Flutter solver

for modes=1:2*length(psi)
    for jj=1:length(Uvec)

        U = Uvec(jj);
        k = Autovalores(modes)*b/U;
        while tol>tolerancia

            Ck = besselh(1, 2, k)/(besselh(1, 2, k) + 1i*besselh(0, 2, k));

            Ma = rho*b^2*[pi,-pi*b*ep,-T(1)*b;
                        ep*pi*b,-pi*b^2*(1/8+ep^2),(T(7)+(0.5079-ep)*T(1))*b^2;
                        T(1)*b,-2*T(13)*b^2,1/pi*T(3)*b^2];
            
            Ca = rho*b^2*[0,pi,-T(4);
                            0,-pi*b*(1/2-ep),-b*(T(1)-T(8)-(0.5079-ep)*T(4)+T(11)/2);
                            0, -b*(-2*T(9)-T(1)+T(4)*(ep-1/2)),b/(2*pi)*T(4)*T(11)];
            Ka = rho*b^2*[0,0,0;
                          0,0,-(T(4)+T(10));
                          0,0,-(1/pi)*(T(5)-T(4)*T(10))];
            
            Bw = rho*b*Cla*[1;b*(ep+1/2);-b*T(12)/(2*pi)];
            
            Cw = [0, 1, T(10)/pi];
            Chatw = [1,b*(1/2-ep),b*T(11)/(2*pi)];
            
            Kaero = (Ka + Bw*Ck*Cw);
            Caero = (Ca+Bw*Ck*Chatw);
            Maero = Ma;
        
            KK = U^2*[-Kaero(1,1)*WW, Kaero(1,2)*WT, Kaero(1,3)*WB; -Kaero(2,1)*TW, Kaero(2,2)*TT, Kaero(2,3)*TB; -Kaero(3, 1)*BW, Kaero(3,2)*BT, Kaero(3,3)*BB];
            CC = U*[-Caero(1,1)*WW, Caero(1,2)*WT, Caero(1,3)*WB; -Caero(2,1)*TW, Caero(2,2)*TT, Caero(2,3)*TB; -Caero(3, 1)*BW, Caero(3,2)*BT, Caero(3,3)*BB];
            MM = [-Maero(1,1)*WW, Maero(1,2)*WT, Maero(1,3)*WB; -Maero(2,1)*TW, Maero(2,2)*TT, Maero(2,3)*TB; -Maero(3, 1)*BW, Maero(3,2)*BT, Maero(3,3)*BB];
    
            TM2 = M -MM;

            TM1 = -CC;
            TM0 = K-KK;
            Amatrix((2*length(psi)+2):end, 1:(2*length(psi)+1)) = -TM2\TM0;
            Amatrix((2*length(psi)+2):end, (2*length(psi)+2):end) = -TM2\TM1;
            
            [AVEC, z] = eig(Amatrix);

            if counter1 == 0
                Autovectores_track = AVEC((2*length(phi)+2):end, :);
                counter1 = 1;
           else 
               Autovectores_track = AVEC;
           end
           % Tracking of the eigenvector -> maximum scalar product of
           % complex numbers

           for tracking=1:length(Autovectores)
               ScalarProducts(tracking) = ps_complex(Autovectores_track(:, tracking), Autovectores(:, modes));
           end

           [aux1, Index2] = max(ScalarProducts); 

           z = diag(z);
           z = z(Index2); % Eigenvalue of the mode we are studying

           diff1= abs(Autovalores(modes) - sqrt(imag(z)^2 + real(z)^2)); % Check for convergence of omega between iterations

           Autovalores(modes) = sqrt(imag(z)^2 + real(z)^2);
           k = Autovalores(modes)*b/U;
           damping = real(z)/Autovalores(modes);
           tol = diff1;
           % Redefinition of eigenvectors for next iteration of convergence
           % (and velocity if last convg iteration)
           Autovectores(1:(2*length(phi)+1), modes) = AVEC(1:(2*length(phi)+1), Index2);
           Autovectores((2*length(phi)+2):(4*length(phi)+2), modes) = AVEC((2*length(phi)+2):end, Index2);
    
        end
        tol = 100;
        k_box (jj, modes) = k;
        damping_box (jj, modes) = damping;
    end
    counter1 = 0;
end

flutternotfound = true;
for jj=1:length(Uvec)
    for kkk =1:2*length(phi)
        if damping_box(jj, kkk)>1e-5 && flutternotfound
            jjj = jj;
            kkkk = kkk;
            flutternotfound = false;
        end
    end
end
if flutternotfound 
    U_flutter = NaN;
elseif jjj==1
    U_flutter = Uvec(1);
else

m = (real(damping_box(jjj-1, kkkk)) - real(damping_box(jjj, kkkk)))/(Uvec(jjj-1) - Uvec(jjj));
n = real (damping_box(jjj-1, kkkk)) - m*Uvec(jjj-1);
U_flutter = -n/m
end

%% Control Reversal 
CRnotfound = true;
DIVnotfound = true;

% Increase the resolution of the velocity vector for better precision
Uvec2 = linspace(0, 300, 500); 

% 1. DISCRETE CALCULATION OF THETA INT
% Reconstruct phi_modes (since 'phi' was overwritten with zeros earlier)
phi_modes = zeros(length(x), orden);
ThetaInt = zeros(1, orden);

for n = 1:orden
    phi_modes(:, n) = sin(n*pi/(2*L)*x');
    % Numerical integration using the trapezoidal rule
    ThetaInt(n) = trapz(x, phi_modes(:, n));
end

% Pre-allocate memory for control effectiveness (MATLAB best practices)
Epsilon = zeros(1, length(Uvec2));

% 2. STATIC RESOLUTION LOOP
for jj = 1:length(Uvec2)
    U = Uvec2(jj);
    Ck1 = 1; % Theodorsen function = 1 for static case (Reversal/Divergence)
    
    % Static aerodynamic matrices
    Kaero = (Ka + Bw*Ck1*Cw);
    
    % Aerodynamic stiffness matrix (KK) for current velocity U
    KK = U^2 * [-Kaero(1,1)*WW,  Kaero(1,2)*WT,  Kaero(1,3)*WB; 
                -Kaero(2,1)*TW,  Kaero(2,2)*TT,  Kaero(2,3)*TB; 
                -Kaero(3,1)*BW,  Kaero(3,2)*BT,  Kaero(3,3)*BB];
    
    % Static system matrix (Structural stiffness - Aerodynamic stiffness)
    TM0_CR = K - KK;
    
    % Generalized force vector for unit control deflection (delta = 1)
    F_mando = [zeros(2*orden, 1); K_beta];
    
    % Solve for static deformations (gamma)
    gamma = TM0_CR \ F_mando;
    
    % 3. CONTROL EFFECTIVENESS CALCULATION (Epsilon)
    % Lift due to elastic torsion along the span
    Ltheta_vec = rho * b * U^2 * Cla * ThetaInt;
    
    % Lift due to flap deflection itself
    Lbeta = rho * b * U^2 * Cla * (T(10)/pi) * (xf-x0);
    
    % Total effectiveness: Flap lift + lift loss due to torsion
    % Note: gamma(orden+1 : 2*orden) extracts the degrees of freedom corresponding to torsion
    Epsilon(jj) = Ltheta_vec * gamma(orden+1 : 2*orden) + Lbeta * gamma(end);
    
    % 4. REVERSAL AND DIVERGENCE DETECTION
    if jj > 1
        % Reversal: Detect zero-crossing of effectiveness
        if CRnotfound && (Epsilon(jj) * Epsilon(jj-1) < 0)
            % Linear interpolation to find the exact crossing point
            U_CR = interp1([Epsilon(jj-1), Epsilon(jj)], [Uvec2(jj-1), Uvec2(jj)], 0);
            CRnotfound = false;
        end
        
        % Divergence: System matrix becomes singular (determinant changes sign)
        if DIVnotfound 
            if det(TM0_CR) < 0
                U_DIV = U;
                DIVnotfound = false;
            end
        end
    end
end

% 5. DISPLAY RESULTS
disp('--- Critical Aeroelastic Results ---');
if ~CRnotfound
    fprintf('Control Reversal Speed (U_CR): %.2f m/s\n', U_CR);
else
    disp('Control Reversal: NOT FOUND in the analyzed range.');
end

if ~DIVnotfound
    fprintf('Divergence Speed (U_DIV): %.2f m/s\n', U_DIV);
else
    disp('Divergence: NOT FOUND in the analyzed range.');
end