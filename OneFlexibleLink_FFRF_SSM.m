clear; close all; clc;

%% 1. Generate FFRF Model and Tensors
[B, A, fnl, Fext] = build_single_link_ffrf();

%% 2. Setup Dynamical System for SSMTool
DS = DynamicalSystem();
set(DS, 'B', B, 'A', A, 'fnl', fnl);
set(DS.Options, 'Emax', 5, 'Nmax', 11, 'notation', 'multiindex');

epsilon = 0.4; 
kappas = [-1; 1];
coeffs = [Fext, Fext] / 2;
DS.add_forcing(coeffs, kappas, epsilon);

%% 3. Linear Modal Analysis
disp('--> Performing Linear Spectral Analysis...');
[V, D, W] = DS.linear_spectral_analysis();
omega_n = abs(imag(diag(D)));
disp('Natural Frequencies (rad/s):');
disp(omega_n(1:2:end)); 

%% 4. Choose Master Subspace
S = SSM(DS);
set(S.Options, 'reltol', 0.5, 'notation', 'multiindex', 'contribNonAuto', true);
master_modes = [1, 2]; % Mode 1 (Rigid swing) and its complex conjugate
S.choose_E(master_modes);
omega0 = imag(DS.spectrum.Lambda(1)); 
freqRange = omega0 * [0.75, 1.25];
%% 5. Extract Backbone Curve (SSM)
disp('--> Extracting Backbone Curve (SSM)...');
set(S.FRCOptions, 'outdof', 1); 
BB_all = S.extract_backbone(master_modes, freqRange, [3,5,7], 3);

%% 5. Continuation via COCO
disp('--> Starting COCO Continuation for FOM validation...');
try
    omega_start = freqRange(1);
    T = 2*pi/omega_start;
    psp = [epsilon; omega_start];
    
    m = 1.0; L = 1.0; g = 9.81; EI = 15.0; c_theta = 0.1; c_eta = 0.05;

    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [~, x0_trans] = ode15s(@(t,x) ffrf_ode(t, x, psp, m, L, g, EI, c_theta, c_eta), [0 150*T], zeros(4,1), options);
    [ttor, xtor]  = ode15s(@(t,x) ffrf_ode(t, x, psp, m, L, g, EI, c_theta, c_eta), [0 T], x0_trans(end,:)', options);

    prob = coco_prob();
    prob = coco_set(prob, 'ode', 'autonomous', false);
    prob = coco_set(prob, 'cont', 'NAdapt', 2, 'h_max', 0.1, 'PtMX', 1000); 
    
    funcs = {@(t,x,p) ffrf_ode(t, x, p, m, L, g, EI, c_theta, c_eta)};
    coll_args = [funcs, {ttor, xtor, {'eps' 'Omega'}, psp}];
    

    prob = ode_isol2po(prob, '', coll_args{:});

    [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
    maps = data.coll_seg.maps;
    prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, [], 'zero', 'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(2))]);
    prob = coco_add_func(prob, 'amp1', @amplitude, struct('dof',1,'zdim',4), 'regular', 'x1', 'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);

    cont_args = {1, {'Omega' 'po.period' 'x1' 'eps'}, freqRange};
    bd_coco = coco(prob, 'fom_freq_resp', [], cont_args{:});

    om_coco = coco_bd_col(bd_coco, 'Omega');
    x1_coco = coco_bd_col(bd_coco, 'x1');
    eig_coco = coco_bd_col(bd_coco, 'eigs');
    stab_coco = all(abs(eig_coco) < 1 + 1e-5, 1);
    coco_success = true;
    disp('COCO Continuation Successful.');
catch ME
    warning('COCO Continuation failed.');
    disp(ME.message);
    coco_success = false;
end

%% 6. Comparison of COCO and SSM FRC
disp('--> Extracting SSM FRCs...');
figure('Name', 'FRC Comparison', 'Position', [100, 100, 800, 600], 'Color', 'w');

% set up FRC options 
set(S.FRCOptions, 'nCycle', 500, 'initialSolver', 'fsolve');
set(S.contOptions, 'PtMX', 300, 'h_max', 0.1);

% FIXED: Removed the invalid property injection.
set(S.FRCOptions, 'omegaSampStyle', 'cocoBD');
set(S.FRCOptions, 'outdof', 1, 'method', 'continuation ep', 'p0', []);
S.extract_FRC('freq', freqRange, [3, 5, 7, 11]);

if coco_success
    hold on; 
    skip = 20; 
    om_c = om_coco(1:skip:end); x1_c = x1_coco(1:skip:end); st_c = stab_coco(1:skip:end);
    
    plot(om_c(st_c), x1_c(st_c), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.0, 'DisplayName', 'FOM (COCO) Stable');
    plot(om_c(~st_c), x1_c(~st_c), 'kx', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'FOM (COCO) Unstable');
    
    legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12);
    
    xlabel('Forcing Frequency $\Omega$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('Amplitude $\|x_1\|_{\infty}$', 'Interpreter', 'latex', 'FontSize', 16);
    title('FRC with COCO (Softening Effect)', 'Interpreter', 'latex', 'FontSize', 16);
    grid on; set(gca, 'Box', 'on');
end

%% 7. 3D Surface of the SSM & Trajectories Overlay
disp('--> Computing 3D Surface of the SSM and Trajectories...');
highest_order = 5; % Using O(5) for the surface geometry
[W0, R0] = S.compute_whisker(highest_order);

% Determine a starting condition deep inside the nonlinear regime
rho_start = 1.0;          % Initial amplitude 
p_init = rho_start * exp(1i * pi/4); % Initial reduced phase
q_init = [p_init; conj(p_init)];     % Initial reduced state
z0_SSM = reduced_to_full_traj(0, q_init, W0); % Map back to 4D physical state

% --- A. Simulate Trajectory ON the Reduced Manifold (SSM) ---
tf = 20; nsteps = 1500;
traj_ssm = transient_traj_on_auto_ssm(DS, master_modes, W0, R0, tf, nsteps, 1:4, [], q_init);

% --- B. Simulate Exact Full Order Model (FOM) Trajectory ---
[tInt, zInt] = ode15s(@(t,x) ffrf_ode(t, x, [0;0], m, L, g, EI, c_theta, c_eta), [0 tf], z0_SSM, options);

% --- Plot the 3D Surface ---
rhosamp = 0:0.05:1.2; 
plotdofs = [1, 2, 3]; % 1: Theta, 2: Eta, 3: dTheta
figure('Name', '3D SSM Surface and Trajectories', 'Position', [150, 150, 900, 700], 'Color', 'w');

% Plot the transparent bowl (SSM surface)
plot_2D_auto_SSM(W0, rhosamp, plotdofs, {'$\theta$ (Rigid Angle)', '$\eta$ (Flexible Deflection)', '$\dot{\theta}$ (Rigid Vel)'});
hold on;

% --- Overlay the Trajectories ---
% We plot the physical coordinates matching plotdofs: (Theta, Eta, dTheta)
h1 = plot3(traj_ssm.phy(:,1), traj_ssm.phy(:,2), traj_ssm.phy(:,3), '-', 'Color', '#D95319', 'LineWidth', 2.5);
h2 = plot3(zInt(:,1), zInt(:,2), zInt(:,3), 'k--', 'LineWidth', 1.5);

% Formatting
title(['2D Spectral Submanifold $\mathcal{O}(', num2str(highest_order), ')$ with Trajectories'], 'Interpreter', 'latex', 'FontSize', 16);
legend([h1, h2], {'$\mathcal{O}(5)$ SSM Trajectory', 'Exact FOM Trajectory'}, 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
set(gca, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', 12);
view([-35, 30]); 
grid on;


%% 8.  Nonlinear Visualization (Subtracting the Linear Plane)
disp('--> Computing Pure Nonlinear Manifold Curvature (The "Bowl")...');

% Generate custom grid for the surface extraction
phi_grid = linspace(0, 2*pi, 60);
[R_grid, PHI_grid] = meshgrid(rhosamp, phi_grid);
X_surf = zeros(size(R_grid)); Y_surf = zeros(size(R_grid)); 
Z_surf = zeros(size(R_grid)); Z_linear = zeros(size(R_grid));

for i = 1:size(R_grid,1)
    for j = 1:size(R_grid,2)
        p_val = R_grid(i,j) * exp(1i * PHI_grid(i,j));
        
        % Full manifold mapping
        z_val = reduced_to_full_traj(0, [p_val; conj(p_val)], W0);
        X_surf(i,j) = z_val(1);             % Theta (Index 1)
        Y_surf(i,j) = z_val(3);             % dTheta (Index 3)
        Z_surf(i,j) = z_val(2) * 1000;      % Eta (Index 2, Magnified by 1000)
        
        % Linear mapping (to extract the tilted plane)
        z_lin = real(V(:, master_modes) * [p_val; conj(p_val)]); 
        Z_linear(i,j) = z_lin(2) * 1000;    % Eta linear part (Magnified by 1000)
    end
end

Z_bowl = Z_surf - Z_linear;

% Isolate nonlinear curvature for the SSM Trajectory
traj_z_lin = real(V(2, master_modes) * traj_ssm.red.'); 
traj_bowl = traj_ssm.phy(:,2)*1000 - traj_z_lin'*1000;

T_proj = V([1, 3], master_modes);
p_fom = T_proj \ [zInt(:,1)'; zInt(:,3)']; % 2xN matrix of equivalent linear modal coords
fom_z_lin = real(V(2, master_modes) * p_fom)'; % Linear part of Eta for the FOM
fom_bowl = zInt(:,2)*1000 - fom_z_lin*1000;    %  nonlinear deflection for FOM

% Plotting & Visualization (Clean Academic Style)
fig4 = figure('Name', 'Pure Nonlinear Bowl', 'Position', [250, 250, 900, 700], 'Color', 'w');

% Plot the Surface (Sleek Steel Blue)
surf(X_surf, Y_surf, Z_bowl, 'FaceAlpha', 0.85, 'FaceColor', '#5B9BD5', ...
     'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.1, 'DisplayName', 'Pure $\mathcal{O}(2)$ Curvature');
hold on; camlight left; lighting gouraud;

% Overlay the SSM Trajectory (Crisp Orange)
h1_bowl = plot3(traj_ssm.phy(:,1), traj_ssm.phy(:,3), traj_bowl, '-', ...
           'Color', '#D95319', 'LineWidth', 2.5);

% Overlay the EXACT FOM Trajectory (Crisp Dashed Black)
h2_bowl = plot3(zInt(:,1), zInt(:,3), fom_bowl, 'k--', ...
           'LineWidth', 1.5);

% Formatting and Styling
xlabel('Rigid Angle $\theta$ [rad]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Rigid Vel $\dot{\theta}$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 15);
zlabel('Pure Nonlinear Deflection $\Delta\eta \times 10^3$', 'Interpreter', 'latex', 'FontSize', 15);

% Add Custom Legend explicitly assigning the handles
legend([h1_bowl, h2_bowl], {'$\mathcal{O}(7)$ SSM Trajectory', 'Exact FOM Trajectory'}, ...
       'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best', 'EdgeColor', 'k');

view([-40, 25]); 
grid on; 
set(gca, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'GridAlpha', 0.2);


%% ========================================================================
%% HELPER FUNCTIONS FOR COCO 
%% ========================================================================

function dz = ffrf_ode(t, z, p, m, L, g, EI, c_theta, c_eta)
    epsilon = p(1,:); 
    Omega = p(2,:);
    
    N = size(z, 2);
    dz = zeros(4, N);
    
    rhoA = m/L;
    int_phi2 = L/5; int_xphi = L^2/4; int_phi = L/3; int_ddphi2 = 4/L^3;
    C_lin = diag([c_theta, c_eta]);
    
    for i = 1:N
        theta = z(1,i); eta = z(2,i); dtheta = z(3,i); deta = z(4,i);
        t_val = t(min(i, end));
        ep = epsilon(min(i, end)); om = Omega(min(i, end));
        
        M = [ rhoA*(L^3/3 + eta^2 * int_phi2),  rhoA * int_xphi;
              rhoA * int_xphi,                  rhoA * int_phi2 ];

        dM_dt = [rhoA * int_phi2 * 2 * eta * deta, 0; 
                 0,                                0];
        dT_dq = [0; 0.5 * dtheta^2 * (rhoA * int_phi2 * 2 * eta)];

        dV_dq = [ rhoA * g * ( (L^2/2)*sin(theta) + eta * int_phi * cos(theta) );
                  EI * int_ddphi2 * eta + rhoA * g * int_phi * sin(theta) ];

        f_cg = dM_dt * [dtheta; deta] - dT_dq + dV_dq;

        tau_ext = [ep * cos(om * t_val); 0];

        ddq = M \ (tau_ext - C_lin * [dtheta; deta] - f_cg);
        dz(:, i) = [dtheta; deta; ddq(1); ddq(2)];
    end
end

function [data, y] = OmegaT(prob, data, u) 
    y = u(1)*u(2) - 2*pi; 
end

function [data, J] = OmegaT_du(prob, data, u)
    J = [u(2), u(1)];
end

function [data, y] = amplitude(prob, data, u)
    xbp = reshape(u, data.zdim, []);
    y = max(abs(xbp(data.dof, :)), [], 2);
end

function [prob, status, xtr] = amplitude_remesh(prob, data, chart, old_u, old_V)
    [colldata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
    maps = colldata.coll_seg.maps;
    xtr    = [];
    prob   = coco_change_func(prob, data, 'uidx', uidx(maps.xbp_idx));
    status = 'success';
end

function [B, A, fnl, Fext] = build_single_link_ffrf()
    % BUILD_SINGLE_LINK_FFRF: Creates the FFRF model of a single flexible link
    % and extracts the A matrix and F2/F3 tensors for SSMTool.
    
    disp('--> Formulating Symbolic FFRF Model...');
    syms theta eta dtheta deta real
    
    % --- 1. System Parameters ---
    m = 1.0;       % Mass (kg)
    L = 1.0;       % Length (m)
    g = 9.81;      % Gravity (m/s^2)
    EI = 15.0;     % Flexural rigidity (N.m^2)
    c_theta = 0.1; % Joint damping
    c_eta = 0.05;  % Structural damping
    
    rhoA = m / L; 
    
    % Shape function for a clamped-free beam 
    % phi(x) = (x/L)^2
    int_phi2 = L / 5;       % int_0^L phi^2 dx
    int_xphi = L^2 / 4;     % int_0^L x*phi dx
    int_phi  = L / 3;       % int_0^L phi dx
    int_ddphi2 = 4 / L^3;   % int_0^L (d^2phi/dx^2)^2 dx
    
    % --- 2. Energies and Mass Matrix ---
    M = [ rhoA*(L^3/3 + eta^2 * int_phi2),  rhoA * int_xphi;
          rhoA * int_xphi,                  rhoA * int_phi2 ];
      
    % Potential Energy (Elastic + Gravity)
    V_el = 0.5 * EI * int_ddphi2 * eta^2;
    V_g  = rhoA * g * ( -(L^2/2)*cos(theta) + eta * int_phi * sin(theta) );
    V = V_el + V_g;
    
    % --- 3. Equations of Motion ---
    q  = [theta; eta];
    dq = [dtheta; deta];
    z  = [q; dq];
    
    dM_dt = [rhoA * int_phi2 * 2 * eta * deta, 0; 
             0,                                0];
    dT_dq = [0; 0.5 * dtheta^2 * (rhoA * int_phi2 * 2 * eta)];
    dV_dq = jacobian(V, q)';
    
    % Nonlinear Force vector components
    f_cg = dM_dt * dq - dT_dq;
    C_lin = diag([c_theta, c_eta]);
    
    % Total RHS force vector h(z) 
    h_sym = -C_lin * dq - f_cg - dV_dq;
    
    
    disp('--> Performing Implicit Tensor Extraction...');
    n_z = 4;
    z0_val = zeros(n_z, 1);
    
    % 1. Extract symbolic derivatives of M and h
    D1_M = zeros(2,2,n_z); D2_M = zeros(2,2,n_z,n_z);
    D1_h = zeros(2,n_z);   D2_h = zeros(2,n_z,n_z);   D3_h = zeros(2,n_z,n_z,n_z);
    
    for i = 1:n_z
        dM_i = diff(M, z(i));
        D1_M(:,:,i) = double(subs(dM_i, z, z0_val));
        
        dh_i = diff(h_sym, z(i));
        D1_h(:,i) = double(subs(dh_i, z, z0_val));
        
        for j = 1:n_z
            D2_M(:,:,i,j) = double(subs(diff(dM_i, z(j)), z, z0_val));
            dh_ij = diff(dh_i, z(j));
            D2_h(:,i,j) = double(subs(dh_ij, z, z0_val));
            
            for k = 1:n_z
                D3_h(:,i,j,k) = double(subs(diff(dh_ij, z(k)), z, z0_val));
            end
        end
    end
    
    % 2. Assemble Exact Acceleration Tensors Algebraically
    M0 = double(subs(M, z, z0_val));
    Ja = M0 \ D1_h; % Acceleration Jacobian (2x4)
    
    Ha = zeros(2,n_z,n_z); % Acceleration Hessian (2x4x4)
    subs2 = []; vals2 = [];
    
    for i = 1:n_z
        for j = 1:n_z
            term = D2_h(:,i,j) - D1_M(:,:,i)*Ja(:,j) - D1_M(:,:,j)*Ja(:,i);
            Ha(:,i,j) = M0 \ term;
            
            % F2 coefficients (scaled by 1/2!)
            val = Ha(:,i,j) / 2; 
            for r = 1:2
                if abs(val(r)) > 1e-10
                    % Offset row by 2 because accelerations are states 3 and 4
                    subs2 = [subs2; r+2, i, j]; 
                    vals2 = [vals2; val(r)];
                end
            end
        end
    end
    
    Ta = zeros(2,n_z,n_z,n_z); % Acceleration 3rd Derivative (2x4x4x4)
    subs3 = []; vals3 = [];
    
    for i = 1:n_z
        for j = 1:n_z
            for k = 1:n_z
                term = D3_h(:,i,j,k) ...
                       - D1_M(:,:,i)*Ha(:,j,k) - D1_M(:,:,j)*Ha(:,k,i) - D1_M(:,:,k)*Ha(:,i,j) ...
                       - D2_M(:,:,i,j)*Ja(:,k) - D2_M(:,:,j,k)*Ja(:,i) - D2_M(:,:,k,i)*Ja(:,j);
                
                Ta(:,i,j,k) = M0 \ term;
                
                % F3 coefficients (scaled by 1/6!)
                val = Ta(:,i,j,k) / 6; 
                for r = 1:2
                    if abs(val(r)) > 1e-10
                        subs3 = [subs3; r+2, i, j, k]; 
                        vals3 = [vals3; val(r)];
                    end
                end
            end
        end
    end
    
    A = [zeros(2,2), eye(2); Ja]; 
    B = eye(n_z);
    
    F2 = sptensor(subs2, vals2, [n_z, n_z, n_z]);
    if isempty(subs3)
        F3 = sptensor([], [], [n_z, n_z, n_z, n_z]);
    else
        F3 = sptensor(subs3, vals3, [n_z, n_z, n_z, n_z]);
    end
    fnl = {F2, F3};
    
    Fext = [0; 0; M0 \ [1; 0]]; 
    
    disp('--> Single Link FFRF Model Build Complete!');
end