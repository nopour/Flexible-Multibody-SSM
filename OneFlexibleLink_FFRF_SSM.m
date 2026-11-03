%% OneFlexibleLink_FFRF_SSM_R2.m
% Single flexible link (FFRF) - SSM reduction
%
% 
%
% Requirements: SSMTool (+ COCO, Tensor Toolbox), Symbolic Math Toolbox.

clear; close all; clc;

%% 0. User settings ======================================================
P = struct('m',1.0,'L',1.0,'g',9.81,'EI',15.0,'c_theta',0.1,'c_eta',0.05);

epsilon        = 0.4;          % forcing amplitude (torque at the joint)
ORDERS         = [3 5 7 11];   % SSM expansion orders
outdof         = 1;            % output DOF: theta
runPO          = false;        % (C) 'continuation po' cross-check 
stateDepForcing = false;       % true: use exact forcing eps*M(q)^{-1}[tau;0]
                               % (x-dependent, Taylor-expanded) instead of
                               % eps*M0^{-1}[tau;0]. See Section 2.
NF             = 4;            % Taylor order of state-dependent forcing
saveFigs       = true;         % export PNG figures
saveData       = true;         % save FRC data to .mat

%% 1. Generate FFRF Model and Tensors (implicit tensor extraction) =======
[B, A, fnl, Fext, FextData] = build_single_link_ffrf(P, NF);

%% 2. Setup Dynamical System for SSMTool =================================
DS = DynamicalSystem();
set(DS, 'B', B, 'A', A, 'fnl', fnl);
set(DS.Options, 'Emax', 5, 'Nmax', 11, 'notation', 'multiindex');

if stateDepForcing
    % Exact joint torque mapped through the configuration-dependent mass
    
    FextData.epsilon = epsilon;
    DS.add_forcing(FextData);
else
    
    kappas = [-1; 1];
    coeffs = [Fext, Fext] / 2;
    DS.add_forcing(coeffs, kappas, epsilon);
end

%% 3. Linear Modal Analysis ==============================================
disp('--> Performing Linear Spectral Analysis...');
[V, D, W] = DS.linear_spectral_analysis();
omega_n = abs(imag(diag(D)));
disp('Natural Frequencies (rad/s):');
disp(omega_n(1:2:end));

%% 4. Master Subspace and Backbone =======================================
S = SSM(DS);
set(S.Options, 'reltol', 0.5, 'notation', 'multiindex', 'contribNonAuto', true);
master_modes = [1, 2];                 % rigid swing mode and its c.c.
S.choose_E(master_modes);
omega0    = imag(DS.spectrum.Lambda(1));
freqRange = omega0 * [0.75, 1.25];

disp('--> Extracting Backbone Curve (SSM)...');
set(S.FRCOptions, 'outdof', outdof);
BB_all = S.extract_backbone(master_modes, freqRange, [3,5,7], 3);

%% 5. FOM: periodic orbits via COCO ======================================
disp('--> Starting COCO Continuation for FOM validation...');
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
try
    omega_start = freqRange(1);
    T   = 2*pi/omega_start;
    psp = [epsilon; omega_start];

    [~, x0_trans] = ode15s(@(t,x) ffrf_ode(t, x, psp, P), [0 150*T], zeros(4,1), options);
    [ttor, xtor]  = ode15s(@(t,x) ffrf_ode(t, x, psp, P), [0 T], x0_trans(end,:)', options);

    prob = coco_prob();
    prob = coco_set(prob, 'ode', 'autonomous', false);
    prob = coco_set(prob, 'cont', 'NAdapt', 2, 'h_max', 0.1, 'PtMX', 1000);

    funcs     = {@(t,x,p) ffrf_ode(t, x, p, P)};
    coll_args = [funcs, {ttor, xtor, {'eps' 'Omega'}, psp}];
    prob = ode_isol2po(prob, '', coll_args{:});

    [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
    maps = data.coll_seg.maps;
    prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, [], 'zero', ...
        'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(2))]);
    prob = coco_add_func(prob, 'amp1', @amplitude, struct('dof',outdof,'zdim',4), ...
        'regular', 'x1', 'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);

    cont_args = {1, {'Omega' 'po.period' 'x1' 'eps'}, freqRange};
    bd_coco = coco(prob, 'fom_freq_resp', [], cont_args{:});

    om_coco   = coco_bd_col(bd_coco, 'Omega');
    x1_coco   = coco_bd_col(bd_coco, 'x1');
    eig_coco  = coco_bd_col(bd_coco, 'eigs');
    stab_coco = all(abs(eig_coco) < 1 + 1e-5, 1);
    coco_success = true;
    disp('COCO Continuation Successful.');
catch ME
    warning('COCO Continuation failed.');
    disp(ME.message);
    coco_success = false;
end

%% 6. SSM-based FRCs: leading-order vs higher-order non-autonomous =======
set(S.FRCOptions, 'nCycle', 500, 'initialSolver', 'fsolve', 'outdof', outdof);
set(S.contOptions, 'PtMX', 300, 'h_max', 0.1);
coord0 = S.FRCOptions.coordinates;

FRC_ep = cell(numel(ORDERS),1);   % (A) leading-order non-autonomous
FRC_ls = cell(numel(ORDERS),1);   % (B) higher-order non-autonomous
FRC_po = cell(numel(ORDERS),1);   % (C) higher-order non-autonomous (optional)
tCPU   = zeros(numel(ORDERS),3);

% ---- (A) continuation ep: leading-order non-autonomous truncation -----
disp('--> (A) continuation ep  [leading-order non-autonomous SSM]');
for j = 1:numel(ORDERS)
    set(S.FRCOptions, 'method', 'continuation ep', 'sampStyle', 'cocoBD', ...
        'coordinates', coord0, 'p0', [], 'z0', []);
    t0 = tic;
    FRC_ep{j} = run_frc(S, freqRange, ORDERS(j));
    tCPU(j,1) = toc(t0);
end

% ---- (B) level set: non-autonomous terms up to O(order-1) --------------
% omDepNonAuto = true recomputes the non-autonomous SSM coefficients

disp('--> (B) level set        [higher-order non-autonomous SSM]');
set(S.FRCOptions, 'method', 'level set', 'omDepNonAuto', true, ...
    'nPar', 200, 'nRho', 400, 'nPsi', 200, 'rhoScale', 2);
for j = 1:numel(ORDERS)
    t0 = tic;
    FRC_ls{j} = run_frc(S, freqRange, ORDERS(j));
    tCPU(j,2) = toc(t0);
end

% ---- (C) continuation po: non-autonomous terms up to O(order-1) --------


if runPO
    disp('--> (C) continuation po  [higher-order non-autonomous SSM]');
    for j = 1:numel(ORDERS)
        set(S.FRCOptions, 'method', 'continuation po', 'coordinates', 'cartesian', ...
            'omDepNonAuto', false, 'omDepNonAutoVal', omega0, 'p0', [], 'z0', []);
        t0 = tic;
        FRC_po{j} = run_frc(S, freqRange, ORDERS(j));
        tCPU(j,3) = toc(t0);
    end
    set(S.FRCOptions, 'coordinates', coord0);
end

%% 7. Quantitative comparison at the resonance peak ======================
fprintf('\n================ Resonance-peak comparison (eps = %.3g) ================\n', epsilon);
if coco_success
    [Apk_fom, ik] = max(x1_coco);  Wpk_fom = om_coco(ik);
    fprintf('FOM (COCO)                     : Omega_pk = %.5f, A_pk = %.6f\n', Wpk_fom, Apk_fom);
else
    Apk_fom = NaN; Wpk_fom = NaN;
end
fprintf('%-6s %-22s %-11s %-11s %-12s %-10s\n','order','method','Omega_pk','A_pk','relErr_A[%]','CPU[s]');
names = {'ep (leading-order)','level set (higher)','po (higher)'};
allF  = {FRC_ep, FRC_ls, FRC_po};
peakTable = nan(numel(ORDERS), 3, 3);   % [Omega_pk, A_pk, relErr]
for j = 1:numel(ORDERS)
    for k = 1:3
        if isempty(allF{k}{j}), continue; end
        [om, Aa] = frc_xy(allF{k}{j});
        [Apk, ii] = max(Aa);
        err = 100*abs(Apk - Apk_fom)/Apk_fom;
        peakTable(j,k,:) = [om(ii), Apk, err];
        fprintf('%-6d %-22s %-11.5f %-11.6f %-12.4f %-10.2f\n', ORDERS(j), names{k}, om(ii), Apk, err, tCPU(j,k));
    end
end
fprintf('=========================================================================\n\n');

%% 8. FRC comparison figure ==============================================
jHi  = numel(ORDERS);                  % highest order
cols = lines(numel(ORDERS));
figFRC = figure('Name','FRC: non-autonomous truncation','Color','w','Position',[80 80 1250 520]);
tl = tiledlayout(figFRC,1,2,'TileSpacing','compact','Padding','compact');

for tile = 1:2
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    h = gobjects(0); lab = {};
    % leading-order (ep) for all orders, lines (solid stable / dashed unstable)
    for j = 1:numel(ORDERS)
        hj = plot_stab_line(ax, FRC_ep{j}, cols(j,:), 1.2);
        h(end+1) = hj; lab{end+1} = sprintf('ep, $\\mathcal{O}(%d)$, leading-order $\\mathbf{X}_0(\\phi)$', ORDERS(j)); %#ok<SAGROW>
    end
    % higher-order (level set) at highest order, circles
    [om, Aa, st] = frc_xy(FRC_ls{jHi});
    h(end+1) = plot(ax, om(st),  Aa(st),  'o', 'Color',[0.85 0.33 0.10], 'MarkerSize',5, 'LineWidth',1.0);
    lab{end+1} = sprintf('level set, $\\mathcal{O}(%d)$, $\\mathbf{X}(\\mathbf{p},\\phi)$ to $\\mathcal{O}(%d)$', ORDERS(jHi), ORDERS(jHi)-1);
    if any(~st)
        h(end+1) = plot(ax, om(~st), Aa(~st), 's', 'Color',[0.85 0.33 0.10], 'MarkerSize',5);
        lab{end+1} = 'level set, unstable';
    end
    % optional po
    if runPO && ~isempty(FRC_po{jHi})
        [om, Aa] = frc_xy(FRC_po{jHi});
        h(end+1) = plot(ax, om, Aa, 'd', 'Color',[0.47 0.67 0.19], 'MarkerSize',5);
        lab{end+1} = sprintf('po, $\\mathcal{O}(%d)$', ORDERS(jHi));
    end
    % FOM
    if coco_success
        skip = 10;
        oc = om_coco(1:skip:end); xc = x1_coco(1:skip:end); sc = stab_coco(1:skip:end);
        h(end+1) = plot(ax, oc(sc),  xc(sc),  'ko', 'MarkerFaceColor','w', 'MarkerSize',6, 'LineWidth',1.0);
        lab{end+1} = 'FOM (COCO), stable';
        if any(~sc)
            h(end+1) = plot(ax, oc(~sc), xc(~sc), 'kx', 'MarkerSize',7, 'LineWidth',1.5);
            lab{end+1} = 'FOM (COCO), unstable';
        end
    end
    xlabel(ax,'Forcing frequency $\Omega$ [rad/s]','Interpreter','latex','FontSize',14);
    ylabel(ax,'$\|\theta\|_{\infty}$ [rad]','Interpreter','latex','FontSize',14);
    set(ax,'TickLabelInterpreter','latex','FontSize',12);
    if tile == 1
        xlim(ax, freqRange);
        title(ax, 'Full frequency range', 'Interpreter','latex','FontSize',14);
        axLeg = ax; hLeg = h; labLeg = lab;
    else
        
        pk = squeeze(peakTable(jHi,2,:));
        if all(isnan(pk)), pk = squeeze(peakTable(jHi,1,:)); end
        if coco_success, Wc = Wpk_fom; Ac = Apk_fom; else, Wc = pk(1); Ac = pk(2); end
        xlim(ax, Wc + [-0.05 0.13]);
        ylim(ax, Ac*[0.90 1.03]);
        title(ax, 'Resonance peak', 'Interpreter','latex','FontSize',14);
    end
end

lgd = legend(axLeg, hLeg, labLeg, 'Interpreter','latex','FontSize',10, ...
    'NumColumns',3,'Box','off');
lgd.Layout.Tile = 'south';
if saveFigs
    exportgraphics(figFRC, 'FRC_nonauto_alink.pdf', 'ContentType', 'vector');   % Fig. 
end

%% 9. Effect of the non-autonomous truncation alone ======================
% Same autonomous order, leading-order vs higher-order non-autonomous part.
figDiff = figure('Name','Leading vs higher-order non-autonomous','Color','w','Position',[120 120 650 480]);
ax = axes(figDiff); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
for j = 1:numel(ORDERS)
    plot_stab_line(ax, FRC_ep{j}, cols(j,:), 1.2);
    [om, Aa] = frc_xy(FRC_ls{j});
    plot(ax, om, Aa, 'o', 'Color', cols(j,:), 'MarkerSize', 4, 'HandleVisibility','off');
end
if coco_success
    plot(ax, om_coco(1:10:end), x1_coco(1:10:end), 'k.', 'MarkerSize', 10, 'DisplayName','FOM (COCO)');
end
legend(ax, [arrayfun(@(o) sprintf('$\\mathcal{O}(%d)$: lines = ep (leading), circles = level set', o), ORDERS, 'UniformOutput', false), ...
    repmat({'FOM (COCO)'}, 1, coco_success)], 'Interpreter','latex','FontSize',9,'Location','northeast');
xlabel(ax,'$\Omega$ [rad/s]','Interpreter','latex','FontSize',14);
ylabel(ax,'$\|\theta\|_{\infty}$ [rad]','Interpreter','latex','FontSize',14);
set(ax,'TickLabelInterpreter','latex','FontSize',12);
xlim(ax, freqRange);
if saveFigs, exportgraphics(figDiff, 'FRC_ep_vs_levelset_alink.pdf', 'ContentType', 'vector'); end   % not used in the paper

if saveData
    save(sprintf('FRC_nonauto_eps%.2g.mat', epsilon), 'ORDERS', 'epsilon', 'freqRange', ...
        'FRC_ep', 'FRC_ls', 'FRC_po', 'peakTable', 'tCPU', 'stateDepForcing');
    if coco_success
        save(sprintf('FRC_nonauto_eps%.2g.mat', epsilon), 'om_coco', 'x1_coco', 'stab_coco', '-append');
    end
end

%% 10. 3D Surface of the autonomous SSM & Trajectories ===================
disp('--> Computing 3D Surface of the SSM and Trajectories...');
highest_order = 5;
[W0, R0] = S.compute_whisker(highest_order);

rho_start = 1.0;
p_init = rho_start * exp(1i * pi/4);
q_init = [p_init; conj(p_init)];
z0_SSM = reduced_to_full_traj(0, q_init, W0);

tf = 20; nsteps = 1500;
traj_ssm = transient_traj_on_auto_ssm(DS, master_modes, W0, R0, tf, nsteps, 1:4, [], q_init);
[tInt, zInt] = ode15s(@(t,x) ffrf_ode(t, x, [0;0], P), [0 tf], z0_SSM, options);

rhosamp  = 0:0.05:1.2;
plotdofs = [1, 2, 3];
plot_2D_auto_SSM(W0, rhosamp, plotdofs, {'$\theta$ (Rigid Angle)', '$\eta$ (Flexible Deflection)', '$\dot{\theta}$ (Rigid Vel)'});
fig3 = gcf;            
set(fig3, 'Name', '3D SSM Surface and Trajectories', 'Position', [150, 150, 900, 700], 'Color', 'w');
hold on;
h1 = plot3(traj_ssm.phy(:,1), traj_ssm.phy(:,2), traj_ssm.phy(:,3), '-', 'Color', '#D95319', 'LineWidth', 2.5);
h2 = plot3(zInt(:,1), zInt(:,2), zInt(:,3), 'k--', 'LineWidth', 1.5);
title(['2D Spectral Submanifold $\mathcal{O}(', num2str(highest_order), ')$ with Trajectories'], 'Interpreter', 'latex', 'FontSize', 16);
legend([h1, h2], {sprintf('$\\mathcal{O}(%d)$ SSM Trajectory', highest_order), 'Exact FOM Trajectory'}, ...
    'Interpreter', 'latex', 'FontSize', 13, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
set(gca, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', 12);
view([-35, 30]); grid on;
if saveFigs, exportgraphics(fig3, 'AlinkManifolds.pdf', 'ContentType', 'image', 'Resolution', 400); end   % Fig. 3a

%% 11. Nonlinear curvature (linear plane subtracted) =====================
disp('--> Computing Pure Nonlinear Manifold Curvature...');
phi_grid = linspace(0, 2*pi, 60);
[R_grid, PHI_grid] = meshgrid(rhosamp, phi_grid);
X_surf = zeros(size(R_grid)); Y_surf = X_surf; Z_surf = X_surf; Z_linear = X_surf;
for i = 1:size(R_grid,1)
    for j = 1:size(R_grid,2)
        p_val = R_grid(i,j) * exp(1i * PHI_grid(i,j));
        z_val = reduced_to_full_traj(0, [p_val; conj(p_val)], W0);
        X_surf(i,j) = z_val(1);
        Y_surf(i,j) = z_val(3);
        Z_surf(i,j) = z_val(2) * 1000;
        z_lin = real(V(:, master_modes) * [p_val; conj(p_val)]);
        Z_linear(i,j) = z_lin(2) * 1000;
    end
end
Z_bowl = Z_surf - Z_linear;

traj_z_lin = real(V(2, master_modes) * traj_ssm.red.');
traj_bowl  = traj_ssm.phy(:,2)*1000 - traj_z_lin'*1000;
T_proj     = V([1, 3], master_modes);
p_fom      = T_proj \ [zInt(:,1)'; zInt(:,3)'];
fom_z_lin  = real(V(2, master_modes) * p_fom)';
fom_bowl   = zInt(:,2)*1000 - fom_z_lin*1000;

fig4 = figure('Name', 'Pure Nonlinear Bowl', 'Position', [250, 250, 900, 700], 'Color', 'w');
surf(X_surf, Y_surf, Z_bowl, 'FaceAlpha', 0.85, 'FaceColor', '#5B9BD5', ...
     'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.1);
hold on; camlight left; lighting gouraud;
h1_bowl = plot3(traj_ssm.phy(:,1), traj_ssm.phy(:,3), traj_bowl, '-', 'Color', '#D95319', 'LineWidth', 2.5);
h2_bowl = plot3(zInt(:,1), zInt(:,3), fom_bowl, 'k--', 'LineWidth', 1.5);
xlabel('Rigid Angle $\theta$ [rad]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Rigid Vel $\dot{\theta}$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 15);
zlabel('Pure Nonlinear Deflection $\Delta\eta \times 10^3$', 'Interpreter', 'latex', 'FontSize', 15);
legend([h1_bowl, h2_bowl], {sprintf('$\\mathcal{O}(%d)$ SSM Trajectory', highest_order), 'Exact FOM Trajectory'}, ...
       'Interpreter', 'latex', 'FontSize', 13, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
view([-40, 25]); grid on;
set(gca, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'GridAlpha', 0.2);
if saveFigs, exportgraphics(fig4, 'ManifoldBowl_alink.pdf', 'ContentType', 'image', 'Resolution', 400); end   % Fig. 3b


%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function FRC = run_frc(S, freqRange, order)
% Runs S.extract_FRC for ONE order, returns the FRC struct array and

    figsBefore = findall(0, 'Type', 'figure');
    FRC = S.extract_FRC('freq', freqRange, order);
    figsNew = setdiff(findall(0, 'Type', 'figure'), figsBefore);
    close(figsNew);
end

function [om, A, st] = frc_xy(FRC)

    if isfield(FRC, 'Omega'), om = [FRC.Omega]; else, om = [FRC.om]; end
    A  = [FRC.Aout];
    st = logical([FRC.stability]);
    om = om(:); A = A(:); st = st(:);
    if numel(A) ~= numel(om), A = A(1:numel(om)); end
end

function h = plot_stab_line(ax, FRC, col, lw)

    [om, A, st] = frc_xy(FRC);
    As = A; As(~st) = NaN;
    Au = A; Au(st)  = NaN;
   
    chg = find(diff(st) ~= 0);
    Au(chg) = A(chg); Au(chg+1) = A(chg+1);
    h = plot(ax, om, As, '-',  'Color', col, 'LineWidth', lw);
    plot(ax, om, Au, '--', 'Color', col, 'LineWidth', lw, 'HandleVisibility', 'off');
end

function dz = ffrf_ode(t, z, p, P)

    epsilon = p(1,:);
    Omega   = p(2,:);
    m = P.m; L = P.L; g = P.g; EI = P.EI;

    N  = size(z, 2);
    dz = zeros(4, N);
    rhoA = m/L;
    int_phi2 = L/5; int_xphi = L^2/4; int_phi = L/3; int_ddphi2 = 4/L^3;
    C_lin = diag([P.c_theta, P.c_eta]);

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
    maps   = colldata.coll_seg.maps;
    xtr    = [];
    prob   = coco_change_func(prob, data, 'uidx', uidx(maps.xbp_idx));
    status = 'success';
end

function [B, A, fnl, Fext, FextData] = build_single_link_ffrf(P, NF)
% BUILD_SINGLE_LINK_FFRF: FFRF model of a single flexible link; implicit
% extraction of A, F2, F3 (no symbolic inversion of M(q)).
% Also returns FextData: Taylor expansion (up to order NF) of the exact
% forcing direction [0; M(q)^{-1}[1;0]] in SSMTool multi-index format.

    disp('--> Formulating Symbolic FFRF Model...');
    syms theta eta dtheta deta real

    m = P.m; L = P.L; g = P.g; EI = P.EI;
    rhoA = m / L;
    int_phi2 = L / 5; int_xphi = L^2 / 4; int_phi = L / 3; int_ddphi2 = 4 / L^3;

    M = [ rhoA*(L^3/3 + eta^2 * int_phi2),  rhoA * int_xphi;
          rhoA * int_xphi,                  rhoA * int_phi2 ];
    V_el = 0.5 * EI * int_ddphi2 * eta^2;
    V_g  = rhoA * g * ( -(L^2/2)*cos(theta) + eta * int_phi * sin(theta) );
    Vpot = V_el + V_g;

    q  = [theta; eta];
    dq = [dtheta; deta];
    z  = [q; dq];

    dM_dt = [rhoA * int_phi2 * 2 * eta * deta, 0;
             0,                                0];
    dT_dq = [0; 0.5 * dtheta^2 * (rhoA * int_phi2 * 2 * eta)];
    dV_dq = jacobian(Vpot, q)';

    f_cg  = dM_dt * dq - dT_dq;
    C_lin = diag([P.c_theta, P.c_eta]);
    h_sym = -C_lin * dq - f_cg - dV_dq;

    disp('--> Performing Implicit Tensor Extraction...');
    n_z = 4;
    z0_val = zeros(n_z, 1);

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

    M0 = double(subs(M, z, z0_val));
    Ja = M0 \ D1_h;

    Ha = zeros(2,n_z,n_z); subs2 = []; vals2 = [];
    for i = 1:n_z
        for j = 1:n_z
            term = D2_h(:,i,j) - D1_M(:,:,i)*Ja(:,j) - D1_M(:,:,j)*Ja(:,i);
            Ha(:,i,j) = M0 \ term;
            val = Ha(:,i,j) / 2;
            for r = 1:2
                if abs(val(r)) > 1e-10
                    subs2 = [subs2; r+2, i, j]; 
                    vals2 = [vals2; val(r)];    
                end
            end
        end
    end

    subs3 = []; vals3 = [];
    for i = 1:n_z
        for j = 1:n_z
            for k = 1:n_z
                term = D3_h(:,i,j,k) ...
                     - D1_M(:,:,i)*Ha(:,j,k) - D1_M(:,:,j)*Ha(:,k,i) - D1_M(:,:,k)*Ha(:,i,j) ...
                     - D2_M(:,:,i,j)*Ja(:,k) - D2_M(:,:,j,k)*Ja(:,i) - D2_M(:,:,k,i)*Ja(:,j);
                val = (M0 \ term) / 6;
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
    fnl  = {F2, F3};
    Fext = [0; 0; M0 \ [1; 0]];

    % ---- State-dependent forcing direction, Taylor-expanded (multi-index)
    % F(z) = [0;0; M(q)^{-1} e1] = sum_k sum_|m|=k  F_m z^m
    bsym = [0; 0; M \ [1; 0]];
    FextData = struct();
    kap = [-1; 1];
    for s = 1:2
        FextData.data(s).kappa = kap(s);
    end
    for k = 0:NF
        I = multi_indices(n_z, k);          % (nTerms x n_z)
        C = zeros(n_z, size(I,1));
        for r = 1:size(I,1)
            Dr = bsym;
            for v = 1:n_z
                if I(r,v) > 0, Dr = diff(Dr, z(v), I(r,v)); end
            end
            C(:,r) = double(subs(Dr, z, z0_val)) / prod(factorial(I(r,:)));
        end
        keep = any(abs(C) > 1e-12, 1);
        if ~any(keep)            
            keep = false(1,size(I,1)); keep(1) = true; C(:,1) = 0;
        end
        for s = 1:2
            FextData.data(s).F_n_k(k+1).coeffs = C(:,keep) / 2;
            FextData.data(s).F_n_k(k+1).ind    = I(keep,:);
        end
    end

    disp('--> Single Link FFRF Model Build Complete!');
end

function I = multi_indices(n, k)

    if k == 0, I = zeros(1, n); return; end
    c = nchoosek(1:(k+n-1), n-1);
    I = diff([zeros(size(c,1),1), c, (k+n)*ones(size(c,1),1)], 1, 2) - 1;
end
