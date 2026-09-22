function R = apollo_lkf_mission()
%APOLLO_LKF_MISSION  Circumlunar free-return navigation with the first
%   spaceflight Kalman filter, reconstructed from the original sources.
%
%   A learning project, built as an extension of an applied Kalman
%   filtering specialization: the linearized Kalman filter (LKF) that was
%   first proposed for onboard cislunar navigation is implemented as
%   published, and flown over a complete Earth -> Moon -> Earth ballistic
%   free-return trajectory.
%
%   The reconstruction follows, and cites throughout:
%
%   [R135]  G. L. Smith, S. F. Schmidt, L. A. McGee, "Application of
%           Statistical Filter Theory to the Optimal Estimation of
%           Position and Velocity on Board a Circumlunar Vehicle,"
%           NASA TR R-135, 1962.  (Dynamics: App. A; perturbation
%           dynamics: App. B; optical angles and partials: App. C;
%           filter equations: eqns (12)-(16).)
%   [BELL]  D. A. Corey et al., "Summary of Apollo Guidance and
%           Navigation Error Analysis," Bellcomm TR-66-310-4, 1966.
%           (Sextant noise: 10 arc-sec nominal; entry interface at
%           400,000 ft.)
%   [A13]   "Apollo 13 Mission Report, Supplement 1: Guidance, Navigation
%           and Control Systems Performance Analysis," NASA MSC-02680
%           Suppl. 1 (TRW), 1970.  (Sec. 3.2: P23 star-horizon sightings
%           showed an ~18 km actual vs ~10 km apparent horizon altitude;
%           sec. 3.1.4: low-level PIPA outputs first read as venting.)
%   [TM847] L. A. McGee, S. F. Schmidt, "Discovery of the Kalman Filter
%           as a Practical Tool for Aerospace and Industry," NASA
%           TM-86847, 1985.  (Historical account; linearizing about the
%           current estimate -> the extended Kalman filter.)
%
%   Two filter variants are run over the same data:
%     'LKF'  linearized about the PRE-COMPUTED NOMINAL trajectory --
%            the textbook linearized Kalman filter;
%     'EKF'  re-linearized about the CURRENT ESTIMATE, the modification
%            R-135 already recommends on p. 15 and which TM-86847
%            identifies as the birth of the extended Kalman filter.
%
%   Usage (MATLAB, no toolboxes):
%       R = apollo_lkf_mission();     % runs everything, saves figures,
%                                     % writes apollo_lkf_results.mat
%   Runtime is a few minutes (dominated by the trajectory targeting).
%   Figures are saved as PNG in the current folder.  The companion file
%   apollo_lkf_animation.m renders the mission animation from the saved
%   results.

% =========================================================================
% 0. CONFIGURATION  (every number is a recorded design decision)
% =========================================================================
C = constants_();

cfg.dt        = 60;                 % integration/filter grid step [s]
cfg.dt_meas   = 900;                % time between sightings [s].  [R135]
                                    % studied 20 obs at 6-min spacing over
                                    % the first 2.5 h; here the cadence is
                                    % extended over the whole mission.
cfg.t_first   = 1800;               % first observation 1/2 h after
                                    % injection, as in [R135] p. 16.
cfg.sig_ang   = 10*C.arcsec;        % 1-sigma angle noise, all angles:
                                    % sextant nominal 10 arc-sec [BELL].
cfg.seed      = 7;                  % reproducibility

% --- targeting: the free-return "figure-8" --------------------------------
tgt.r_perilune = C.R_M + 111;       % 111-km pericynthion class pass
tgt.r_perigee  = C.R0  + 40;        % vacuum perigee below the 400,000-ft
                                    % entry interface [BELL] -> re-entry
tgt.tol        = 0.5;               % km
% Warm start for the targeting iteration: recorded output of the coarse
% corridor scan (set cfg.do_scan = true to redo the scan from scratch,
% ~2 min: it sweeps v0 x th_v and keeps the best RETROGRADE lunar pass).
cfg.do_scan  = false;
cfg.tli_warm = [10.9649; deg2rad(-132.77)];

% --- truth-model errors (what the filter does not know) -------------------
% Injection error: the sampled TLI error of the R-135 simulation study
% (p. 19), scaled by 0.05.  The full R-135 sample (~2 m/s) produces a
% multi-thousand-km deviation at the Moon (R-135 reports a 4528-km miss
% on their trajectory) -- a real mission trims this with midcourse
% corrections, which this ballistic simulation deliberately omits.  The
% 0.05 scale represents post-midcourse dispersion; its sign is chosen so
% the dispersion raises the lunar pass (true perilune ~244 km vs the
% nominal 111 km) -- the uncorrected trajectory stays safely clear of
% the Moon.  Recorded design decisions.
cfg.dx0 = -0.05*[0.495; -0.886; -1.001; 0.281e-3; 1.999e-3; 0.194e-3];

% Unmodeled constant acceleration on the true vehicle: a venting /
% outgassing-type disturbance of ~0.7 micro-g.  Motivated by [A13]
% sec. 3.1.4 (PIPA outputs after the LOX-tank incident first interpreted
% as venting).  The filter carries it as an augmented state.
cfg.a_vent = [-3; 6; -2]*1e-9;      % [km/s^2]

% Earth-horizon bias: [A13] sec. 3.2 -- P23 processing showed the actual
% horizon altitude was ~18 km while the apparent (sighted) altitude was
% ~10 km.  The ~8 km discrepancy is taken as the true horizon-altitude
% bias on the subtended-Earth-angle measurement.  R-135 (p. 21) found
% gamma_e to be precisely the angle most sensitive to bias, and suggests
% estimating such biases as additional state variables -- done here.
cfg.b_horizon = 8.0;                % [km]

% --- filter statistics -----------------------------------------------------
% P0: injection covariance of [R135] p. 18: 1 km and 1 m/s per axis,
% plus the augmented parameters.
P0 = blkdiag( eye(3)*1^2, ...               % position [km^2]
              eye(3)*(1e-3)^2, ...          % velocity [(km/s)^2]
              eye(3)*(1e-8)^2, ...          % accel    [(km/s^2)^2]
              10^2 );                       % horizon bias [km^2]
% Q: small white acceleration noise covering linearization residue and
% integration error; the venting acceleration itself is a STATE, so Q on
% the acceleration states is a slow random walk.
cfg.q_acc  = (4e-11)^2;             % PSD on velocity states [(km/s^2)^2*s]
cfg.q_adot = (1e-14)^2;             % random walk of a_u [(km/s^3)^2*s]

% =========================================================================
% 1. NOMINAL TRAJECTORY -- target the free return
% =========================================================================
fprintf('[1/5] targeting the ballistic free-return trajectory...\n');
[p_tli, tinfo] = shoot_free_return_(C, tgt, cfg);
fprintf('      TLI: v0 = %.6f km/s at injection angle %.4f deg\n', ...
        p_tli(1), rad2deg(p_tli(2)));
fprintf('      perilune alt %.1f km at t = %.2f h; return perigee alt %.1f km at t = %.2f h\n', ...
        tinfo.detail.r_pl - C.R_M, tinfo.detail.t_pl/C.hr, ...
        tinfo.detail.r_pe - C.R0,  tinfo.detail.t_pe/C.hr);

% mission grid: TLI to entry interface (400,000 ft, [BELL])
X0n = tli_state_(p_tli, C);
t_end  = tinfo.detail.t_pe;                       % provisional
t      = 0:cfg.dt:(t_end + 2*C.hr);
Xnom   = rk4_prop_(t, X0n, C, [0;0;0]);
kEnd   = find(sqrt(sum(Xnom(1:3,:).^2,1)) < C.r_entry & t > tinfo.detail.t_pl, 1);
if isempty(kEnd), kEnd = numel(t); end
t    = t(1:kEnd);   Xnom = Xnom(:,1:kEnd);   N = numel(t);
fprintf('      mission duration TLI -> entry interface: %.2f h (%.2f days)\n', ...
        t(end)/C.hr, t(end)/86400);

% =========================================================================
% 2. TRUE TRAJECTORY -- injection error + unmodeled venting acceleration
% =========================================================================
fprintf('[2/5] propagating the true trajectory...\n');
Xtrue = rk4_prop_(t, X0n + cfg.dx0, C, cfg.a_vent);

rM = ephem_moon_(t, C);
devP = sqrt(sum((Xtrue(1:3,:) - Xnom(1:3,:)).^2, 1));
dMoonTrue = sqrt(sum((Xtrue(1:3,:) - rM).^2, 1));
fprintf('      deviation from nominal: %.2f km at TLI -> %.0f km at perilune -> %.0f km at entry\n', ...
        devP(1), interp1(t, devP, tinfo.detail.t_pl), devP(end));
fprintf('      true closest lunar approach: %.1f km altitude (nominal 111 km)\n', ...
        min(dMoonTrue) - C.R_M);

% =========================================================================
% 3. MEASUREMENTS -- R-135 optical angles, generated from the truth
% =========================================================================
fprintf('[3/5] generating optical sightings...\n');
seed_(cfg.seed);
MEAS = make_meas_(t, Xtrue, C, cfg);
fprintf('      %d sighting epochs, %d scalar angles (%d during Earth occultation dropped)\n', ...
        MEAS.n_epoch, MEAS.n_scalar, MEAS.n_occulted);

% =========================================================================
% 4. FILTERS -- LKF about the nominal; EKF variant about the estimate
% =========================================================================
fprintf('[4/5] running the filters...\n');
FL = run_filter_('LKF', t, Xnom, MEAS, P0, C, cfg);
FE = run_filter_('EKF', t, Xnom, MEAS, P0, C, cfg);

errL = FL.Xest(1:3,:) - Xtrue(1:3,:);
errE = FE.Xest(1:3,:) - Xtrue(1:3,:);
fprintf('      final position error: LKF %.2f km | EKF %.2f km (open-loop %.0f km)\n', ...
        norm(errL(:,end)), norm(errE(:,end)), devP(end));
fprintf('      RMS position error over mission: LKF %.2f km | EKF %.2f km\n', ...
        sqrt(mean(sum(errL.^2,1))), sqrt(mean(sum(errE.^2,1))));

% =========================================================================
% 5. RESULTS -- figures + saved state for the animation
% =========================================================================
fprintf('[5/5] plotting and saving...\n');
R = struct('t',t, 'Xnom',Xnom, 'Xtrue',Xtrue, 'rM',rM, 'C',C, 'cfg',cfg, ...
           'tgt',tgt, 'p_tli',p_tli, 'tinfo',tinfo, 'MEAS',MEAS, ...
           'FL',FL, 'FE',FE);
save('apollo_lkf_results.mat', '-struct', 'R');
make_figures_(R);
fprintf('done.  Results in apollo_lkf_results.mat, figures fig1..fig6*.png\n');
end

% =========================================================================
%  CONSTANTS AND EPHEMERIDES
% =========================================================================
function C = constants_()
% Values used in the original study: NASA TR R-135, App. A ([R135]).
% (R-135's "J" is the coefficient of its second-harmonic term.)
C.mu_E = 3.986135e5;      % km^3/s^2
C.mu_M = 4.89820e3;       % km^3/s^2
C.mu_S = 1.3253e11;       % km^3/s^2
C.a_E  = 6378.26;         % km, Earth equatorial radius (J2 term)
C.J    = 1.6246e-3;       % second-harmonic coefficient
C.R0   = 6371.0;          % km, horizon reference radius (R-135 App. C)
C.R_M  = 1737.4;          % km, Moon radius
C.D_EM = 384400;  C.D_SE = 1.496e8;                 % km
C.om_M = 2*pi/(27.321661*86400);                    % rad/s (sidereal)
C.om_S = 2*pi/(365.25636*86400);                    % rad/s
C.thM0 = 0.0;  C.thS0 = 2.0;                        % phases at t=0
C.SOI_M   = C.D_EM*(C.mu_M/C.mu_E)^(2/5);           % lunar SOI, km
C.r_park  = C.R0 + 185;                             % 185-km parking orbit
C.r_entry = C.R0 + 121.9;   % entry interface, 400,000 ft [BELL glossary]
C.arcsec  = pi/(180*3600);
C.hr      = 3600;
end

function [rM, vM] = ephem_moon_(tv, C)
% Circular kinematic lunar ephemeris in the XY plane (design decision:
% R-135 used a stored ephemeris; the circular approximation is ours).
th = C.thM0 + C.om_M*tv(:).';
rM = C.D_EM*[cos(th); sin(th); zeros(1,numel(tv))];
if nargout > 1
    vM = C.D_EM*C.om_M*[-sin(th); cos(th); zeros(1,numel(tv))];
end
end

function X0 = tli_state_(p, C)
% Tangential prograde injection from the parking orbit at angle p(2).
v0 = p(1); thv = p(2);
X0 = [C.r_park*[cos(thv); sin(thv); 0];
      v0*     [-sin(thv); cos(thv); 0]];
end

% =========================================================================
%  DYNAMICS -- R-135 Appendix A, eqns (A1)-(A3)
% =========================================================================
function X = rk4_prop_(t_grid, X0, C, a_extra)
% Fixed-step RK4 with proximity substepping (10 substeps inside
% 15,000 km of Earth or Moon, 4 inside 40,000 km): accuracy through the
% perigee/perilune passages on a fixed filter grid.  RHS inline for speed.
muE = C.mu_E;  muM = C.mu_M;  muS = C.mu_S;
Ja2 = C.J*C.a_E^2;
DEM = C.D_EM;  DSE = C.D_SE;
omM = C.om_M;  omS = C.om_S;  thM0 = C.thM0;  thS0 = C.thS0;

N = numel(t_grid);  X = zeros(6, N);  X(:,1) = X0;  x = X0;
for k = 1:N-1
    tk = t_grid(k);  dt = t_grid(k+1) - tk;
    thM = thM0 + omM*tk;
    dE2 = x(1)^2 + x(2)^2 + x(3)^2;
    dM2 = (x(1)-DEM*cos(thM))^2 + (x(2)-DEM*sin(thM))^2 + x(3)^2;
    d2 = min(dE2, dM2);
    if d2 < 2.25e8, nsub = 10; elseif d2 < 1.6e9, nsub = 4; else, nsub = 1; end
    h = dt/nsub;
    for s = 1:nsub
        ts = tk + (s-1)*h;
        k1 = eom_(ts,     x,          muE,muM,muS,Ja2,DEM,DSE,omM,omS,thM0,thS0,a_extra);
        k2 = eom_(ts+h/2, x + h/2*k1, muE,muM,muS,Ja2,DEM,DSE,omM,omS,thM0,thS0,a_extra);
        k3 = eom_(ts+h/2, x + h/2*k2, muE,muM,muS,Ja2,DEM,DSE,omM,omS,thM0,thS0,a_extra);
        k4 = eom_(ts+h,   x + h  *k3, muE,muM,muS,Ja2,DEM,DSE,omM,omS,thM0,thS0,a_extra);
        x = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
    end
    X(:,k+1) = x;
end
end

function dx = eom_(t, x, muE,muM,muS,Ja2,DEM,DSE,omM,omS,thM0,thS0,ae)
% Equations of motion of [R135] App. A, (A1)-(A3): oblate Earth (second
% harmonic only) + spherical Moon and Sun in differential ("tidal") form.
rx = x(1); ry = x(2); rz = x(3);
r2 = rx*rx + ry*ry + rz*rz;  rn = sqrt(r2);  r3 = r2*rn;
Jf = Ja2/r2;  zf = 5*rz*rz/r2;
cE = -muE/r3;
aX = cE*rx*(1 + Jf*(1 - zf));
aY = cE*ry*(1 + Jf*(1 - zf));
aZ = cE*rz*(1 + Jf*(3 - zf));
thM = thM0 + omM*t;  mx = DEM*cos(thM);  my = DEM*sin(thM);
d1x = rx-mx; d1y = ry-my;
dm3 = (d1x*d1x + d1y*d1y + rz*rz)^1.5;  rm3 = DEM^3;
aX = aX - muM*(d1x/dm3 + mx/rm3);
aY = aY - muM*(d1y/dm3 + my/rm3);
aZ = aZ - muM*(rz /dm3);
thS = thS0 + omS*t;  sx = DSE*cos(thS);  sy = DSE*sin(thS);
d2x = rx-sx; d2y = ry-sy;
ds3 = (d2x*d2x + d2y*d2y + rz*rz)^1.5;  rs3 = DSE^3;
aX = aX - muS*(d2x/ds3 + sx/rs3);
aY = aY - muS*(d2y/ds3 + sy/rs3);
aZ = aZ - muS*(rz /ds3);
dx = [x(4); x(5); x(6); aX+ae(1); aY+ae(2); aZ+ae(3)];
end

function G = grav_grad_(tk, r, C)
% Gravity gradient d(accel)/d(pos), Earth + Moon point masses: the block
% that populates F of [R135] App. B, eqn (B7).  (J2/Sun gradient terms are
% second order here and omitted from F; they remain in the propagation.)
rM = ephem_moon_(tk, C);
G  = pmg_(r, C.mu_E) + pmg_(r - rM, C.mu_M);
end
function Gb = pmg_(d, mu)
dn = norm(d);
Gb = mu*(3*(d*d.')/dn^5 - eye(3)/dn^3);
end

% =========================================================================
%  FREE-RETURN TARGETING
% =========================================================================
function [p, info] = shoot_free_return_(C, tgt, cfg)
% Ballistic circumlunar free-return targeting: a coarse scan locates the
% corridor (or a recorded warm start skips it); nested secant iterations
% then converge it.  Control separation: injection speed v0 sets the
% perilune height (thousands of km per m/s of injection speed), while
% the injection angle th_v shapes the return perigee.  Nested 1-D
% secants are far more robust here than a raw 2-D Newton in this stiff,
% narrow corridor.
if cfg.do_scan
    best = struct('cost', inf, 'p', [nan;nan]);
    for v0 = 10.86:0.010:11.00
        for thv = deg2rad(-165:2.5:-95)
            d = fly_([v0; thv], C, 600);
            % require a RETROGRADE pass about the Moon (hz < 0): the
            % branch that bends the trajectory back toward Earth -- the
            % figure-8.  Prograde passes sling the vehicle outward.
            if d.ok_pl && d.hz < 0
                cost = abs(d.r_pl - tgt.r_perilune)/1000;
                if d.ok_pe, cost = cost + abs(d.r_pe - tgt.r_perigee)/1000;
                else,       cost = cost + 50; end
                if cost < best.cost, best.cost = cost; best.p = [v0; thv]; end
            end
        end
    end
    if ~isfinite(best.cost)
        error('apollo_lkf:scan', 'no free-return candidate in scan range');
    end
    info.scan_best = best;
    p = best.p;
else
    info.scan_best = struct('cost', nan, 'p', cfg.tli_warm);
    p = cfg.tli_warm;
end
[p, h1] = nested_(p, C, tgt, 120);
[p, h2] = nested_(p, C, tgt, 60);
info.history  = [h1; h2];
info.detail   = fly_(p, C, 60);
info.residual = [info.detail.r_pl - tgt.r_perilune;
                 info.detail.r_pe - tgt.r_perigee];
end

function [p, hist] = nested_(p, C, tgt, dt)
% Outer secant on th_v (return perigee); inner secant on v0 (perilune).
% On the retrograde branch the return perigee GROWS with th_v, so an
% iterate that loses its return leg (or its perilune) recovers by
% stepping th_v downward, back into the corridor.
v0 = p(1);  thv = p(2);  hist = [];
g_prev = [];  thv_prev = [];
for it = 1:25
    v0 = solve_v0_(thv, v0, C, tgt, dt);
    d  = fly_([v0; thv], C, dt);
    if ~d.ok_pe || abs(d.r_pl - tgt.r_perilune) > 1000
        thv = thv - deg2rad(1.0);               % recovery step
        g_prev = [];  thv_prev = [];
        continue;
    end
    g = d.r_pe - tgt.r_perigee;
    hist(end+1,:) = [dt, v0, thv, d.r_pl-tgt.r_perilune, g]; %#ok<AGROW>
    fprintf('      [targeting dt=%3ds] v0=%.5f  th=%.3f deg  perilune miss %+8.2f km  perigee miss %+9.1f km\n', ...
            dt, v0, rad2deg(thv), d.r_pl-tgt.r_perilune, g);
    if abs(g) < tgt.tol, break; end
    if isempty(g_prev)
        step = deg2rad(0.2)*sign(-g);
    else
        step = -g*(thv - thv_prev)/(g - g_prev);
        step = max(min(step, deg2rad(2)), -deg2rad(2));
    end
    g_prev = g;  thv_prev = thv;
    thv = thv + step;
end
p = [v0; thv];
end

function v0 = solve_v0_(thv, v0, C, tgt, dt)
f  = @(v) pl_of_(v, thv, C, dt) - tgt.r_perilune;
f0 = f(v0);  v1 = v0 + 2e-4;  f1 = f(v1);
for it = 1:10
    if abs(f1) < 0.2, break; end
    vn = v1 - f1*(v1 - v0)/(f1 - f0);
    vn = max(min(vn, v1 + 5e-3), v1 - 5e-3);
    v0 = v1;  f0 = f1;  v1 = vn;  f1 = f(v1);
end
v0 = v1;
end

function r = pl_of_(v0, thv, C, dt)
d = fly_([v0; thv], C, dt);
if d.ok_pl, r = d.r_pl; else, r = 1e6; end
end

function d = fly_(p, C, dt)
% Propagate injection conditions; locate perilune and return perigee,
% and record the pass direction about the Moon (z angular momentum hz
% in Moon-relative coordinates: hz < 0 = retrograde = figure-8 branch).
tg = 0:dt:175*C.hr;
X  = rk4_prop_(tg, tli_state_(p, C), C, [0;0;0]);
[rM, vM] = ephem_moon_(tg, C);
dMn = sqrt(sum((X(1:3,:) - rM).^2, 1));
dEn = sqrt(sum(X(1:3,:).^2, 1));
d = struct('ok_pl',false,'ok_pe',false,'t_pl',nan,'r_pl',nan, ...
           't_pe',nan,'r_pe',nan,'hz',nan);
[~, kM] = min(dMn);
if kM <= 2 || kM >= numel(tg)-2 || dMn(kM) > 40000, return; end
[d.t_pl, d.r_pl] = qmin_(tg(kM-1:kM+1), dMn(kM-1:kM+1));
d.ok_pl = true;
dr = X(1:3,kM) - rM(:,kM);  dv = X(4:6,kM) - vM(:,kM);
d.hz = dr(1)*dv(2) - dr(2)*dv(1);
seg = kM+3:numel(tg);
[~, kE] = min(dEn(seg));  kE = seg(1) - 1 + kE;
if kE <= kM+3 || kE >= numel(tg)-1 || dEn(kE) > 3e5, return; end
[d.t_pe, d.r_pe] = qmin_(tg(kE-1:kE+1), dEn(kE-1:kE+1));
d.ok_pe = true;
end

function [tm, fm] = qmin_(t3, f3)
c  = polyfit(t3 - t3(2), f3, 2);
tm = t3(2) - c(2)/(2*c(1));
fm = polyval(c, tm - t3(2));
end

% =========================================================================
%  MEASUREMENTS -- R-135 Appendix C space angles
% =========================================================================
function MEAS = make_meas_(t, Xtrue, C, cfg)
% Observables ([R135] App. C, eqns (C1)): the direction of the
% vehicle-Earth line (angles alpha_e, beta_e) and the Earth subtended
% angle (gamma_e = asin(R0/R)); the same triple is formed for the Moon
% when the vehicle is inside the lunar sphere of influence.  The true
% Earth gamma is generated with the horizon-altitude bias b_horizon
% ([A13] sec. 3.2), and every angle carries independent Gaussian noise
% (sigma = 10 arc-sec [BELL]) as assumed in [R135] p. 17.
%
% Two authentic outage effects:
%   - first sighting 1/2 h after injection ([R135] p. 16 schedule);
%   - Earth angles are unavailable while the Moon occults the Earth
%     line-of-sight (near perilune) -- checked geometrically.
N = numel(t);
k_meas = find(mod(t, cfg.dt_meas) == 0 & t >= cfg.t_first);
MEAS = struct('k', {{}}, 'z', {{}}, 'body', {{}}, 'idx', k_meas, ...
              'n_epoch', 0, 'n_scalar', 0, 'n_occulted', 0);
MEAS.k = cell(1, N);  MEAS.z = cell(1, N);  MEAS.body = cell(1, N);
rM = ephem_moon_(t, C);
for k = k_meas
    r  = Xtrue(1:3,k);
    zb = [];  bb = [];
    % --- Earth angles (drop when Moon occults the Earth LOS) -------------
    if ~occulted_(r, rM(:,k), C)
        hE = angles_(r, C.R0 + cfg.b_horizon);
        zb = [zb; hE];  bb = [bb; [1;1;1]];      
    else
        MEAS.n_occulted = MEAS.n_occulted + 3;
    end
    % --- Moon angles inside the sphere of influence -----------------------
    if norm(r - rM(:,k)) < C.SOI_M
        hM = angles_(r - rM(:,k), C.R_M);
        zb = [zb; hM];  bb = [bb; [2;2;2]];     
    end
    if ~isempty(zb)
        MEAS.k{k}    = k;
        MEAS.z{k}    = zb + cfg.sig_ang*randn(size(zb));
        MEAS.body{k} = bb;
        MEAS.n_epoch  = MEAS.n_epoch + 1;
        MEAS.n_scalar = MEAS.n_scalar + numel(zb);
    end
end
end

function h = angles_(d, Rref)
% [R135] App. C, eqns (C1): alpha = -asin(Z/R); beta from -Y,-X (the
% atan2 form reproduces the tabulated partials); gamma = asin(Rref/R).
Rn = norm(d);
h  = [ -asin(d(3)/Rn);
        atan2(-d(2), -d(1));
        asin(min(Rref/Rn, 1)) ];
end

function [H, hpred] = angles_H_(d, Rref, bias_col)
% Measurement partials, [R135] App. C table (p. 34), plus the horizon-
% bias partial d(gamma)/d(b) = 1/sqrt(R^2 - Rref^2) when bias_col is true.
% Rows: [alpha; beta; gamma]; columns of H: 10-state
% [dr(3), dv(3), a_u(3), b_h].
X = d(1); Y = d(2); Z = d(3);
R2 = X*X + Y*Y + Z*Z;  Rn = sqrt(R2);
rho = sqrt(X*X + Y*Y);
Ha = [ X*Z/(R2*rho),  Y*Z/(R2*rho),  (Z*Z - R2)/(R2*rho) ];
Hb = [ -Y/rho^2,      X/rho^2,       0                    ];
sg = sqrt(max(R2 - Rref^2, 1e-6));
Hg = [ -Rref*X/(R2*sg), -Rref*Y/(R2*sg), -Rref*Z/(R2*sg) ];
H  = zeros(3, 10);
H(1,1:3) = Ha;  H(2,1:3) = Hb;  H(3,1:3) = Hg;
if bias_col
    H(3,10) = 1/sg;                 % d(gamma_e)/d(b_horizon)
end
hpred = [ -asin(Z/Rn); atan2(-Y, -X); asin(min(Rref/Rn, 1)) ];
end

function tf = occulted_(r, rMk, C)
% Earth LOS from the vehicle: blocked when the Moon lies between the
% vehicle and the Earth within one lunar radius of the line of sight.
u  = -r/norm(r);              % vehicle -> Earth-center direction
m  = rMk - r;                 % vehicle -> Moon
s  = dot(m, u);
tf = (s > 0) && (s < norm(r)) && (norm(m - s*u) < C.R_M);
end

% =========================================================================
%  THE FILTER -- [R135] eqns (12)-(16), Joseph-stabilized
% =========================================================================
function F = run_filter_(mode, t, Xnom, MEAS, P0, C, cfg)
% mode 'LKF': linearize F and H about the precomputed NOMINAL trajectory
%             and estimate the deviation state -- the filter exactly as
%             formulated in [R135].
% mode 'EKF': linearize about the CURRENT ESTIMATE and propagate the
%             estimate through the full nonlinear dynamics -- the
%             modification recommended on p. 15 of [R135] ("linearize
%             around the estimated rather than the reference trajectory"),
%             identified in [TM847] as the extended Kalman filter.
%
% State (10): [dr (km); dv (km/s); a_u (km/s^2); b_h (km)]
% Updates are processed one scalar angle at a time: with independent
% angle errors the innovation "matrix to be inverted is 1x1 ... the
% ultimate in calculation simplicity" ([R135] p. 23).  Covariance uses
% the Joseph form with explicit re-symmetrization.
n  = 10;
N  = numel(t);
dt = t(2) - t(1);
x  = zeros(n,1);                    % deviation estimate (LKF frame)
P  = P0;
Xe = zeros(6,N);                    % total position/velocity estimate
sig = zeros(n,N);
au_hist = zeros(4,N);               % [a_u; b_h] estimates
Prr6 = zeros(6,N);                  % packed position covariance (NEES)
innov = nan(2,N);                   % [gammaE innovation; sqrt(S)] record

Xtot = Xnom(:,1);                   % EKF total state (starts on nominal)
R1 = cfg.sig_ang^2;

for k = 1:N
    % ---------------- time update ---------------------------------------
    if k > 1
        if strcmp(mode, 'LKF')
            rlin = Xnom(1:3,k-1);
        else
            rlin = Xtot(1:3);
        end
        G = grav_grad_(t(k-1), rlin, C);
        A = [zeros(3), eye(3),   zeros(3),   zeros(3,1);
             G,        zeros(3), eye(3),     zeros(3,1);
             zeros(3), zeros(3), zeros(3),   zeros(3,1);
             zeros(1,9),                     0          ];
        Phi = expm(A*dt);                       % [R135] App. E role
        Qd  = qd_(dt, cfg);
        P   = Phi*P*Phi.' + Qd;
        P   = (P + P.')/2;
        if strcmp(mode, 'LKF')
            x = Phi*x;
        else
            % propagate the total estimate through the nonlinear dynamics,
            % including the currently estimated disturbance acceleration
            Xseg = rk4_prop_([t(k-1) t(k)], Xtot, C, x(7:9));
            Xtot = Xseg(:,2);
        end
    end

    % ---------------- measurement update (scalar, sequential) -----------
    if ~isempty(MEAS.k{k})
        zb = MEAS.z{k};  bb = MEAS.body{k};
        rM = ephem_moon_(t(k), C);
        for j = 1:3:numel(zb)
            body = bb(j);
            if strcmp(mode, 'LKF')
                % pure linearized filter: h and H evaluated on the
                % precomputed nominal, bias handled linearly through H
                rr = Xnom(1:3,k);
                if body == 1, [H, hp] = angles_H_(rr, C.R0, true);
                else,         [H, hp] = angles_H_(rr - rM, C.R_M, false);
                end
            else
                % estimate-linearized: h and H at the current estimate,
                % including the current horizon-bias estimate
                rr = Xtot(1:3);
                if body == 1, [H, hp] = angles_H_(rr, C.R0 + x(10), true);
                else,         [H, hp] = angles_H_(rr - rM, C.R_M, false);
                end
            end
            for i = 0:2
                zi = zb(j+i);  Hi = H(i+1,:);
                if strcmp(mode, 'LKF')
                    ino = wrap_(zi - hp(i+1) - Hi*x);
                else
                    ino = wrap_(zi - hp(i+1));
                end
                S = Hi*P*Hi.' + R1;
                K = (P*Hi.')/S;
                upd = K*ino;
                if strcmp(mode, 'LKF')
                    x = x + upd;
                else
                    Xtot = Xtot + upd(1:6);
                    x(7:10) = x(7:10) + upd(7:10);
                end
                IKH = eye(n) - K*Hi;
                P = IKH*P*IKH.' + K*R1*K.';     % Joseph form
                P = (P + P.')/2;                % symmetry maintenance
                if body == 1 && i == 2 && isnan(innov(1,k))
                    innov(:,k) = [ino; sqrt(S)];
                end
            end
        end
    end

    % ---------------- bookkeeping ----------------------------------------
    if strcmp(mode, 'LKF')
        Xe(:,k) = Xnom(:,k) + x(1:6);
        au_hist(:,k) = x(7:10);
    else
        Xe(:,k) = Xtot;
        au_hist(:,k) = x(7:10);
    end
    sig(:,k) = sqrt(max(diag(P), 0));
    Prr6(:,k) = [P(1,1); P(2,1); P(2,2); P(3,1); P(3,2); P(3,3)];
end

F = struct('mode',mode, 'Xest',Xe, 'sig',sig, 'au',au_hist, ...
           'innov',innov, 'Prr6',Prr6);
end

function Qd = qd_(dt, cfg)
% Discrete process noise: white acceleration PSD q_acc on the velocity
% states (standard PWN discretization) + slow random walk q_adot on the
% disturbance-acceleration states; the horizon bias is a constant.
qa = cfg.q_acc;   qw = cfg.q_adot;
Qrr = qa*dt^3/3*eye(3);  Qrv = qa*dt^2/2*eye(3);  Qvv = qa*dt*eye(3);
Qd = [Qrr,      Qrv,      zeros(3), zeros(3,1);
      Qrv,      Qvv,      zeros(3), zeros(3,1);
      zeros(3), zeros(3), qw*dt*eye(3), zeros(3,1);
      zeros(1,9),                       0         ];
end

function a = wrap_(a)
a = atan2(sin(a), cos(a));
end

function seed_(s)
% reproducible noise in both MATLAB and Octave
try
    rng(s);
catch
    randn('state', s);  rand('state', s);      
end
end

% =========================================================================
%  FIGURES
% =========================================================================
function make_figures_(R)
t = R.t; C = R.C; hr = C.hr;
tpl = R.tinfo.detail.t_pl;
th  = t/hr;

% ---- fig 1: the free-return figure-8 --------------------------------------
f = dark_fig_([1100 520]);
subplot(1,2,1); hold on; axis equal;
h1 = plot(R.rM(1,:), R.rM(2,:), '--', 'Color', [0.45 0.45 0.45]);
h2 = plot(R.Xnom(1,:), R.Xnom(2,:), '-', 'Color', [0.3 1 0.4], 'LineWidth', 1.4);
draw_body_(0, 0, C.R0*4, [0.25 0.55 1]);
kpl = round(tpl/(t(2)-t(1))) + 1;
draw_body_(R.rM(1,kpl), R.rM(2,kpl), C.R_M*4, [0.75 0.75 0.75]);
h3 = plot(R.Xnom(1,1), R.Xnom(2,1), 'w^', 'MarkerFaceColor', 'w', 'MarkerSize', 5);
title('Earth-centered inertial frame', 'Color', 'w');
xlabel('X [km]'); ylabel('Y [km]'); dark_ax_();
legend([h1 h2 h3], {'Moon orbit','free-return trajectory','TLI'}, ...
       'TextColor','w','Color','none','Location','northwest');
subplot(1,2,2); hold on; axis equal;
% Moon-corotating frame: rotate by -theta_M(t); the figure-8 appears
Xr = rot_frame_(R.Xnom(1:3,:), t, C);
Xt = rot_frame_(R.Xtrue(1:3,:), t, C);
h1 = plot(Xr(1,:), Xr(2,:), '-', 'Color', [0.3 1 0.4], 'LineWidth', 1.4);
h2 = plot(Xt(1,:), Xt(2,:), '-', 'Color', [1 0.45 0.35], 'LineWidth', 0.8);
draw_body_(0, 0, C.R0*4, [0.25 0.55 1]);
draw_body_(C.D_EM, 0, C.R_M*4, [0.75 0.75 0.75]);
xlim([-0.6e5 4.4e5]); ylim([-2.0e5 2.0e5]);
title('Moon-corotating frame: the figure-8', 'Color', 'w');
xlabel('X [km]'); ylabel('Y [km]'); dark_ax_();
legend([h1 h2], {'nominal','true (uncorrected)'}, 'TextColor','w','Color','none','Location','northwest');
print(f, '-dpng', '-r150', 'fig1_trajectory.png'); close(f);

% ---- fig 2: error amplification vs the filter ----------------------------
f = light_fig_([900 420]);
dev  = sqrt(sum((R.Xtrue(1:3,:) - R.Xnom(1:3,:)).^2,1));
eL   = sqrt(sum((R.FL.Xest(1:3,:) - R.Xtrue(1:3,:)).^2,1));
eE   = sqrt(sum((R.FE.Xest(1:3,:) - R.Xtrue(1:3,:)).^2,1));
semilogy(th, dev, 'k-', 'LineWidth', 1.4); hold on; grid on;
semilogy(th, eL, 'b-', 'LineWidth', 1.1);
semilogy(th, eE, 'r-', 'LineWidth', 1.1);
xline_(tpl/hr);
xlabel('mission time [h]'); ylabel('position magnitude [km]');
title('Deviation from nominal (open loop) vs filter estimation error');
legend({'|r_{true} - r_{nom}| (no navigation)', ...
        '|r_{est} - r_{true}|  LKF (nominal-linearized)', ...
        '|r_{est} - r_{true}|  EKF variant (estimate-linearized)', ...
        'perilune'}, 'Location', 'northwest');
print(f, '-dpng', '-r150', 'fig2_amplification.png'); close(f);

% ---- fig 3 / fig 4: errors vs 3-sigma, both filters ------------------------
% (the LKF panels are clipped to the pre-breakdown scale; its final
%  divergence is documented by fig 2)
plot_err3sig_(R, 'FL', 'fig3_lkf_errors.png', ...
    'LKF (linearized about the nominal): error vs \pm3\sigma  --  diverges off-scale near entry', 130);
plot_err3sig_(R, 'FE', 'fig4_ekf_errors.png', ...
    'Estimate-linearized variant (EKF): error vs \pm3\sigma', inf);

% ---- fig 5: augmented parameter states ------------------------------------
f = light_fig_([950 620]);
lbl = {'a_{u,x}','a_{u,y}','a_{u,z}'};
for i = 1:3
    subplot(2,2,i); hold on; grid on;
    ylim([-40 40]);  xline_(tpl/hr);
    hf = fill_sig_(th, R.FE.au(i,:), 3*R.FE.sig(6+i,:), [1 0.85 0.8], 1e9);
    h1 = plot(th, R.FE.au(i,:)*1e9, 'r-', 'LineWidth', 1.1);
    h2 = plot(th, R.FL.au(i,:)*1e9, 'b-', 'LineWidth', 0.9);
    h3 = plot(th([1 end]), R.cfg.a_vent(i)*[1 1]*1e9, 'k--');
    xlabel('t [h]'); ylabel([lbl{i} ' [10^{-9} km/s^2]']);
    ylim([-40 40]);
    if i == 1
        legend([hf h1 h2 h3], {'\pm3\sigma (EKF)','EKF est.','LKF est.','truth'}, ...
               'Location','northeast');
    end
end
subplot(2,2,4); hold on; grid on;
ylim([-10 30]);  xline_(tpl/hr);
fill_sig_(th, R.FE.au(4,:), 3*R.FE.sig(10,:), [1 0.85 0.8]);
plot(th, R.FE.au(4,:), 'r-', 'LineWidth', 1.1);
plot(th, R.FL.au(4,:), 'b-', 'LineWidth', 0.9);
plot(th([1 end]), R.cfg.b_horizon*[1 1], 'k--');
xlabel('t [h]'); ylabel('horizon bias b_h [km]');
ylim([-10 30]);
sgtitle_('Augmented states: venting acceleration and Earth-horizon bias (LKF runs off-scale after the flyby)');
print(f, '-dpng', '-r150', 'fig5_parameters.png'); close(f);

% ---- fig 6: innovations ----------------------------------------------------
f = light_fig_([900 420]);
kk = find(~isnan(R.FE.innov(1,:)));
plot(th(kk), R.FE.innov(1,kk)./R.FE.innov(2,kk), '.', 'MarkerSize', 5); hold on; grid on;
plot(th([1 end]), [ 3  3; -3 -3].', 'r--');
xline_(tpl/hr);
ylim([-6 6]);
xlabel('mission time [h]'); ylabel('\gamma_E innovation / \surdS');
title('Normalized subtended-Earth-angle innovations (EKF variant)');
print(f, '-dpng', '-r150', 'fig6_innovations.png'); close(f);
end

function plot_err3sig_(R, fld, fname, ttl, tcap)
t = R.t; th = t/R.C.hr;  F = R.(fld);
err = F.Xest - R.Xtrue;
lbl = {'x [km]','y [km]','z [km]','v_x [m/s]','v_y [m/s]','v_z [m/s]'};
scl = [1 1 1 1e3 1e3 1e3];
kc  = th <= tcap;
f = light_fig_([1000 620]);
for i = 1:6
    subplot(2,3,i); hold on; grid on;
    cap = 1.2*max( max(abs(err(i,kc)))*scl(i), max(3*F.sig(i,kc))*scl(i) );
    ylim([-cap cap]);  xline_(R.tinfo.detail.t_pl/R.C.hr);
    fill_sig_(th, 0*th, 3*F.sig(i,:)*scl(i), [0.82 0.88 1], 1);
    plot(th, err(i,:)*scl(i), 'b-', 'LineWidth', 0.7);
    ylim([-cap cap]);
    xlabel('t [h]'); ylabel(lbl{i});
end
sgtitle_(ttl);
print(f, '-dpng', '-r150', fname); close(f);
end

% ---- small plotting utilities ---------------------------------------------
function f = dark_fig_(sz)
f = figure('Visible','off','Color','k','Position',[50 50 sz], ...
           'InvertHardcopy','off');
end
function f = light_fig_(sz)
f = figure('Visible','off','Color','w','Position',[50 50 sz]);
end
function dark_ax_()
set(gca, 'Color','k','XColor','w','YColor','w','GridColor',[0.35 0.35 0.35]);
grid on;
end
function draw_body_(x, y, rad, col)
th = linspace(0, 2*pi, 60);
fill(x + rad*cos(th), y + rad*sin(th), col, 'EdgeColor', 'none');
end
function Xr = rot_frame_(Xp, t, C)
th = C.thM0 + C.om_M*t;
Xr = zeros(2, numel(t));
Xr(1,:) =  cos(th).*Xp(1,:) + sin(th).*Xp(2,:);
Xr(2,:) = -sin(th).*Xp(1,:) + cos(th).*Xp(2,:);
end
function h = fill_sig_(th, mid, s3, col, scl)
if nargin < 5, scl = 1; end
xx = [th, fliplr(th)];
yy = [(mid + s3)*scl, fliplr((mid - s3)*scl)];
h = fill(xx, yy, col, 'EdgeColor', 'none');
end
function xline_(xv)
yl = ylim;
plot([xv xv], yl, ':', 'Color', [0.4 0.4 0.4]);
end
function sgtitle_(s)
try
    sgtitle(s);
catch
    axes('Position',[0 0.96 1 0.04],'Visible','off');
    text(0.5, 0.5, s, 'HorizontalAlignment','center','FontWeight','bold');
end
end
