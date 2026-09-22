function apollo_lkf_animation(resfile, outname)
%APOLLO_LKF_ANIMATION  Animated reconstruction of the free-return mission.
%
%   apollo_lkf_animation()                       uses apollo_lkf_results.mat
%   apollo_lkf_animation(resfile, outname)       custom input / output
%
%   Renders the mission produced by APOLLO_LKF_MISSION as an MP4 video
%   (MATLAB: VideoWriter; if no video backend is available, individual
%   frames are written to ./anim_frames for external assembly).
%
%   Three synchronized views:
%     LEFT    Earth-centered inertial frame, top-down: the Moon moves on
%             its orbit while the spacecraft flies the ballistic
%             free-return "figure-8" (trajectory shape per the equations
%             of motion of NASA TR R-135, App. A).
%     TOP R.  Moon-corotating frame: the figure-8 draws itself.
%     BOT R.  Navigation view, centered on the true spacecraft: nominal,
%             true, LKF and EKF estimates with the EKF 3-sigma position
%             ellipse (enlarged 30x) -- the filter's own statement of
%             its uncertainty, cf. R-135's error-ellipsoid discussion.
%
%   Mission clock and phase captions follow the style of the author's
%   earlier animation; bodies are enlarged 3x for visibility.

if nargin < 1 || isempty(resfile), resfile = 'apollo_lkf_results.mat'; end
if nargin < 2 || isempty(outname), outname = 'Apollo_FreeReturn_LKF.mp4'; end
R = load(resfile);
C = R.C;  t = R.t;  hr = C.hr;  N = numel(t);

% ---- frame schedule: one frame per 15 min of mission time ----------------
stride = max(1, round(900/(t(2)-t(1))));
kk = 1:stride:N;  if kk(end) ~= N, kk(end+1) = N; end

% ---- derived signals -------------------------------------------------------
dMoon  = sqrt(sum((R.Xtrue(1:3,:) - R.rM).^2, 1));
dev    = sqrt(sum((R.Xtrue(1:3,:) - R.Xnom(1:3,:)).^2, 1));
errE   = sqrt(sum((R.FE.Xest(1:3,:) - R.Xtrue(1:3,:)).^2, 1));
s3pos  = 3*sqrt(sum(R.FE.sig(1:3,:).^2, 1));      % 3-sigma position norm
inSOI  = dMoon < C.SOI_M;
k_in   = find(inSOI, 1);  k_out = find(inSOI, 1, 'last');
Xr_nom = corot_(R.Xnom(1:3,:), t, C);
Xr_tru = corot_(R.Xtrue(1:3,:), t, C);

% ---- occultation flags (Earth LOS blocked by the Moon) --------------------
occ = false(1, N);
for k = 1:N
    r = R.Xtrue(1:3,k);  m = R.rM(:,k) - r;  u = -r/norm(r);
    s = dot(m, u);
    occ(k) = (s > 0) && (s < norm(r)) && (norm(m - s*u) < C.R_M);
end

% ---- figure and static scenery --------------------------------------------
f = figure('Visible','off','Color','k','Position',[40 40 1280 720], ...
           'InvertHardcopy','off');

axL = axes('Position',[0.045 0.07 0.56 0.86]);  hold on; axis equal;
set(axL,'Color','k','XColor',[0.6 0.6 0.6],'YColor',[0.6 0.6 0.6], ...
        'GridColor',[0.25 0.25 0.25]); grid on;
plot(axL, C.D_EM*cos(0:0.02:2*pi), C.D_EM*sin(0:0.02:2*pi), ':', ...
     'Color',[0.35 0.35 0.35]);
plot(axL, R.Xnom(1,:), R.Xnom(2,:), '--', 'Color',[0.30 0.45 0.32], 'LineWidth',0.6);
body_(axL, 0, 0, C.R0*3, [0.25 0.55 1]);
hMoonL  = body_(axL, R.rM(1,1), R.rM(2,1), C.R_M*3, [0.72 0.72 0.72]);
hTrail  = plot(axL, nan, nan, '-',  'Color',[0.35 1 0.45], 'LineWidth',1.6);
hEstL   = plot(axL, nan, nan, '--', 'Color',[0.4 0.85 1], 'LineWidth',0.9);
hShip   = plot(axL, nan, nan, 'o', 'MarkerSize',5, ...
               'MarkerFaceColor',[1 0.35 0.3], 'MarkerEdgeColor','w');
lim = 4.6e5;  xlim(axL, [-lim lim]);  ylim(axL, [-lim lim]);
xlabel(axL,'X [km]'); ylabel(axL,'Y [km]');
hClock = text(axL, -lim*0.95, lim*0.92, '', 'Color','w', ...
              'FontSize',13, 'FontWeight','bold');
hPhase = text(axL, -lim*0.95, lim*0.83, '', 'Color','c', 'FontSize',11);
hInfo  = text(axL, -lim*0.95, -lim*0.97, '', 'Color',[0.8 0.8 0.8], ...
              'FontSize',9, 'VerticalAlignment','bottom');
title(axL, 'Apollo free-return mission -- LKF navigation reconstruction', ...
      'Color','w', 'FontSize', 11);

axR1 = axes('Position',[0.66 0.55 0.32 0.38]);  hold on; axis equal;
set(axR1,'Color','k','XColor',[0.6 0.6 0.6],'YColor',[0.6 0.6 0.6], ...
         'GridColor',[0.25 0.25 0.25]); grid on;
body_(axR1, 0, 0, C.R0*3, [0.25 0.55 1]);
body_(axR1, C.D_EM, 0, C.R_M*3, [0.72 0.72 0.72]);
plot(axR1, Xr_nom(1,:), Xr_nom(2,:), ':', 'Color',[0.32 0.42 0.34]);
hFig8 = plot(axR1, nan, nan, '-', 'Color',[0.35 1 0.45], 'LineWidth',1.2);
hDot8 = plot(axR1, nan, nan, 'o', 'MarkerSize',4, ...
             'MarkerFaceColor',[1 0.35 0.3], 'MarkerEdgeColor','none');
xlim(axR1, [-0.7e5 4.5e5]); ylim(axR1, [-1.9e5 1.9e5]);
title(axR1, 'Moon-corotating frame', 'Color','w', 'FontSize',9);

axR2 = axes('Position',[0.66 0.08 0.32 0.38]);  hold on; axis equal;
set(axR2,'Color','k','XColor',[0.6 0.6 0.6],'YColor',[0.6 0.6 0.6], ...
         'GridColor',[0.25 0.25 0.25]); grid on;
hZnom = plot(axR2, R.Xnom(1,:),   R.Xnom(2,:),   '--', 'Color',[0.45 0.45 0.45]);
hZtru = plot(axR2, R.Xtrue(1,:),  R.Xtrue(2,:),  '-',  'Color',[0.35 1 0.45]);
hZlkf = plot(axR2, R.FL.Xest(1,:),R.FL.Xest(2,:),'-',  'Color',[1 0.7 0.2]);
hZekf = plot(axR2, R.FE.Xest(1,:),R.FE.Xest(2,:),'-',  'Color',[0.4 0.85 1]);
hZell = plot(axR2, nan, nan, '-', 'Color',[0.4 0.85 1], 'LineWidth',1.1);
hZs   = plot(axR2, nan, nan, 'o', 'MarkerSize',5, ...
             'MarkerFaceColor',[1 0.35 0.3], 'MarkerEdgeColor','w');
hZbar = plot(axR2, nan, nan, '-', 'Color','w', 'LineWidth',2);
hZbtx = text(axR2, 0, 0, '', 'Color','w', 'FontSize',8);
set(axR2, 'XTickLabel', [], 'YTickLabel', []);
hZttl = title(axR2, 'navigation view', 'Color','w', 'FontSize',9);
lg = legend(axR2, [hZnom hZtru hZlkf hZekf], ...
       {'nominal','true','LKF','EKF'}, 'TextColor','w', 'Color','none', ...
       'EdgeColor',[0.4 0.4 0.4], 'FontSize',7, 'Location','northeast');
try, set(lg, 'AutoUpdate', 'off'); catch, end

% ---- video sink ------------------------------------------------------------
use_vw = true;
try
    vw = VideoWriter(outname, 'MPEG-4');  vw.FrameRate = 30;  vw.Quality = 95;
    open(vw);
catch
    try
        vw = VideoWriter(strrep(outname, '.mp4', '.avi'), 'Motion JPEG AVI');
        vw.FrameRate = 30;  open(vw);
    catch
        use_vw = false;
        if ~exist('anim_frames', 'dir'), mkdir('anim_frames'); end
        fprintf('VideoWriter unavailable: writing PNG frames to ./anim_frames\n');
    end
end

% ---- render loop -----------------------------------------------------------
nf = 0;
for k = kk
    nf = nf + 1;
    ud = get(hMoonL, 'UserData');
    set(hMoonL, 'XData', ud.x + R.rM(1,k), 'YData', ud.y + R.rM(2,k));
    set(hTrail, 'XData', R.Xtrue(1,1:k), 'YData', R.Xtrue(2,1:k));
    set(hEstL,  'XData', R.FE.Xest(1,1:k), 'YData', R.FE.Xest(2,1:k));
    set(hShip,  'XData', R.Xtrue(1,k), 'YData', R.Xtrue(2,k));
    set(hFig8,  'XData', Xr_tru(1,1:k), 'YData', Xr_tru(2,1:k));
    set(hDot8,  'XData', Xr_tru(1,k), 'YData', Xr_tru(2,k));

    set(hClock, 'String', sprintf('Mission time: %6.1f h  (day %.1f)', ...
                                  t(k)/hr, t(k)/86400));
    if k < k_in
        set(hPhase, 'String', 'Phase: OUTBOUND  (Earth \rightarrow Moon)', ...
            'Color', 'c');
    elseif k <= k_out
        if occ(k)
            set(hPhase, 'String', ...
                'Phase: LUNAR FLYBY -- Earth occulted, Moon angles only', ...
                'Color', [1 0.8 0.3]);
        else
            set(hPhase, 'String', ...
                'Phase: LUNAR FLYBY  (Moon sightings active)', 'Color', 'g');
        end
    else
        set(hPhase, 'String', 'Phase: INBOUND  (return to Earth)', ...
            'Color', [1 0.4 0.35]);
    end
    set(hInfo, 'String', sprintf( ...
        ['open-loop deviation from nominal: %8.1f km\n' ...
         'EKF position error: %6.2f km   (3\\sigma: %.1f km)'], ...
        dev(k), errE(k), s3pos(k)));

    % navigation zoom: follow the true spacecraft
    w = max(2500, 1.8*dev(k));
    xlim(axR2, R.Xtrue(1,k) + [-w w]);  ylim(axR2, R.Xtrue(2,k) + [-w w]);
    P2 = [R.FE.Prr6(1,k), R.FE.Prr6(2,k); R.FE.Prr6(2,k), R.FE.Prr6(3,k)];
    [V, D] = eig((P2 + P2.')/2);
    a  = 3*sqrt(max(D(1,1), 0));  b = 3*sqrt(max(D(2,2), 0));
    esc = min(300, 0.55*w/max(max(a, b), 1e-9));  % adaptive ellipse scale
    ph = 0:0.15:2*pi+0.15;
    E  = V*[a*cos(ph); b*sin(ph)]*esc;
    set(hZell, 'XData', R.FE.Xest(1,k) + E(1,:), 'YData', R.FE.Xest(2,k) + E(2,:));
    set(hZs,   'XData', R.Xtrue(1,k), 'YData', R.Xtrue(2,k));
    set(hZttl, 'String', sprintf('navigation view (EKF 3\\sigma ellipse \\times%.0f)', esc));
    % scale bar: a round-number length ~ a third of the window
    bl = 10^floor(log10(0.7*w));
    if 0.7*w/bl >= 5, bl = 5*bl; elseif 0.7*w/bl >= 2, bl = 2*bl; end
    bx = R.Xtrue(1,k) - 0.92*w;  by = R.Xtrue(2,k) - 0.90*w;
    set(hZbar, 'XData', [bx bx+bl], 'YData', [by by]);
    set(hZbtx, 'Position', [bx + bl/2, by + 0.07*w], ...
               'String', fmt_km_(bl), 'HorizontalAlignment', 'center');

    drawnow;
    if use_vw
        writeVideo(vw, getframe(f));
    else
        print(f, '-dpng', '-r110', sprintf('anim_frames/frame_%04d.png', nf));
    end
end
if use_vw
    close(vw);
    fprintf('animation written: %s (%d frames)\n', outname, nf);
else
    fprintf(['frames written (%d).  Assemble with e.g.:\n' ...
             '  ffmpeg -r 30 -i anim_frames/frame_%%04d.png -pix_fmt yuv420p %s\n'], ...
            nf, outname);
end
close(f);
end

% ---------------------------------------------------------------------------
function h = body_(ax, x, y, rad, col)
ph = linspace(0, 2*pi, 48);
h = fill(ax, x + rad*cos(ph), y + rad*sin(ph), col, 'EdgeColor', 'none');
set(h, 'UserData', struct('x', rad*cos(ph), 'y', rad*sin(ph)));
end

function Xr = corot_(Xp, t, C)
th = C.thM0 + C.om_M*t;
Xr = [ cos(th).*Xp(1,:) + sin(th).*Xp(2,:);
      -sin(th).*Xp(1,:) + cos(th).*Xp(2,:)];
end

function s = fmt_km_(v)
if v >= 1000, s = sprintf('%g,000 km', v/1000); else, s = sprintf('%g km', v); end
end
