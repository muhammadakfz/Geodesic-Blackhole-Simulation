function geodesic_blackhole_simulation()

clc; close all;

p = defaultParameters();
[p.kappa, particleLabel] = particleTypeToKappa(p.particleType);

if p.kappa == 1 && p.E <= 1
    error('For proton mode, energy parameter E must satisfy E > 1.');
end

fprintf('=== %s Geodesic in Schwarzschild Spacetime ===\n', particleLabel);
fprintf('Black hole mass M        : %.6f\n', p.M);
fprintf('Energy parameter E       : %.6f\n', p.E);
fprintf('Horizon radius r_h       : %.6f\n', p.r_horizon);

if p.kappa == 0
    bAnalytic = 3 * sqrt(3) * p.M;
    bGuess = bAnalytic;
    fprintf('Analytic b_crit (photon) : %.10f\n', bAnalytic);
else
    bAnalytic = NaN;
    bGuess = p.bGuessProton;
    fprintf('No simple closed-form b_crit for chosen proton E in this script.\n');
end

[bLow, bHigh] = autoBracketCritical(p, bGuess);
[bCritNum, bisectLog] = bisectCriticalImpact(p, bLow, bHigh, p.bisectTol, p.bisectMaxIter);

fprintf('\nBisection bracket start   : [%.10f, %.10f]\n', bLow, bHigh);
fprintf('Bisection iterations      : %d\n', size(bisectLog, 1));
fprintf('Numerical b_crit          : %.10f\n', bCritNum);

if p.kappa == 0
    relErr = abs(bCritNum - bAnalytic) / bAnalytic;
    fprintf('Relative error vs analytic: %.3e\n', relErr);

    arg = bCritNum * sqrt(1 - 2 * p.M / p.rObserver) / p.rObserver;
    if abs(arg) <= 1
        alphaCrit = asin(arg);
        fprintf('Critical local angle at observer (r = %.2f M): %.6f rad (%.4f deg)\n', ...
            p.rObserver / p.M, alphaCrit, rad2deg(alphaCrit));
    else
        fprintf('Observer radius too small for real local critical angle in this setup.\n');
    end
end

if strcmpi(p.visualStyle, 'interactive')
    if p.kappa ~= 0
        warning('Interactive drag mode currently supports photon dynamics only. Falling back to dashboard mode.');
    else
        runInteractivePhotonSandbox(p, particleLabel, bCritNum, bAnalytic);
        return;
    end
end

if strcmpi(p.visualStyle, 'beam')
    if p.kappa ~= 0
        warning('Beam visual style is optimized for photon mode. Falling back to dashboard mode.');
    else
        beam2D = generateBeamTrajectories(p, bCritNum, p.beamNumRays2D, 0.0);
        plotBeamStyle2D(p, beam2D, particleLabel, bCritNum);

        if p.enableAnimation
            animateBeamStyle2D(p, beam2D, particleLabel, bCritNum);
        end

        if p.enable3D
            zOffsets = linspace(-p.beam3DZSpread, p.beam3DZSpread, p.beam3DLayers);
            beam3D = generateBeamTrajectories(p, bCritNum, p.beamNumRays3D, zOffsets);
            plotBeamStyle3D(p, beam3D, particleLabel, bCritNum);

            if p.enableAnimation3D
                animateBeamStyle3D(p, beam3D, particleLabel, bCritNum);
            end
        end
        return;
    end
end

sampleFactors = [0.80, 0.90, 0.97, 0.995, 1.005, 1.03, 1.10];
bSamples = sampleFactors * bCritNum;
trajectories = cell(numel(bSamples), 1);

fprintf('\n--- Sample trajectories around b_crit ---\n');
for i = 1:numel(bSamples)
    trajectories{i} = integrateTrajectory(p, bSamples(i));
    fprintf('b = %.10f -> %s\n', trajectories{i}.b, trajectories{i}.status);
end

styles = buildTrajectoryStyles(trajectories);

if p.useSingleWindow
    dashboard = createDashboardFigure(p);
    populateInfoPanel(dashboard.axInfo, p, particleLabel, bCritNum, bAnalytic, trajectories, styles);
    ax2DStatic = dashboard.ax2DStatic;
    ax2DAnim = dashboard.ax2DAnim;
    ax3DStatic = dashboard.ax3DStatic;
    ax3DAnim = dashboard.ax3DAnim;
    axConv = dashboard.axConv;
    loopFigure = dashboard.fig;
else
    ax2DStatic = [];
    ax2DAnim = [];
    ax3DStatic = [];
    ax3DAnim = [];
    axConv = [];
    loopFigure = [];
end

plotTrajectories2D(p, trajectories, styles, particleLabel, bCritNum, bAnalytic, ax2DStatic);

if p.enable3D
    plotTrajectories3D(p, trajectories, styles, particleLabel, ax3DStatic);
elseif p.useSingleWindow
    showDisabledPanel(ax3DStatic, '3D Trajectories', 'Static 3D plot is disabled');
end

if p.showConvergence
    plotBisectionConvergence(bisectLog, particleLabel, axConv);
elseif p.useSingleWindow
    populateTrajectoryStatsPanel(axConv, trajectories);
end

if p.useSingleWindow && p.animationLoop && isFigureVisible(loopFigure)
    while isFigureVisible(loopFigure)
        runAnimationPass(p, trajectories, styles, particleLabel, bCritNum, bAnalytic, ax2DAnim, ax3DAnim);
        if p.animationLoopDelay > 0 && isFigureVisible(loopFigure)
            pause(p.animationLoopDelay);
        end
    end
else
    runAnimationPass(p, trajectories, styles, particleLabel, bCritNum, bAnalytic, ax2DAnim, ax3DAnim);
end

end


function p = defaultParameters()
    p.visualStyle = 'interactive'; % 'interactive', 'beam' (cinematic), or 'dashboard'

    p.useSingleWindow = true;    % true -> all charts in one figure window
    p.dashboardScale = 1.0;      % Relative size multiplier for single-window dashboard

    p.particleType = 'photon';   % 'photon' or 'proton'

    p.M = 1.0;                   % Black hole mass
    p.E = 1.0;                   % Photon: any positive E, Proton: E > 1 (gamma-like)

    p.r_horizon = 2 * p.M;
    p.r0 = 80 * p.M;             % Initial radial position
    p.r_escape = 70 * p.M;       % Outgoing radius to classify as escaped
    p.r_max = 250 * p.M;         % Safety bound for runaway integration
    p.rObserver = 20 * p.M;      % Static observer radius for angle readout

    p.lambdaMax = 7000 * p.M;
    p.maxStep = 0.5 * p.M;
    p.relTol = 1e-9;
    p.absTol = 1e-10;

    p.bisectTol = 1e-6 * p.M;
    p.bisectMaxIter = 60;

    p.bGuessProton = 6.0 * p.M;

    p.enableAnimation = true;
    p.animationRealtime = true;
    p.animationSpeed = 500.0;    % Affine-time speed-up factor (larger is faster)
    p.animationFps = 35;
    p.animationMaxFrames = 420;
    p.animationTrailFrames = 160; % <= 0 keeps full trajectory trail
    p.animationLineWidth = 2.2;
    p.animationMarkerSize = 7.5;
    p.animationLoop = true;       % Loop animations continuously (single-window mode)
    p.animationLoopDelay = 0.0;   % Pause between cycles in seconds

    p.staticLineWidth = 2.4;
    p.staticMarkerSize = 8.5;
    p.showTrajectoryIdLabels = true;

    p.enableAnimation3D = true;
    p.animation3DRealtime = true;
    p.animation3DSpeed = p.animationSpeed;
    p.animation3DFps = 30;
    p.animation3DMaxFrames = 380;
    p.animation3DTrailFrames = 150;
    p.animation3DLineWidth = 2.1;
    p.animation3DMarkerSize = 7.5;
    p.animation3DAzimuth0 = 34;
    p.animation3DElevation = 24;
    p.animation3DAutoRotateDegPerSec = 26;

    p.enable3D = true;
    p.showConvergence = false;    % Optional convergence plot (off by default)

    p.beamNumRays2D = 18;
    p.beamNumRays3D = 12;
    p.beamBMinFactor = 0.50;      % Relative to b_crit
    p.beamBMaxFactor = 2.70;      % Relative to b_crit
    p.beamXStart = 42 * p.M;
    p.beamXMin = -11 * p.M;
    p.beamYMargin = 1.08;

    p.beamBgColor = [0.0, 0.0, 0.0];
    p.beamEscapedColor = [0.92, 0.83, 0.65];
    p.beamCapturedColor = [0.72, 0.61, 0.45];
    p.beamCriticalRingColor = [0.60, 0.60, 0.65];
    p.beamHorizonEdgeColor = [0.95, 0.95, 0.95];
    p.beamLineWidth = 2.0;
    p.beamShowArrow = true;
    p.beamShowLabel = true;

    p.beamAnimationFrames = 340;
    p.beamAnimationTrailFrames = -1; % <=0 means full drawn path
    p.beamAnimationLoop = true;
    p.beamAnimationLoopDelay = 0.0;

    p.beam3DLayers = 6;
    p.beam3DZSpread = 3.5 * p.M;
    p.beam3DViewAz = 38;
    p.beam3DViewEl = 18;

    p.interactiveXLim = [-16, 36] * p.M;
    p.interactiveYLim = [-16, 16] * p.M;
    p.interactiveSourceMode = 'click'; % 'click' (default), 'line', or 'point'
    p.interactiveSourceX = 30 * p.M;
    p.interactiveSourceCount = 11;
    p.interactiveSourceHalfSpan = 7.0 * p.M;
    p.interactiveStartPos = [p.interactiveSourceX, 0.0]; % used when sourceMode='point' or 'click'
    p.interactivePickupRadius = 1.6 * p.M;
    p.interactiveSourceMarkerSize = 6.0;
    p.interactiveActiveMarkerSize = 10.5;
    p.interactiveMinDrag = 0.25 * p.M;
    p.interactiveLineWidth = 2.6;
    p.interactiveGlowWidth = 6.0;
    p.interactiveGlowMix = 0.72; % 0 -> same as main color, 1 -> closer to background
    p.interactiveHeadSize = 8.0;
    p.interactiveHistoryLimit = 40;
    p.interactiveFps = 48;
    p.interactiveMaxFrames = 220;
    p.interactiveUsePathSmoothing = true;
    p.interactiveSmoothSamples = 1400;
    p.interactiveUseEasing = false;
    p.interactiveSolverMaxStep = 0.20 * p.M;
    p.interactiveFramePause = true;
    p.interactiveShortDragAutoLeft = true;
    p.interactiveEqualizePlaybackSpeed = true;
    p.interactiveLaunchDurationSec = 2.2;
    p.interactiveSnapToGrid = true;
    p.interactiveGridStep = 1.0 * p.M;
    p.interactiveGridMajorEvery = 2;
    p.interactiveAnimateLaunch = true;
end


function [kappa, label] = particleTypeToKappa(particleType)
    switch lower(strtrim(particleType))
        case 'photon'
            kappa = 0;
            label = 'Photon';
        case 'proton'
            kappa = 1;
            label = 'Proton';
        otherwise
            error('Unknown particleType: %s. Use ''photon'' or ''proton''.', particleType);
    end
end


function L = impactToAngularMomentum(p, b)
    if p.kappa == 0
        L = p.E * b;
    else
        L = b * sqrt(p.E^2 - 1);
    end
end


function status = classifyImpactParameter(p, b)
    traj = integrateTrajectory(p, b);
    status = traj.status;
end


function traj = integrateTrajectory(p, b)
    L = impactToAngularMomentum(p, b);

    V0 = (1 - 2 * p.M / p.r0) * (p.kappa + L^2 / p.r0^2);
    vr0sq = p.E^2 - V0;

    if vr0sq <= 0
        error('Invalid initial state: imaginary radial speed for b = %.6f.', b);
    end

    y0 = [p.r0; -sqrt(vr0sq); pi];

    opts = odeset(...
        'RelTol', p.relTol, ...
        'AbsTol', p.absTol, ...
        'MaxStep', p.maxStep, ...
        'Events', @(lambda, y) geodesicEvents(lambda, y, p));

    [lambda, y, ~, ~, ie] = ode45(@(lambda, y) geodesicODE(lambda, y, p, L), ...
        [0, p.lambdaMax], y0, opts);

    r = y(:, 1);
    vr = y(:, 2);
    phi = y(:, 3);

    traj.lambda = lambda;
    traj.r = r;
    traj.vr = vr;
    traj.phi = phi;
    traj.x = r .* cos(phi);
    traj.y = r .* sin(phi);
    traj.b = b;
    traj.L = L;
    traj.status = classifyByEvents(p, r, vr, ie);
end


function dy = geodesicODE(~, y, p, L)
    r = max(y(1), p.r_horizon * (1 + 1e-8));
    vr = y(2);

    dy = zeros(3, 1);
    dy(1) = vr;
    dy(2) = -(p.M * p.kappa) / r^2 + (L^2) / r^3 - (3 * p.M * L^2) / r^4;
    dy(3) = L / r^2;
end


function [value, isterminal, direction] = geodesicEvents(~, y, p)
    r = y(1);
    vr = y(2);

    value = zeros(3, 1);
    isterminal = [1; 1; 1];
    direction = [-1; 1; 1];

    value(1) = r - p.r_horizon * (1 + 1e-5);

    if vr > 0
        value(2) = r - p.r_escape;
    else
        value(2) = 1;
    end

    value(3) = r - p.r_max;
end


function status = classifyByEvents(p, r, vr, ie)
    if any(ie == 1)
        status = 'captured';
        return;
    end

    if any(ie == 2) || any(ie == 3)
        status = 'escaped';
        return;
    end

    if r(end) <= p.r_horizon * 1.02
        status = 'captured';
    elseif vr(end) > 0 && r(end) > 4 * p.M
        status = 'escaped';
    else
        status = 'captured';
    end
end


function [bLow, bHigh] = autoBracketCritical(p, bGuess)
    bLow = 0.5 * bGuess;
    bHigh = 1.5 * bGuess;

    for k = 1:35
        sLow = classifyImpactParameter(p, bLow);
        sHigh = classifyImpactParameter(p, bHigh);

        if strcmp(sLow, 'captured') && strcmp(sHigh, 'escaped')
            return;
        end

        if ~strcmp(sLow, 'captured')
            bLow = 0.7 * bLow;
        end

        if ~strcmp(sHigh, 'escaped')
            bHigh = 1.3 * bHigh;
        end
    end

    error('Could not find [captured, escaped] bracket for bisection.');
end


function [bCrit, logTable] = bisectCriticalImpact(p, bLow, bHigh, tol, maxIter)
    logTable = zeros(maxIter, 6);

    for iter = 1:maxIter
        bMid = 0.5 * (bLow + bHigh);
        sMid = classifyImpactParameter(p, bMid);

        if strcmp(sMid, 'escaped')
            bHigh = bMid;
            escapedFlag = 1;
        else
            bLow = bMid;
            escapedFlag = 0;
        end

        width = abs(bHigh - bLow);
        logTable(iter, :) = [iter, bLow, bMid, bHigh, width, escapedFlag];

        if width < tol
            logTable = logTable(1:iter, :);
            bCrit = 0.5 * (bLow + bHigh);
            return;
        end
    end

    logTable = logTable(1:maxIter, :);
    bCrit = 0.5 * (bLow + bHigh);
end


function beamTrajectories = generateBeamTrajectories(p, bCrit, nRaysPerLayer, zOffsets)
    if nargin < 4 || isempty(zOffsets)
        zOffsets = 0.0;
    end

    bVals = linspace(p.beamBMinFactor * bCrit, p.beamBMaxFactor * bCrit, nRaysPerLayer);
    nZ = numel(zOffsets);
    beamTrajectories = cell(nRaysPerLayer * nZ, 1);

    idx = 0;
    for iz = 1:nZ
        for ib = 1:nRaysPerLayer
            idx = idx + 1;
            t = integrateBeamPhotonRay(p, bVals(ib));
            t.z = zOffsets(iz) * ones(size(t.x));
            t.beamImpact = bVals(ib);
            t.layer = iz;
            beamTrajectories{idx} = t;
        end
    end
end


function traj = integrateBeamPhotonRay(p, beamImpact)
    L = p.E * beamImpact;

    x0 = p.beamXStart;
    y0cart = beamImpact;
    r0 = hypot(x0, y0cart);
    phi0 = atan2(y0cart, x0);

    V0 = (1 - 2 * p.M / r0) * (L^2 / r0^2);
    vr0sq = p.E^2 - V0;
    vr0sq = max(vr0sq, eps);

    y0 = [r0; -sqrt(vr0sq); phi0];

    opts = odeset(...
        'RelTol', p.relTol, ...
        'AbsTol', p.absTol, ...
        'MaxStep', p.maxStep, ...
        'Events', @(lambda, y) geodesicEvents(lambda, y, p));

    [lambda, y, ~, ~, ie] = ode45(@(lambda, y) geodesicODE(lambda, y, p, L), ...
        [0, p.lambdaMax], y0, opts);

    r = y(:, 1);
    vr = y(:, 2);
    phi = y(:, 3);

    traj.lambda = lambda;
    traj.r = r;
    traj.vr = vr;
    traj.phi = phi;
    traj.x = r .* cos(phi);
    traj.y = r .* sin(phi);
    traj.b = beamImpact;
    traj.L = L;
    traj.status = classifyByEvents(p, r, vr, ie);
end


function plotBeamStyle2D(p, beamTrajectories, particleLabel, bCrit)
    fig = figure('Color', p.beamBgColor, 'Name', 'Photon Beam Lensing 2D', 'NumberTitle', 'off');
    ax = axes('Parent', fig);
    hold(ax, 'on');

    applyBeamAxesStyle2D(ax, p, bCrit, particleLabel, false);
    drawBeamRings2D(ax, p, bCrit);

    for i = 1:numel(beamTrajectories)
        t = beamTrajectories{i};
        [x, y] = clipBeamTrajectoryXY(t, p);
        plot(ax, x, y, 'Color', beamColorFromStatus(p, t.status), 'LineWidth', p.beamLineWidth);
    end

    drawBeamArrowLabel2D(ax, p, bCrit);
end


function animateBeamStyle2D(p, beamTrajectories, particleLabel, bCrit)
    fig = figure('Color', p.beamBgColor, 'Name', 'Photon Beam Lensing 2D Animation', 'NumberTitle', 'off');
    ax = axes('Parent', fig);
    hold(ax, 'on');

    applyBeamAxesStyle2D(ax, p, bCrit, particleLabel, true);
    drawBeamRings2D(ax, p, bCrit);

    nTraj = numel(beamTrajectories);
    lambdaEnd = zeros(nTraj, 1);
    for i = 1:nTraj
        lambdaEnd(i) = beamTrajectories{i}.lambda(end);
    end

    nFrames = max(90, p.beamAnimationFrames);
    lambdaFrames = linspace(0, max(lambdaEnd), nFrames).';
    [xFrame, yFrame, activeFrame] = precomputeAnimationFrames2D(p, beamTrajectories, lambdaFrames);

    rayHandles = zeros(nTraj, 1);
    for i = 1:nTraj
        t = beamTrajectories{i};
        rayHandles(i) = plot(ax, NaN, NaN, ...
            'Color', beamColorFromStatus(p, t.status), ...
            'LineWidth', p.beamLineWidth);
    end

    drawBeamArrowLabel2D(ax, p, bCrit);

    yLim = get(ax, 'YLim');
    yText = yLim(2) - 0.08 * (yLim(2) - yLim(1));
    infoText = text(ax, p.beamXMin + 0.8 * p.M, yText, '', ...
        'Color', [0.95, 0.95, 0.95], ...
        'FontName', 'Consolas', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    keepRunning = true;
    while keepRunning && isFigureVisible(fig)
        for f = 1:nFrames
            if ~isFigureVisible(fig) || ~ishghandle(ax)
                return;
            end

            for i = 1:nTraj
                if p.beamAnimationTrailFrames > 0
                    startIdx = max(1, f - p.beamAnimationTrailFrames + 1);
                else
                    startIdx = 1;
                end

                xSeg = xFrame(startIdx:f, i);
                ySeg = yFrame(startIdx:f, i);
                inWindow = (xSeg >= p.beamXMin) & (xSeg <= p.beamXStart * 1.02);
                xSeg(~inWindow) = NaN;
                ySeg(~inWindow) = NaN;

                if activeFrame(f, i)
                    set(rayHandles(i), 'XData', xSeg, 'YData', ySeg, 'Visible', 'on');
                else
                    set(rayHandles(i), 'Visible', 'off');
                end
            end

            set(infoText, 'String', sprintf('frame %d/%d', f, nFrames));
            drawnow;
        end

        if p.beamAnimationLoop
            if p.beamAnimationLoopDelay > 0 && isFigureVisible(fig)
                pause(p.beamAnimationLoopDelay);
            end
        else
            keepRunning = false;
        end
    end
end


function plotBeamStyle3D(p, beamTrajectories, particleLabel, bCrit)
    fig = figure('Color', p.beamBgColor, 'Name', 'Photon Beam Lensing 3D', 'NumberTitle', 'off');
    ax = axes('Parent', fig);
    hold(ax, 'on');

    applyBeamAxesStyle3D(ax, p, bCrit, particleLabel, false);
    drawBeamRings3D(ax, p, bCrit);

    for i = 1:numel(beamTrajectories)
        t = beamTrajectories{i};
        [x, y] = clipBeamTrajectoryXY(t, p);
        z = t.z(1) * ones(size(x));
        plot3(ax, x, y, z, 'Color', beamColorFromStatus(p, t.status), 'LineWidth', 1.4);
    end

    view(ax, p.beam3DViewAz, p.beam3DViewEl);
end


function animateBeamStyle3D(p, beamTrajectories, particleLabel, bCrit)
    fig = figure('Color', p.beamBgColor, 'Name', 'Photon Beam Lensing 3D Animation', 'NumberTitle', 'off');
    ax = axes('Parent', fig);
    hold(ax, 'on');

    applyBeamAxesStyle3D(ax, p, bCrit, particleLabel, true);
    drawBeamRings3D(ax, p, bCrit);

    nTraj = numel(beamTrajectories);
    lambdaEnd = zeros(nTraj, 1);
    zFrame = zeros(max(90, p.beamAnimationFrames), nTraj);
    for i = 1:nTraj
        lambdaEnd(i) = beamTrajectories{i}.lambda(end);
        zFrame(:, i) = beamTrajectories{i}.z(1);
    end

    nFrames = max(90, p.beamAnimationFrames);
    lambdaFrames = linspace(0, max(lambdaEnd), nFrames).';
    [xFrame, yFrame, activeFrame] = precomputeAnimationFrames2D(p, beamTrajectories, lambdaFrames);

    rayHandles = zeros(nTraj, 1);
    for i = 1:nTraj
        t = beamTrajectories{i};
        rayHandles(i) = plot3(ax, NaN, NaN, NaN, ...
            'Color', beamColorFromStatus(p, t.status), ...
            'LineWidth', 1.3);
    end

    infoText = text(ax, p.beamXMin, 0.85 * p.beamBMaxFactor * bCrit, 0.88 * p.beam3DZSpread, '', ...
        'Color', [0.95, 0.95, 0.95], ...
        'FontName', 'Consolas', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    keepRunning = true;
    while keepRunning && isFigureVisible(fig)
        for f = 1:nFrames
            if ~isFigureVisible(fig) || ~ishghandle(ax)
                return;
            end

            for i = 1:nTraj
                if p.beamAnimationTrailFrames > 0
                    startIdx = max(1, f - p.beamAnimationTrailFrames + 1);
                else
                    startIdx = 1;
                end

                xSeg = xFrame(startIdx:f, i);
                ySeg = yFrame(startIdx:f, i);
                zSeg = zFrame(startIdx:f, i);
                inWindow = (xSeg >= p.beamXMin) & (xSeg <= p.beamXStart * 1.02);
                xSeg(~inWindow) = NaN;
                ySeg(~inWindow) = NaN;
                zSeg(~inWindow) = NaN;

                if activeFrame(f, i)
                    set(rayHandles(i), 'XData', xSeg, 'YData', ySeg, 'ZData', zSeg, 'Visible', 'on');
                else
                    set(rayHandles(i), 'Visible', 'off');
                end
            end

            view(ax, p.beam3DViewAz + 18 * (f - 1) / max(nFrames - 1, 1), p.beam3DViewEl);
            set(infoText, 'String', sprintf('frame %d/%d', f, nFrames));
            drawnow;
        end

        if p.beamAnimationLoop
            if p.beamAnimationLoopDelay > 0 && isFigureVisible(fig)
                pause(p.beamAnimationLoopDelay);
            end
        else
            keepRunning = false;
        end
    end
end


function applyBeamAxesStyle2D(ax, p, bCrit, particleLabel, isAnimation)
    rs = 2 * p.M;
    yMax = p.beamBMaxFactor * bCrit * p.beamYMargin;

    set(ax, 'Color', p.beamBgColor, 'XColor', [0.88, 0.88, 0.88], 'YColor', [0.88, 0.88, 0.88]);
    axis(ax, 'equal');
    xlim(ax, [p.beamXMin, p.beamXStart * 1.02]);
    ylim(ax, [-1.4 * rs, yMax]);
    axis(ax, 'off');

    if isAnimation
        title(ax, sprintf('%s beam lensing animation (2D)', particleLabel), 'Color', [0.95, 0.95, 0.95]);
    else
        title(ax, sprintf('%s beam lensing trajectories (2D)', particleLabel), 'Color', [0.95, 0.95, 0.95]);
    end
end


function applyBeamAxesStyle3D(ax, p, bCrit, particleLabel, isAnimation)
    yMax = p.beamBMaxFactor * bCrit * p.beamYMargin;
    set(ax, 'Color', p.beamBgColor, 'XColor', [0.88, 0.88, 0.88], 'YColor', [0.88, 0.88, 0.88], 'ZColor', [0.88, 0.88, 0.88]);
    axis(ax, 'equal');
    xlim(ax, [p.beamXMin, p.beamXStart * 1.02]);
    ylim(ax, [-1.2 * yMax, yMax]);
    zlim(ax, [-1.2 * p.beam3DZSpread, 1.2 * p.beam3DZSpread]);
    axis(ax, 'off');

    if isAnimation
        title(ax, sprintf('%s beam lensing animation (3D)', particleLabel), 'Color', [0.95, 0.95, 0.95]);
    else
        title(ax, sprintf('%s beam lensing trajectories (3D)', particleLabel), 'Color', [0.95, 0.95, 0.95]);
    end
end


function drawBeamRings2D(ax, p, bCrit)
    th = linspace(0, 2 * pi, 500);
    rs = 2 * p.M;
    fill(ax, rs * cos(th), rs * sin(th), p.beamBgColor, ...
        'EdgeColor', p.beamHorizonEdgeColor, 'LineWidth', 2.8);
    plot(ax, bCrit * cos(th), bCrit * sin(th), ...
        'Color', p.beamCriticalRingColor, 'LineWidth', 2.2);
end


function drawBeamRings3D(ax, p, bCrit)
    [sx, sy, sz] = sphere(80);
    rs = 2 * p.M;
    surf(ax, rs * sx, rs * sy, rs * sz, ...
        'FaceColor', [0.03, 0.03, 0.03], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 1.0);

    surf(ax, bCrit * sx, bCrit * sy, bCrit * sz, ...
        'FaceColor', p.beamCriticalRingColor, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.10);
end


function drawBeamArrowLabel2D(ax, p, bCrit)
    if ~p.beamShowArrow && ~p.beamShowLabel
        return;
    end

    rs = 2 * p.M;
    xArrow = p.beamXStart * 0.83;
    arrowColor = [0.95, 0.95, 0.95];

    if p.beamShowArrow
        quiver(ax, xArrow, 0, 0, bCrit, 0, ...
            'Color', arrowColor, 'LineWidth', 4.0, 'MaxHeadSize', 0.35);
    end

    if p.beamShowLabel
        text(ax, xArrow - 0.28 * bCrit, 0.55 * bCrit, sprintf('%.1f r_s', bCrit / rs), ...
            'Color', arrowColor, ...
            'FontSize', 20, ...
            'FontWeight', 'bold');
    end
end


function [x, y] = clipBeamTrajectoryXY(traj, p)
    mask = (traj.x >= p.beamXMin) & (traj.x <= p.beamXStart * 1.02);
    x = traj.x(mask);
    y = traj.y(mask);
    if isempty(x)
        x = NaN;
        y = NaN;
    end
end


function c = beamColorFromStatus(p, status)
    if strcmp(status, 'captured')
        c = p.beamCapturedColor;
    else
        c = p.beamEscapedColor;
    end
end


function plotTrajectories2D(p, trajectories, styles, particleLabel, bCritNum, bAnalytic, axIn)
    [ax, ~] = prepareAxes(axIn, 'Geodesic 2D');
    hold(ax, 'on');

    th = linspace(0, 2 * pi, 500);
    fill(ax, p.r_horizon * cos(th), p.r_horizon * sin(th), [0.08, 0.08, 0.08], ...
        'FaceAlpha', 0.95, 'EdgeColor', 'none');

    plot(ax, 3 * p.M * cos(th), 3 * p.M * sin(th), 'k--', 'LineWidth', 1.1);

    for i = 1:numel(trajectories)
        t = trajectories{i};
        st = styles(i);
        plot(ax, t.x, t.y, ...
            'LineStyle', st.lineStyle, ...
            'LineWidth', p.staticLineWidth, ...
            'Color', st.color);

        if strcmp(st.marker, 'o')
            plot(ax, t.x(end), t.y(end), 'o', ...
                'MarkerSize', p.staticMarkerSize, ...
                'MarkerFaceColor', st.color, ...
                'MarkerEdgeColor', 'k', ...
                'LineWidth', 1.0);
        else
            plot(ax, t.x(end), t.y(end), 'x', ...
                'MarkerSize', p.staticMarkerSize, ...
                'Color', st.color, ...
                'LineWidth', 1.8);
        end

        if p.showTrajectoryIdLabels
            text(ax, t.x(end), t.y(end), sprintf(' %d', i), ...
                'Color', st.color, ...
                'FontWeight', 'bold', ...
                'FontSize', 10, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle');
        end
    end

    axis(ax, 'equal');
    rPlot = max([p.r0, p.r_escape]) * 1.05;
    xlim(ax, [-rPlot, rPlot]);
    ylim(ax, [-rPlot, rPlot]);
    grid(ax, 'on');
    xlabel(ax, 'x (M)');
    ylabel(ax, 'y (M)');
    set(ax, 'FontSize', 11, 'LineWidth', 1.0);

    if isnan(bAnalytic)
        title(ax, sprintf('%s trajectories, b_{crit,num}=%.6f', particleLabel, bCritNum));
    else
        title(ax, sprintf('%s trajectories, b_{crit,num}=%.6f, b_{crit,ana}=%.6f', ...
            particleLabel, bCritNum, bAnalytic));
    end

    if ~p.useSingleWindow
        legend(ax, 'Location', 'northeast', 'Interpreter', 'none');
    end
end


function plotTrajectories3D(p, trajectories, styles, particleLabel, axIn)
    [ax, ~] = prepareAxes(axIn, 'Geodesic 3D');
    hold(ax, 'on');

    [sx, sy, sz] = sphere(80);
    surf(ax, p.r_horizon * sx, p.r_horizon * sy, p.r_horizon * sz, ...
        'FaceColor', [0.08, 0.08, 0.08], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.85);

    for i = 1:numel(trajectories)
        t = trajectories{i};
        st = styles(i);
        plot3(ax, t.x, t.y, zeros(size(t.x)), ...
            'LineStyle', st.lineStyle, ...
            'LineWidth', p.staticLineWidth, ...
            'Color', st.color);
    end

    axis(ax, 'equal');
    rPlot = max([p.r0, p.r_escape]) * 1.05;
    xlim(ax, [-rPlot, rPlot]);
    ylim(ax, [-rPlot, rPlot]);
    zlim(ax, [-0.35 * p.r0, 0.35 * p.r0]);

    xlabel(ax, 'x (M)');
    ylabel(ax, 'y (M)');
    zlabel(ax, 'z (M)');
    title(ax, sprintf('%s trajectories in 3D view (equatorial plane)', particleLabel));
    grid(ax, 'on');
    set(ax, 'FontSize', 11, 'LineWidth', 1.0);
    view(ax, 38, 24);
end


function animateTrajectories3D(p, trajectories, styles, particleLabel, bCritNum, bAnalytic, axIn)
    [ax, fig] = prepareAxes(axIn, 'Affine-time animation 3D');
    hold(ax, 'on');

    [sx, sy, sz] = sphere(80);
    surf(ax, p.r_horizon * sx, p.r_horizon * sy, p.r_horizon * sz, ...
        'FaceColor', [0.08, 0.08, 0.08], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.85);

    nTraj = numel(trajectories);
    lambdaEnd = zeros(nTraj, 1);
    for i = 1:nTraj
        lambdaEnd(i) = trajectories{i}.lambda(end);
    end

    lambdaMax = max(lambdaEnd);
    if p.animation3DRealtime
        totalWallTime = lambdaMax / max(p.animation3DSpeed, eps);
        nFrames = ceil(totalWallTime * p.animation3DFps);
    else
        nFrames = p.animation3DMaxFrames;
    end

    nFrames = max(90, min(p.animation3DMaxFrames, nFrames));
    lambdaFrames = linspace(0, lambdaMax, nFrames).';

    [xFrame, yFrame, zFrame, activeFrame] = precomputeAnimationFrames3D(p, trajectories, lambdaFrames);

    trailHandles = zeros(nTraj, 1);
    markerHandles = zeros(nTraj, 1);

    for i = 1:nTraj
        st = styles(i);

        trailHandles(i) = plot3(ax, NaN, NaN, NaN, '-', ...
            'LineWidth', p.animation3DLineWidth, ...
            'LineStyle', st.lineStyle, ...
            'Color', st.color);

        if strcmp(st.marker, 'o')
            markerHandles(i) = plot3(ax, NaN, NaN, NaN, 'o', ...
                'MarkerSize', p.animation3DMarkerSize, ...
                'MarkerFaceColor', st.color, ...
                'MarkerEdgeColor', 'k', ...
                'LineWidth', 1.0, ...
                'Color', st.color);
        else
            markerHandles(i) = plot3(ax, NaN, NaN, NaN, 'x', ...
                'MarkerSize', p.animation3DMarkerSize, ...
                'LineWidth', 1.6, ...
                'Color', st.color);
        end
    end

    axis(ax, 'equal');
    rPlot = max([p.r0, p.r_escape]) * 1.05;
    xlim(ax, [-rPlot, rPlot]);
    ylim(ax, [-rPlot, rPlot]);
    zlim(ax, [-0.35 * p.r0, 0.35 * p.r0]);

    xlabel(ax, 'x (M)');
    ylabel(ax, 'y (M)');
    zlabel(ax, 'z (M)');
    set(ax, 'FontSize', 11, 'LineWidth', 1.0);

    if isnan(bAnalytic)
        title(ax, sprintf('%s affine-time animation 3D, b_{crit,num}=%.6f', particleLabel, bCritNum));
    else
        title(ax, sprintf('%s affine-time animation 3D, b_{crit,num}=%.6f, b_{crit,ana}=%.6f', ...
            particleLabel, bCritNum, bAnalytic));
    end

    grid(ax, 'on');
    if ~p.useSingleWindow
        legend(ax, 'Location', 'northeast', 'Interpreter', 'none');
    end
    view(ax, p.animation3DAzimuth0, p.animation3DElevation);

    infoText = text(ax, -0.98 * rPlot, 0.95 * rPlot, 0.32 * rPlot, '', ...
        'FontName', 'Consolas', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    doRealtime = p.animation3DRealtime && isFigureVisible(fig);
    tRef = tic;

    for f = 1:nFrames
        if ~ishghandle(ax) || ~isFigureVisible(fig)
            return;
        end

        for i = 1:nTraj
            if p.animation3DTrailFrames > 0
                startIdx = max(1, f - p.animation3DTrailFrames + 1);
            else
                startIdx = 1;
            end

            set(trailHandles(i), ...
                'XData', xFrame(startIdx:f, i), ...
                'YData', yFrame(startIdx:f, i), ...
                'ZData', zFrame(startIdx:f, i));

            if activeFrame(f, i)
                set(markerHandles(i), ...
                    'XData', xFrame(f, i), ...
                    'YData', yFrame(f, i), ...
                    'ZData', 0, ...
                    'Visible', 'on');
            else
                set(markerHandles(i), 'Visible', 'off');
            end
        end

        if doRealtime
            camTime = lambdaFrames(f) / max(p.animation3DSpeed, eps);
            waitTime = camTime - toc(tRef);
            if waitTime > 0
                pause(waitTime);
            end
        else
            camTime = (f - 1) / max(p.animation3DFps, eps);
        end

        az = p.animation3DAzimuth0 + p.animation3DAutoRotateDegPerSec * camTime;
        view(ax, az, p.animation3DElevation);

        set(infoText, 'String', sprintf('lambda = %.3f M | frame %d/%d | az = %.1f deg', ...
            lambdaFrames(f) / p.M, f, nFrames, az));

        drawnow;
    end
end


function animateTrajectories2D(p, trajectories, styles, particleLabel, bCritNum, bAnalytic, axIn)
    [ax, fig] = prepareAxes(axIn, 'Affine-time animation 2D');
    hold(ax, 'on');

    th = linspace(0, 2 * pi, 500);
    fill(ax, p.r_horizon * cos(th), p.r_horizon * sin(th), [0.08, 0.08, 0.08], ...
        'FaceAlpha', 0.95, 'EdgeColor', 'none');

    plot(ax, 3 * p.M * cos(th), 3 * p.M * sin(th), 'k--', 'LineWidth', 1.1);

    axis(ax, 'equal');
    rPlot = max([p.r0, p.r_escape]) * 1.05;
    xlim(ax, [-rPlot, rPlot]);
    ylim(ax, [-rPlot, rPlot]);
    grid(ax, 'on');
    xlabel(ax, 'x (M)');
    ylabel(ax, 'y (M)');
    set(ax, 'FontSize', 11, 'LineWidth', 1.0);

    if isnan(bAnalytic)
        title(ax, sprintf('%s affine-time animation, b_{crit,num}=%.6f', particleLabel, bCritNum));
    else
        title(ax, sprintf('%s affine-time animation, b_{crit,num}=%.6f, b_{crit,ana}=%.6f', ...
            particleLabel, bCritNum, bAnalytic));
    end

    nTraj = numel(trajectories);
    lambdaEnd = zeros(nTraj, 1);
    for i = 1:nTraj
        lambdaEnd(i) = trajectories{i}.lambda(end);
    end

    lambdaMax = max(lambdaEnd);
    if p.animationRealtime
        totalWallTime = lambdaMax / max(p.animationSpeed, eps);
        nFrames = ceil(totalWallTime * p.animationFps);
    else
        nFrames = p.animationMaxFrames;
    end

    nFrames = max(90, min(p.animationMaxFrames, nFrames));
    lambdaFrames = linspace(0, lambdaMax, nFrames).';

    [xFrame, yFrame, activeFrame] = precomputeAnimationFrames2D(p, trajectories, lambdaFrames);

    trailHandles = zeros(nTraj, 1);
    markerHandles = zeros(nTraj, 1);

    for i = 1:nTraj
        st = styles(i);

        trailHandles(i) = plot(ax, NaN, NaN, '-', ...
            'LineWidth', p.animationLineWidth, ...
            'LineStyle', st.lineStyle, ...
            'Color', st.color);

        if strcmp(st.marker, 'o')
            markerHandles(i) = plot(ax, NaN, NaN, 'o', ...
                'MarkerSize', p.animationMarkerSize, ...
                'MarkerFaceColor', st.color, ...
                'MarkerEdgeColor', 'k', ...
                'LineWidth', 1.0, ...
                'Color', st.color);
        else
            markerHandles(i) = plot(ax, NaN, NaN, 'x', ...
                'MarkerSize', p.animationMarkerSize, ...
                'LineWidth', 1.6, ...
                'Color', st.color);
        end
    end

    if ~p.useSingleWindow
        legend(ax, 'Location', 'northeast', 'Interpreter', 'none');
    end

    infoText = text(ax, -0.98 * rPlot, 0.95 * rPlot, '', ...
        'FontName', 'Consolas', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    doRealtime = p.animationRealtime && isFigureVisible(fig);
    tRef = tic;

    for f = 1:nFrames
        if ~ishghandle(ax) || ~isFigureVisible(fig)
            return;
        end

        for i = 1:nTraj
            if p.animationTrailFrames > 0
                startIdx = max(1, f - p.animationTrailFrames + 1);
            else
                startIdx = 1;
            end

            set(trailHandles(i), 'XData', xFrame(startIdx:f, i), 'YData', yFrame(startIdx:f, i));

            if activeFrame(f, i)
                set(markerHandles(i), 'XData', xFrame(f, i), 'YData', yFrame(f, i), 'Visible', 'on');
            else
                set(markerHandles(i), 'Visible', 'off');
            end
        end

        set(infoText, 'String', sprintf('lambda = %.3f M | frame %d/%d', ...
            lambdaFrames(f) / p.M, f, nFrames));

        if doRealtime
            targetTime = lambdaFrames(f) / max(p.animationSpeed, eps);
            waitTime = targetTime - toc(tRef);
            if waitTime > 0
                pause(waitTime);
            end
        end

        drawnow;
    end
end


function tf = isFigureVisible(fig)
    tf = false;
    if ishghandle(fig)
        tf = strcmpi(get(fig, 'Visible'), 'on');
    end
end


function plotBisectionConvergence(logTable, particleLabel, axIn)
    [ax, ~] = prepareAxes(axIn, 'Bisection convergence');
    semilogy(ax, logTable(:, 1), logTable(:, 5), 'o-', 'LineWidth', 1.4, 'MarkerSize', 5);
    grid(ax, 'on');
    xlabel(ax, 'Iteration');
    ylabel(ax, '|b_{high} - b_{low}|');
    title(ax, sprintf('Bisection convergence for %s critical impact parameter', particleLabel));
end


function runAnimationPass(p, trajectories, styles, particleLabel, bCritNum, bAnalytic, ax2DAnim, ax3DAnim)
    if p.enableAnimation
        animateTrajectories2D(p, trajectories, styles, particleLabel, bCritNum, bAnalytic, ax2DAnim);
    elseif p.useSingleWindow
        showDisabledPanel(ax2DAnim, '2D Affine-Time Animation', '2D animation is disabled');
    end

    if p.enableAnimation3D
        animateTrajectories3D(p, trajectories, styles, particleLabel, bCritNum, bAnalytic, ax3DAnim);
    elseif p.useSingleWindow
        showDisabledPanel(ax3DAnim, '3D Affine-Time Animation', '3D animation is disabled');
    end
end


function dashboard = createDashboardFigure(p)
    screenRect = get(0, 'ScreenSize');
    scale = max(0.75, min(1.35, p.dashboardScale));
    figW = round(screenRect(3) * 0.95 * scale);
    figH = round(screenRect(4) * 0.90 * scale);
    figW = min(max(figW, 1200), max(920, screenRect(3) - 20));
    figH = min(max(figH, 760), max(620, screenRect(4) - 60));
    figX = max(5, round((screenRect(3) - figW) / 2));
    figY = max(30, round((screenRect(4) - figH) / 2));

    dashboard.fig = figure('Color', 'w', 'Name', 'Geodesic Dashboard', ...
        'NumberTitle', 'off', 'Position', [figX, figY, figW, figH]);

    dashboard.ax2DStatic = subplot(2, 3, 1);
    dashboard.ax2DAnim = subplot(2, 3, 2);
    dashboard.axConv = subplot(2, 3, 3);
    dashboard.ax3DStatic = subplot(2, 3, 4);
    dashboard.ax3DAnim = subplot(2, 3, 5);
    dashboard.axInfo = subplot(2, 3, 6);

    set(dashboard.ax2DStatic, 'Position', [0.04, 0.55, 0.56, 0.40]);
    set(dashboard.ax2DAnim, 'Position', [0.04, 0.07, 0.56, 0.40]);
    set(dashboard.ax3DAnim, 'Position', [0.62, 0.40, 0.36, 0.55]);
    set(dashboard.axInfo, 'Position', [0.62, 0.24, 0.36, 0.14]);
    set(dashboard.ax3DStatic, 'Position', [0.62, 0.07, 0.18, 0.14]);
    set(dashboard.axConv, 'Position', [0.81, 0.07, 0.17, 0.14]);

    title(dashboard.ax2DStatic, '2D trajectories (larger view)');
    title(dashboard.ax2DAnim, '2D animation (larger view)');
    title(dashboard.axConv, 'Trajectory stats');
    title(dashboard.ax3DStatic, '3D trajectories');
    title(dashboard.ax3DAnim, '3D animation');
end


function populateInfoPanel(ax, p, particleLabel, bCritNum, bAnalytic, trajectories, styles)
    cla(ax);
    axis(ax, 'off');
    hold(ax, 'on');

    if isnan(bAnalytic)
        analyticText = 'N/A';
    else
        analyticText = sprintf('%.6f', bAnalytic);
    end

    infoLines = {
        sprintf('%s Simulation Dashboard', particleLabel)
        sprintf('M = %.3f, E = %.3f', p.M, p.E)
        sprintf('b_{crit,num} = %.6f', bCritNum)
        sprintf('b_{crit,ana} = %s', analyticText)
        'Captured: solid red + x marker'
        'Escaped: dashed blue + o marker'
    };

    text(ax, 0.02, 0.98, strjoin(infoLines, '\n'), ...
        'Units', 'normalized', ...
        'FontName', 'Consolas', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    statusLines = composeStatusListLines(trajectories);
    text(ax, 0.02, 0.60, strjoin(statusLines, '\n'), ...
        'Units', 'normalized', ...
        'FontName', 'Consolas', ...
        'FontSize', 8.5, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    hEvent = plot(ax, NaN, NaN, '-', 'Color', [0.08, 0.08, 0.08], 'LineWidth', 4);
    hPhoton = plot(ax, NaN, NaN, 'k--', 'LineWidth', 1.2);
    hCap = plot(ax, NaN, NaN, '-x', 'Color', [0.85, 0.18, 0.18], ...
        'LineWidth', 2.0, 'MarkerSize', 7);
    hEsc = plot(ax, NaN, NaN, '--o', 'Color', [0.10, 0.45, 0.90], ...
        'LineWidth', 2.0, 'MarkerSize', 6, ...
        'MarkerFaceColor', [0.10, 0.45, 0.90], 'MarkerEdgeColor', 'k');

    lgd = legend(ax, [hEvent, hPhoton, hCap, hEsc], ...
        {'Event horizon', 'Photon sphere r=3M', 'Captured trajectory', 'Escaped trajectory'}, ...
        'Location', 'south', 'Interpreter', 'none');
    set(lgd, 'FontSize', 8, 'Box', 'off');

    if ~isempty(styles)
        title(ax, 'Run Summary and Visual Key');
    else
        title(ax, 'Run Summary');
    end
end


function statusLines = composeStatusListLines(trajectories)
    statusLines = cell(numel(trajectories) + 1, 1);
    statusLines{1} = 'Trajectory IDs (lihat angka di plot 2D):';

    for i = 1:numel(trajectories)
        t = trajectories{i};
        if strcmp(t.status, 'captured')
            flag = 'C';
        else
            flag = 'E';
        end
        statusLines{i + 1} = sprintf('#%d [%s] b=%.5f', i, flag, t.b);
    end
end


function populateTrajectoryStatsPanel(ax, trajectories)
    if isempty(ax) || ~ishghandle(ax)
        return;
    end

    cla(ax);
    axis(ax, 'off');

    nTraj = numel(trajectories);
    bVals = zeros(nTraj, 1);
    nCaptured = 0;
    nEscaped = 0;

    for i = 1:nTraj
        bVals(i) = trajectories{i}.b;
        if strcmp(trajectories{i}.status, 'captured')
            nCaptured = nCaptured + 1;
        else
            nEscaped = nEscaped + 1;
        end
    end

    if nTraj > 0
        capPct = 100 * nCaptured / nTraj;
        escPct = 100 * nEscaped / nTraj;
        bMin = min(bVals);
        bMax = max(bVals);
    else
        capPct = 0;
        escPct = 0;
        bMin = NaN;
        bMax = NaN;
    end

    statLines = {
        'Trajectory Stats'
        sprintf('Total trajectories: %d', nTraj)
        sprintf('Captured: %d (%.1f%%)', nCaptured, capPct)
        sprintf('Escaped:  %d (%.1f%%)', nEscaped, escPct)
        ''
        sprintf('b range: [%.5f, %.5f]', bMin, bMax)
        ''
        'Kode: C = captured, E = escaped'
    };

    text(ax, 0.04, 0.96, strjoin(statLines, '\n'), ...
        'Units', 'normalized', ...
        'FontName', 'Consolas', ...
        'FontSize', 9.5, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');
    title(ax, 'Trajectory Stats');
end


function showDisabledPanel(ax, panelTitle, message)
    if isempty(ax) || ~ishghandle(ax)
        return;
    end

    cla(ax);
    axis(ax, 'off');
    title(ax, panelTitle);
    text(ax, 0.5, 0.5, message, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', 'Consolas', ...
        'FontSize', 10);
end


function styles = buildTrajectoryStyles(trajectories)
    n = numel(trajectories);
    styles = repmat(struct('color', [0.45, 0.45, 0.45], 'lineStyle', '-.', 'marker', 's'), n, 1);

    statusList = cell(n, 1);
    for i = 1:n
        statusList{i} = trajectories{i}.status;
    end

    capIdx = find(strcmp(statusList, 'captured'));
    escIdx = find(strcmp(statusList, 'escaped'));

    for k = 1:numel(capIdx)
        idx = capIdx(k);
        t = (k - 1) / max(numel(capIdx) - 1, 1);
        c0 = [0.95, 0.35, 0.20];
        c1 = [0.58, 0.07, 0.07];
        styles(idx).color = (1 - t) * c0 + t * c1;
        styles(idx).lineStyle = '-';
        styles(idx).marker = 'x';
    end

    for k = 1:numel(escIdx)
        idx = escIdx(k);
        t = (k - 1) / max(numel(escIdx) - 1, 1);
        c0 = [0.25, 0.75, 0.98];
        c1 = [0.04, 0.30, 0.78];
        styles(idx).color = (1 - t) * c0 + t * c1;
        styles(idx).lineStyle = '--';
        styles(idx).marker = 'o';
    end
end


function [xFrame, yFrame, activeFrame] = precomputeAnimationFrames2D(p, trajectories, lambdaFrames)
    [xFrame, yFrame, activeFrame] = precomputeAnimationFramesXY(p, trajectories, lambdaFrames);
end


function [xFrame, yFrame, zFrame, activeFrame] = precomputeAnimationFrames3D(p, trajectories, lambdaFrames)
    [xFrame, yFrame, activeFrame] = precomputeAnimationFramesXY(p, trajectories, lambdaFrames);
    zFrame = zeros(size(xFrame));
end


function [xFrame, yFrame, activeFrame] = precomputeAnimationFramesXY(p, trajectories, lambdaFrames)
    nFrames = numel(lambdaFrames);
    nTraj = numel(trajectories);

    xFrame = zeros(nFrames, nTraj);
    yFrame = zeros(nFrames, nTraj);
    activeFrame = false(nFrames, nTraj);

    for i = 1:nTraj
        t = trajectories{i};
        q = min(lambdaFrames, t.lambda(end));
        xFrame(:, i) = interp1(t.lambda, t.x, q, 'linear');
        yFrame(:, i) = interp1(t.lambda, t.y, q, 'linear');
        activeFrame(:, i) = lambdaFrames <= t.lambda(end);
    end
end


function [ax, fig] = prepareAxes(axIn, figName)
    if nargin < 1
        axIn = [];
    end

    if isempty(axIn) || ~ishghandle(axIn)
        fig = figure('Color', 'w', 'Name', figName, 'NumberTitle', 'off');
        ax = axes('Parent', fig);
    else
        ax = axIn;
        fig = ancestorFigure(ax);
        cla(ax);
    end
end


function fig = ancestorFigure(h)
    fig = [];
    if isempty(h) || ~ishghandle(h)
        return;
    end

    fig = h;
    while ishghandle(fig) && ~strcmpi(get(fig, 'Type'), 'figure')
        fig = get(fig, 'Parent');
    end
end


function runInteractivePhotonSandbox(p, particleLabel, bCritNum, bAnalytic)
    fig = figure('Color', [0.01, 0.01, 0.02], ...
        'Name', 'Interactive Photon Sandbox', ...
        'NumberTitle', 'off');

    axMain = axes('Parent', fig, 'Position', [0.05, 0.08, 0.70, 0.86]);
    axInfo = axes('Parent', fig, 'Position', [0.78, 0.08, 0.20, 0.86]);

    hold(axMain, 'on');
    set(axMain, 'Color', [0.01, 0.01, 0.02], ...
        'XColor', [0.88, 0.88, 0.90], ...
        'YColor', [0.88, 0.88, 0.90], ...
        'FontSize', 11, ...
        'LineWidth', 1.0);
    axis(axMain, 'equal');
    xlim(axMain, p.interactiveXLim);
    ylim(axMain, p.interactiveYLim);
    grid(axMain, 'on');
    set(axMain, 'Layer', 'top');

    gridStep = max(eps, p.interactiveGridStep);
    majorEvery = max(1, round(p.interactiveGridMajorEvery));
    majorStep = majorEvery * gridStep;
    xTickStart = ceil(p.interactiveXLim(1) / majorStep) * majorStep;
    xTickEnd = floor(p.interactiveXLim(2) / majorStep) * majorStep;
    yTickStart = ceil(p.interactiveYLim(1) / majorStep) * majorStep;
    yTickEnd = floor(p.interactiveYLim(2) / majorStep) * majorStep;
    if xTickEnd >= xTickStart
        set(axMain, 'XTick', xTickStart:majorStep:xTickEnd);
    end
    if yTickEnd >= yTickStart
        set(axMain, 'YTick', yTickStart:majorStep:yTickEnd);
    end

    xlabel(axMain, 'x (M)');
    ylabel(axMain, 'y (M)');
    title(axMain, sprintf('%s Trajectories in Curved Spacetime around a Black Hole', particleLabel), ...
        'Color', [0.96, 0.96, 0.98]);

    th = linspace(0, 2 * pi, 500);

    guideRadii = (4:4:16) * p.M;
    for i = 1:numel(guideRadii)
        rg = guideRadii(i);
        plot(axMain, rg * cos(th), rg * sin(th), ':', ...
            'Color', [0.13, 0.13, 0.16], 'LineWidth', 0.9);
    end

    plot(axMain, p.interactiveXLim, [0, 0], ':', ...
        'Color', [0.14, 0.14, 0.18], 'LineWidth', 0.9);
    plot(axMain, [0, 0], p.interactiveYLim, ':', ...
        'Color', [0.14, 0.14, 0.18], 'LineWidth', 0.9);

    fill(axMain, p.r_horizon * cos(th), p.r_horizon * sin(th), [0.00, 0.00, 0.00], ...
        'EdgeColor', [0.97, 0.97, 0.99], 'LineWidth', 2.6);
    plot(axMain, 3 * p.M * cos(th), 3 * p.M * sin(th), '--', ...
        'Color', [0.55, 0.55, 0.62], 'LineWidth', 1.5);
    plot(axMain, bCritNum * cos(th), bCritNum * sin(th), '-', ...
        'Color', [0.75, 0.75, 0.80], 'LineWidth', 1.6);

    rs = 2 * p.M;
    text(axMain, -0.55 * bCritNum, 0.82 * bCritNum, ...
        sprintf('b_{crit} = %.2f M (%.2f r_s)', bCritNum / p.M, bCritNum / rs), ...
        'Color', [0.88, 0.88, 0.92], ...
        'FontName', 'Consolas', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');

    sourceMode = getInteractiveSourceMode(p);
    sourcePoints = buildInteractiveSourcePoints(p);
    nSource = size(sourcePoints, 1);
    activeEmitterIdx = ceil(nSource / 2);

    if nSource > 1
        plot(axMain, sourcePoints(:, 1), sourcePoints(:, 2), '-', ...
            'Color', [0.34, 0.30, 0.18], 'LineWidth', 1.0);
    end

    sourceHandle = plot(axMain, sourcePoints(:, 1), sourcePoints(:, 2), 'o', ...
        'LineStyle', 'none', ...
        'MarkerSize', p.interactiveSourceMarkerSize, ...
        'MarkerFaceColor', [0.99, 0.86, 0.35], ...
        'MarkerEdgeColor', [0.10, 0.10, 0.10], ...
        'LineWidth', 1.0);

    activeEmitterHandle = plot(axMain, ...
        sourcePoints(activeEmitterIdx, 1), sourcePoints(activeEmitterIdx, 2), 'o', ...
        'MarkerSize', p.interactiveActiveMarkerSize, ...
        'MarkerFaceColor', [1.00, 0.95, 0.60], ...
        'MarkerEdgeColor', [0.15, 0.15, 0.15], ...
        'LineWidth', 1.3);

    previewHandle = plot(axMain, NaN, NaN, '--', ...
        'Color', [0.95, 0.95, 0.95], ...
        'LineWidth', 1.5, ...
        'Visible', 'off');

    rayGlowHandle = plot(axMain, NaN, NaN, '-', ...
        'Color', [0.11, 0.30, 0.36], ...
        'LineWidth', p.interactiveGlowWidth, ...
        'Visible', 'off');

    rayHandle = plot(axMain, NaN, NaN, '-', ...
        'Color', [0.28, 0.86, 0.98], ...
        'LineWidth', p.interactiveLineWidth);

    headHandle = plot(axMain, NaN, NaN, 'o', ...
        'Color', [0.28, 0.86, 0.98], ...
        'MarkerSize', p.interactiveHeadSize, ...
        'MarkerFaceColor', [0.28, 0.86, 0.98], ...
        'MarkerEdgeColor', 'k', ...
        'Visible', 'off');

    if strcmp(sourceMode, 'click')
        bannerText = '';
    else
        bannerText = '';
    end

    bannerHandle = text(axMain, p.interactiveXLim(1) + 0.8 * p.M, p.interactiveYLim(2) - 0.9 * p.M, ...
        bannerText, ...
        'Color', [0.96, 0.96, 0.98], ...
        'FontName', 'Consolas', ...
        'FontSize', 11, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    axis(axInfo, 'off');
    hold(axInfo, 'on');

    infoHandle = text(axInfo, 0.03, 0.98, ...
        composeInteractiveInfoText(p, bCritNum, bAnalytic, 'READY', [NaN, NaN], NaN, NaN, NaN), ...
        'Units', 'normalized', ...
        'Color', [0.20, 0.22, 0.25], ...
        'FontName', 'Consolas', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    data = struct();
    data.p = p;
    data.sourceMode = sourceMode;
    data.axMain = axMain;
    data.sourcePoints = sourcePoints;
    data.activeEmitterIdx = activeEmitterIdx;
    data.isDragging = false;
    data.previewHandle = previewHandle;
    data.rayGlowHandle = rayGlowHandle;
    data.rayHandle = rayHandle;
    data.headHandle = headHandle;
    data.sourceHandle = sourceHandle;
    data.activeEmitterHandle = activeEmitterHandle;
    data.bannerHandle = bannerHandle;
    data.bannerBaseText = bannerText;
    data.infoHandle = infoHandle;
    data.bCritNum = bCritNum;
    data.bAnalytic = bAnalytic;
    data.historyHandles = [];

    setappdata(fig, 'interactiveData', data);

    set(fig, 'WindowButtonDownFcn', @interactivePhotonMouseDown);
    set(fig, 'WindowButtonMotionFcn', @interactivePhotonMouseMove);
    set(fig, 'WindowButtonUpFcn', @interactivePhotonMouseUp);
end


function interactivePhotonMouseDown(fig, ~)
    data = getappdata(fig, 'interactiveData');
    if isempty(data) || ~isstruct(data) || ~ishghandle(data.axMain)
        return;
    end

    xy = interactiveCurrentPoint(data.axMain);
    if isempty(xy)
        return;
    end

    xLim = get(data.axMain, 'XLim');
    yLim = get(data.axMain, 'YLim');
    if xy(1) < xLim(1) || xy(1) > xLim(2) || xy(2) < yLim(1) || xy(2) > yLim(2)
        return;
    end

    xy = snapPointToInteractiveGrid(data.p, xy, xLim, yLim);

    if strcmp(data.sourceMode, 'click')
        startPos = xy;
        data.sourcePoints = startPos;
        data.activeEmitterIdx = 1;

        if ishghandle(data.sourceHandle)
            set(data.sourceHandle, 'XData', startPos(1), 'YData', startPos(2));
        end

        if ishghandle(data.activeEmitterHandle)
            set(data.activeEmitterHandle, 'XData', startPos(1), 'YData', startPos(2));
        end

        data.isDragging = true;
        set(data.previewHandle, ...
            'XData', [startPos(1), xy(1)], ...
            'YData', [startPos(2), xy(2)], ...
            'Visible', 'on');
        refreshInteractiveInfoHandle(data, 'AIMING ...', [0.90, 0.60, 0.20], startPos, NaN, NaN, NaN);

        if ishghandle(data.bannerHandle)
            set(data.bannerHandle, ...
                'String', '', ...
                'Color', [0.98, 0.96, 0.88]);
        end

        setappdata(fig, 'interactiveData', data);
        return;
    end

    dAll = hypot(xy(1) - data.sourcePoints(:, 1), xy(2) - data.sourcePoints(:, 2));
    [distToPhoton, idxNearest] = min(dAll);
    data.activeEmitterIdx = idxNearest;
    startPos = data.sourcePoints(idxNearest, :);

    if ishghandle(data.activeEmitterHandle)
        set(data.activeEmitterHandle, 'XData', startPos(1), 'YData', startPos(2));
    end

    if distToPhoton <= data.p.interactivePickupRadius
        data.isDragging = true;
        set(data.previewHandle, ...
            'XData', [startPos(1), xy(1)], ...
            'YData', [startPos(2), xy(2)], ...
            'Visible', 'on');
        refreshInteractiveInfoHandle(data, 'AIMING ...', [0.90, 0.60, 0.20], startPos, NaN, NaN, NaN);

        if ishghandle(data.bannerHandle)
            set(data.bannerHandle, ...
                'String', 'Arahkan drag, lalu lepas untuk launch', ...
                'Color', [0.98, 0.96, 0.88]);
        end

        setappdata(fig, 'interactiveData', data);
        return;
    end

    data = interactiveLaunchPhoton(fig, data, startPos, xy);
    setappdata(fig, 'interactiveData', data);
end


function interactivePhotonMouseMove(fig, ~)
    data = getappdata(fig, 'interactiveData');
    if isempty(data) || ~isstruct(data) || ~ishghandle(data.axMain)
        return;
    end

    if ~data.isDragging
        if ishghandle(data.previewHandle)
            set(data.previewHandle, 'Visible', 'off');
        end
        return;
    end

    xy = interactiveCurrentPoint(data.axMain);
    if isempty(xy)
        return;
    end

    xLim = get(data.axMain, 'XLim');
    yLim = get(data.axMain, 'YLim');
    xy(1) = min(max(xy(1), xLim(1)), xLim(2));
    xy(2) = min(max(xy(2), yLim(1)), yLim(2));
    xy = snapPointToInteractiveGrid(data.p, xy, xLim, yLim);

    startPos = data.sourcePoints(data.activeEmitterIdx, :);

    set(data.previewHandle, ...
        'XData', [startPos(1), xy(1)], ...
        'YData', [startPos(2), xy(2)], ...
        'Visible', 'on');

    dragVec = xy - startPos;
    ang = atan2d(dragVec(2), dragVec(1));
    refreshInteractiveInfoHandle(data, 'AIMING ...', [0.90, 0.60, 0.20], startPos, ang, NaN, NaN);

    setappdata(fig, 'interactiveData', data);
end


function interactivePhotonMouseUp(fig, ~)
    data = getappdata(fig, 'interactiveData');
    if isempty(data) || ~isstruct(data) || ~ishghandle(data.axMain)
        return;
    end

    if ~data.isDragging
        return;
    end

    data.isDragging = false;
    xy = interactiveCurrentPoint(data.axMain);
    if isempty(xy)
        set(data.previewHandle, 'Visible', 'off');
        setappdata(fig, 'interactiveData', data);
        return;
    end

    xLim = get(data.axMain, 'XLim');
    yLim = get(data.axMain, 'YLim');
    xy(1) = min(max(xy(1), xLim(1)), xLim(2));
    xy(2) = min(max(xy(2), yLim(1)), yLim(2));
    xy = snapPointToInteractiveGrid(data.p, xy, xLim, yLim);

    startPos = data.sourcePoints(data.activeEmitterIdx, :);
    data = interactiveLaunchPhoton(fig, data, startPos, xy);
    setappdata(fig, 'interactiveData', data);
end


function data = interactiveLaunchPhoton(fig, data, startPos, targetPoint)
    data.isDragging = false;
    [traj, ok, msg] = integrateDirectedTrajectory(data.p, startPos, targetPoint);

    if ~ok
        set(data.previewHandle, 'Visible', 'off');
        dragVec = targetPoint - startPos;
        if all(isfinite(dragVec)) && norm(dragVec) > eps
            launchAngleDeg = atan2d(dragVec(2), dragVec(1));
        else
            launchAngleDeg = NaN;
        end
        refreshInteractiveInfoHandle(data, msg, [0.78, 0.22, 0.22], startPos, launchAngleDeg, NaN, NaN);
        if ishghandle(data.bannerHandle)
            set(data.bannerHandle, 'String', data.bannerBaseText, 'Color', [0.96, 0.96, 0.98]);
        end
        return;
    end

    oldX = get(data.rayHandle, 'XData');
    oldY = get(data.rayHandle, 'YData');
    if ~isempty(oldX) && any(~isnan(oldX))
        oldMainColor = get(data.rayHandle, 'Color');
        oldGlowColor = get(data.rayGlowHandle, 'Color');

        hOldGlow = plot(data.axMain, oldX, oldY, '-', ...
            'Color', oldGlowColor, ...
            'LineWidth', max(1.2, 0.75 * data.p.interactiveGlowWidth));

        hOldMain = plot(data.axMain, oldX, oldY, '-', ...
            'Color', oldMainColor, ...
            'LineWidth', max(1.0, 0.82 * data.p.interactiveLineWidth));

        oldHeadX = get(data.headHandle, 'XData');
        oldHeadY = get(data.headHandle, 'YData');
        oldHeadValid = ~isempty(oldHeadX) && ~isempty(oldHeadY) && any(~isnan(oldHeadX)) && any(~isnan(oldHeadY));

        if oldHeadValid
            hOldHead = plot(data.axMain, oldHeadX(end), oldHeadY(end), 'o', ...
                'Color', oldMainColor, ...
                'MarkerFaceColor', oldMainColor, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerSize', max(4.5, 0.78 * data.p.interactiveHeadSize));
            data.historyHandles = [data.historyHandles; hOldGlow; hOldMain; hOldHead];
        else
            data.historyHandles = [data.historyHandles; hOldGlow; hOldMain];
        end

        if data.p.interactiveHistoryLimit > 0
            maxHandles = 3 * data.p.interactiveHistoryLimit;
            while numel(data.historyHandles) > maxHandles
                if ishghandle(data.historyHandles(1))
                    delete(data.historyHandles(1));
                end
                data.historyHandles(1) = [];
            end
        end
    end

    [statusLabel, statusColor] = interactiveStatusVisual(traj.status);
    if data.p.interactiveUsePathSmoothing
        [xDraw, yDraw] = resampleTrajectoryPath(traj.x, traj.y, data.p.interactiveSmoothSamples);
    else
        xDraw = traj.x;
        yDraw = traj.y;
    end

    [xAll, yAll] = clipPathToLimits(xDraw, yDraw, data.p.interactiveXLim, data.p.interactiveYLim);
    glowColor = blendColor(statusColor, [0.02, 0.03, 0.05], data.p.interactiveGlowMix);

    set(data.previewHandle, 'Visible', 'off');
    set(data.rayGlowHandle, ...
        'XData', NaN, ...
        'YData', NaN, ...
        'Color', glowColor, ...
        'LineWidth', data.p.interactiveGlowWidth, ...
        'Visible', 'on');

    set(data.rayHandle, ...
        'XData', NaN, ...
        'YData', NaN, ...
        'Color', statusColor, ...
        'LineWidth', data.p.interactiveLineWidth, ...
        'Visible', 'on');

    set(data.headHandle, ...
        'XData', NaN, ...
        'YData', NaN, ...
        'Color', statusColor, ...
        'MarkerFaceColor', statusColor, ...
        'Visible', 'off');

    if data.p.interactiveAnimateLaunch
        nPts = numel(xAll);
        if nPts < 2
            idx = 1;
        else
            if data.p.interactiveEqualizePlaybackSpeed
                frameBudget = round(max(0.4, data.p.interactiveLaunchDurationSec) * max(data.p.interactiveFps, eps));
                frameBudget = max(24, min(data.p.interactiveMaxFrames, frameBudget));
                idx = buildFrameIndicesByPathLength(xAll, yAll, frameBudget, data.p.interactiveUseEasing);
            else
                frameBudget = max(24, min(data.p.interactiveMaxFrames, ceil(0.65 * nPts)));
                idx = buildFrameIndices(nPts, frameBudget, data.p.interactiveUseEasing);
            end
        end

        tRef = tic;

        for k = 1:numel(idx)
            if ~ishghandle(fig) || ~ishghandle(data.axMain)
                return;
            end

            j = idx(k);
            xSeg = xAll(1:j);
            ySeg = yAll(1:j);
            set(data.rayGlowHandle, 'XData', xSeg, 'YData', ySeg);
            set(data.rayHandle, 'XData', xSeg, 'YData', ySeg);

            lastValid = find(~isnan(xSeg), 1, 'last');
            if isempty(lastValid)
                set(data.headHandle, 'Visible', 'off');
            else
                set(data.headHandle, ...
                    'XData', xSeg(lastValid), ...
                    'YData', ySeg(lastValid), ...
                    'Visible', 'on');
            end

            drawnow;
            if data.p.interactiveFramePause
                targetTime = (k - 1) / max(data.p.interactiveFps, eps);
                waitTime = targetTime - toc(tRef);
                if waitTime > 0
                    pause(waitTime);
                end
            end
        end
    else
        set(data.rayGlowHandle, 'XData', xAll, 'YData', yAll);
        set(data.rayHandle, 'XData', xAll, 'YData', yAll);
        lastValid = find(~isnan(xAll), 1, 'last');
        if ~isempty(lastValid)
            set(data.headHandle, 'XData', xAll(lastValid), 'YData', yAll(lastValid), 'Visible', 'on');
        end
    end

    refreshInteractiveInfoHandle(data, statusLabel, statusColor, startPos, traj.launchAngleDeg, traj.impactEquivalent, traj.r(end));

    if ishghandle(data.bannerHandle)
        set(data.bannerHandle, 'String', data.bannerBaseText, 'Color', [0.96, 0.96, 0.98]);
    end
end


function [traj, ok, msg] = integrateDirectedTrajectory(p, startPos, targetPoint)
    traj = struct();
    ok = false;
    msg = 'DRAG TOO SHORT';

    dragVec = targetPoint - startPos;
    dragNorm = norm(dragVec);
    usedAutoLeft = false;
    if dragNorm < p.interactiveMinDrag
        if isfield(p, 'interactiveShortDragAutoLeft') && p.interactiveShortDragAutoLeft
            u = [-1, 0];
            usedAutoLeft = true;
            msg = 'AUTO LEFT (SHORT DRAG)';
        else
            return;
        end
    else
        u = dragVec / dragNorm;
        msg = 'OK';
    end

    x0 = startPos(1);
    y0cart = startPos(2);
    r0 = hypot(x0, y0cart);

    if r0 <= p.r_horizon * (1 + 1e-3)
        msg = 'START TOO CLOSE TO HORIZON';
        return;
    end

    phi0 = atan2(y0cart, x0);

    vrUnit = (x0 * u(1) + y0cart * u(2)) / r0;
    LUnit = x0 * u(2) - y0cart * u(1);

    metricFactor = 1 - 2 * p.M / r0;
    baseScale = vrUnit^2 + metricFactor * (LUnit^2 / r0^2);
    energyAvail = p.E^2 - metricFactor * p.kappa;

    if energyAvail <= 0 || baseScale <= eps
        msg = 'UNPHYSICAL DIRECTION';
        return;
    end

    scale = sqrt(energyAvail / baseScale);
    vr0 = scale * vrUnit;
    L = scale * LUnit;

    y0 = [r0; vr0; phi0];
    maxStepUse = p.maxStep;
    if isfield(p, 'interactiveSolverMaxStep')
        maxStepUse = min(p.maxStep, p.interactiveSolverMaxStep);
    end

    opts = odeset(...
        'RelTol', p.relTol, ...
        'AbsTol', p.absTol, ...
        'MaxStep', maxStepUse, ...
        'Events', @(lambda, y) geodesicEvents(lambda, y, p));

    try
        [lambda, y, ~, ~, ie] = ode45(@(lambda, y) geodesicODE(lambda, y, p, L), ...
            [0, p.lambdaMax], y0, opts);
    catch me
        msg = sprintf('SOLVER ERROR: %s', me.message);
        return;
    end

    r = y(:, 1);
    vr = y(:, 2);
    phi = y(:, 3);

    traj.lambda = lambda;
    traj.r = r;
    traj.vr = vr;
    traj.phi = phi;
    traj.x = r .* cos(phi);
    traj.y = r .* sin(phi);
    traj.L = L;
    traj.status = classifyByEvents(p, r, vr, ie);
    traj.impactEquivalent = L / max(p.E, eps);
    traj.launchAngleDeg = atan2d(u(2), u(1));
    traj.usedAutoLeft = usedAutoLeft;

    ok = true;
end


function sourcePoints = buildInteractiveSourcePoints(p)
    modeStr = getInteractiveSourceMode(p);

    if strcmp(modeStr, 'line')
        nSource = max(2, round(p.interactiveSourceCount));
        yMin = -abs(p.interactiveSourceHalfSpan);
        yMax = abs(p.interactiveSourceHalfSpan);
        yVals = linspace(yMin, yMax, nSource).';
        xVals = p.interactiveSourceX * ones(nSource, 1);
        sourcePoints = [xVals, yVals];
    else
        sourcePoints = p.interactiveStartPos(:).';
    end
end


function modeStr = getInteractiveSourceMode(p)
    modeStr = 'point';
    if isfield(p, 'interactiveSourceMode')
        modeStr = lower(strtrim(p.interactiveSourceMode));
    end

    if ~strcmp(modeStr, 'line') && ~strcmp(modeStr, 'point') && ~strcmp(modeStr, 'click')
        modeStr = 'point';
    end
end


function [xOut, yOut] = resampleTrajectoryPath(xIn, yIn, targetSamples)
    xOut = xIn;
    yOut = yIn;

    nPts = numel(xIn);
    if nPts < 3
        return;
    end

    targetSamples = max(40, round(targetSamples));
    if targetSamples <= nPts
        return;
    end

    ds = hypot(diff(xIn), diff(yIn));
    s = [0; cumsum(ds)];

    keep = [true; ds > 1e-12];
    s = s(keep);
    xKeep = xIn(keep);
    yKeep = yIn(keep);

    if numel(s) < 3 || s(end) <= 0
        return;
    end

    sQ = linspace(0, s(end), targetSamples).';

    xOut = interp1(s, xKeep, sQ, 'pchip');
    yOut = interp1(s, yKeep, sQ, 'pchip');
end


function idx = buildFrameIndices(nPts, nFrames, useEasing)
    if nPts < 1
        idx = 1;
        return;
    end

    if nFrames <= 1
        idx = nPts;
        return;
    end

    t = linspace(0, 1, nFrames).';
    if useEasing
        tUse = 0.5 - 0.5 * cos(pi * t);
    else
        tUse = t;
    end

    idx = 1 + (nPts - 1) * tUse;
    idx = unique(max(1, min(nPts, round(idx))));

    if idx(1) ~= 1
        idx = [1; idx];
    end
    if idx(end) ~= nPts
        idx = [idx; nPts];
    end
end


function idx = buildFrameIndicesByPathLength(xPath, yPath, nFrames, useEasing)
    nPts = numel(xPath);
    if nPts < 2
        idx = 1;
        return;
    end

    if nFrames <= 1
        idx = nPts;
        return;
    end

    xWork = xPath(:);
    yWork = yPath(:);
    valid = isfinite(xWork) & isfinite(yWork);

    if ~any(valid)
        idx = buildFrameIndices(nPts, nFrames, useEasing);
        return;
    end

    firstValid = find(valid, 1, 'first');
    for i = 1:(firstValid - 1)
        xWork(i) = xWork(firstValid);
        yWork(i) = yWork(firstValid);
        valid(i) = true;
    end

    for i = 2:nPts
        if ~valid(i)
            xWork(i) = xWork(i - 1);
            yWork(i) = yWork(i - 1);
            valid(i) = true;
        end
    end

    ds = hypot(diff(xWork), diff(yWork));
    s = [0; cumsum(ds)];
    [sUnique, ia] = unique(s, 'stable');

    if numel(sUnique) < 2 || sUnique(end) <= 0
        idx = buildFrameIndices(nPts, nFrames, useEasing);
        return;
    end

    idxBase = (1:nPts).';
    idxUnique = idxBase(ia);

    t = linspace(0, 1, nFrames).';
    if useEasing
        t = 0.5 - 0.5 * cos(pi * t);
    end

    sQuery = t * sUnique(end);
    idxInterp = interp1(sUnique, idxUnique, sQuery, 'linear');
    idx = unique(max(1, min(nPts, round(idxInterp))));

    if idx(1) ~= 1
        idx = [1; idx];
    end
    if idx(end) ~= nPts
        idx = [idx; nPts];
    end
end


function c = blendColor(cA, cB, mixVal)
    mixVal = max(0, min(1, mixVal));
    cA = max(0, min(1, cA));
    cB = max(0, min(1, cB));
    c = (1 - mixVal) * cA + mixVal * cB;
end


function xyOut = snapPointToInteractiveGrid(p, xyIn, xLim, yLim)
    xyOut = xyIn;
    if isfield(p, 'interactiveSnapToGrid') && p.interactiveSnapToGrid
        h = max(eps, p.interactiveGridStep);
        xyOut = h * round(xyOut / h);
    end

    xyOut(1) = min(max(xyOut(1), xLim(1)), xLim(2));
    xyOut(2) = min(max(xyOut(2), yLim(1)), yLim(2));
end


function refreshInteractiveInfoHandle(data, statusText, statusColor, sourcePos, launchAngleDeg, impactEquivalent, finalRadius)
    if ~isstruct(data) || ~isfield(data, 'infoHandle') || ~ishghandle(data.infoHandle)
        return;
    end

    infoText = composeInteractiveInfoText( ...
        data.p, data.bCritNum, data.bAnalytic, statusText, sourcePos, launchAngleDeg, impactEquivalent, finalRadius);
    set(data.infoHandle, 'String', infoText, 'Color', statusColor);
end


function infoText = composeInteractiveInfoText(p, bCritNum, bAnalytic, statusText, sourcePos, launchAngleDeg, impactEquivalent, finalRadius)
    if nargin < 4 || isempty(statusText)
        statusText = 'READY';
    end

    if isnan(bAnalytic)
        analyticText = 'N/A';
    else
        analyticText = sprintf('%.6f', bAnalytic);
    end

    if nargin < 5 || numel(sourcePos) < 2 || any(~isfinite(sourcePos))
        sourceXText = '--';
        sourceYText = '--';
    else
        sourceXText = sprintf('%+8.3f M', sourcePos(1) / p.M);
        sourceYText = sprintf('%+8.3f M', sourcePos(2) / p.M);
    end

    if nargin < 6 || ~isfinite(launchAngleDeg)
        launchText = '--';
    else
        launchText = sprintf('%+5.0f deg', round(launchAngleDeg));
    end

    if nargin < 7 || ~isfinite(impactEquivalent)
        impactText = '--';
    else
        impactText = sprintf('%+8.4f', impactEquivalent);
    end

    if nargin < 8 || ~isfinite(finalRadius)
        finalRadiusText = '--';
    else
        finalRadiusText = sprintf('%8.4f M', finalRadius / p.M);
    end

    infoText = sprintf([ ...
        'b_{crit,num} = %.6f\n' ...
        'b_{crit,ana} = %s\n' ...
        'r_h          = %.3f M\n' ...
        'STATUS       : %s\n' ...
        '\n' ...
        'Source x     : %s\n' ...
        'Source y     : %s\n' ...
        'Launch angle : %s\n' ...
        'Impact b~    : %s\n' ...
        'Final radius : %s'], ...
        bCritNum, analyticText, p.r_horizon / p.M, statusText, ...
        sourceXText, sourceYText, launchText, impactText, finalRadiusText);
end


function [label, color] = interactiveStatusVisual(status)
    if strcmp(status, 'escaped')
        label = 'ESCAPED (SAFE)';
        color = [0.24, 0.86, 0.98];
    else
        label = 'CAPTURED BY BLACK HOLE';
        color = [0.96, 0.45, 0.26];
    end
end


function xy = interactiveCurrentPoint(ax)
    xy = [];
    if isempty(ax) || ~ishghandle(ax)
        return;
    end

    cp = get(ax, 'CurrentPoint');
    if isempty(cp) || size(cp, 2) < 2
        return;
    end

    xy = cp(1, 1:2);
end


function [xOut, yOut] = clipPathToLimits(xIn, yIn, xLim, yLim)
    xOut = xIn;
    yOut = yIn;

    mask = (xIn >= xLim(1)) & (xIn <= xLim(2)) & ...
           (yIn >= yLim(1)) & (yIn <= yLim(2));

    xOut(~mask) = NaN;
    yOut(~mask) = NaN;
end