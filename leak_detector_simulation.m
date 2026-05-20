%% Leak detector simulation in a closed rectangular room
% Topic: concentration fields from a leak in an enclosed space.
% The model solves a 2D advection-diffusion equation with a point leak source
% and emulates a mobile leak detector that scans the room and records
% concentration versus time.

clear; clc; close all;

%% Geometry and grid
Lx = 10;                  % room length, m
Ly = 6;                   % room width, m
Nx = 140;                 % grid points along x
Ny = 84;                  % grid points along y
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
[x, y] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));

%% Physical parameters
D  = 1.4e-2;              % effective diffusion coefficient, m^2/s
uX = 0.02;                % weak recirculation along x, m/s
uY = 0.00;                % m/s

% Stable time step (explicit scheme)
dt_diff = 0.25*min(dx^2, dy^2)/D;
dt_adv  = 0.45/min(max(abs(uX)/dx, 1e-12), max(abs(uY)/dy, 1e-12));
dt = min(dt_diff, dt_adv);
Tend = 240;               % simulation time, s
Nt = ceil(Tend/dt);

%% Leak source
xLeak = 2.4; yLeak = 1.7; % leak location, m
Q = 8.0;                  % source strength, concentration units / s
sigma = 0.22;             % source spread, m
src = exp(-((x-xLeak).^2 + (y-yLeak).^2)/(2*sigma^2));
src = src / sum(src(:));  % normalize for mesh-independent injection

%% Detector trajectory (serpentine scan)
scanLines = 8;
yLines = linspace(0.5, Ly-0.5, scanLines);
pathX = [];
pathY = [];
for k = 1:scanLines
    if mod(k,2)==1
        pathX = [pathX, linspace(0.4, Lx-0.4, 220)]; %#ok<AGROW>
    else
        pathX = [pathX, linspace(Lx-0.4, 0.4, 220)]; %#ok<AGROW>
    end
    pathY = [pathY, yLines(k)*ones(1,220)]; %#ok<AGROW>
end
pathT = linspace(0, Tend, numel(pathX));

%% Simulation fields
C = zeros(Ny, Nx);               % concentration map
sensorSignal = zeros(1, Nt);     % detector output
sensorX = interp1(pathT, pathX, linspace(0, Tend, Nt), 'linear', 'extrap');
sensorY = interp1(pathT, pathY, linspace(0, Tend, Nt), 'linear', 'extrap');

%% Time integration
for n = 1:Nt
    % Spatial derivatives
    dCdx = zeros(size(C)); dCdy = zeros(size(C));
    d2Cdx2 = zeros(size(C)); d2Cdy2 = zeros(size(C));

    dCdx(:,2:end-1) = (C(:,3:end)-C(:,1:end-2))/(2*dx);
    dCdy(2:end-1,:) = (C(3:end,:)-C(1:end-2,:))/(2*dy);

    d2Cdx2(:,2:end-1) = (C(:,3:end)-2*C(:,2:end-1)+C(:,1:end-2))/dx^2;
    d2Cdy2(2:end-1,:) = (C(3:end,:)-2*C(2:end-1,:)+C(1:end-2,:))/dy^2;

    % PDE: dC/dt = D*Laplace(C) - u·grad(C) + Q*src
    C = C + dt*(D*(d2Cdx2+d2Cdy2) - uX*dCdx - uY*dCdy + Q*src);

    % Impermeable walls (zero normal gradient)
    C(:,1)   = C(:,2);
    C(:,end) = C(:,end-1);
    C(1,:)   = C(2,:);
    C(end,:) = C(end-1,:);

    % Detector reading (bilinear interpolation)
    sx = sensorX(n); sy = sensorY(n);
    i = min(max(floor(sx/dx)+1,1),Nx-1);
    j = min(max(floor(sy/dy)+1,1),Ny-1);
    tx = (sx - (i-1)*dx)/dx;
    ty = (sy - (j-1)*dy)/dy;

    c00 = C(j,i);   c10 = C(j,i+1);
    c01 = C(j+1,i); c11 = C(j+1,i+1);
    sensorSignal(n) = (1-tx)*(1-ty)*c00 + tx*(1-ty)*c10 + (1-tx)*ty*c01 + tx*ty*c11;
end

time = (0:Nt-1)*dt;

%% Visualisation
figure('Color','w','Position',[80 80 1200 520]);

subplot(1,2,1);
imagesc([0 Lx],[0 Ly],C);
set(gca,'YDir','normal');
axis equal tight;
hold on;
plot(pathX, pathY, 'w--', 'LineWidth', 1.0);
plot(sensorX(end), sensorY(end), 'wo', 'MarkerFaceColor','k', 'MarkerSize', 6);
plot(xLeak, yLeak, 'rp', 'MarkerFaceColor','r', 'MarkerSize', 14);
colorbar;
xlabel('x, m'); ylabel('y, m');
title('Concentration field at final time');
legend({'Detector trajectory','Detector current position','Leak point'}, 'Location','northoutside');

subplot(1,2,2);
plot(time, sensorSignal, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time, s'); ylabel('Measured concentration, a.u.');
title('Leak detector response during scan');

sgtitle('Simulation of leak localisation in enclosed space');
