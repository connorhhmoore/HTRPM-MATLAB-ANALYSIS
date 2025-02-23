%% Analysis of Sphere Packing in the HTR-PM Reactor
% Connor Moore, <connor.moore@ontariotechu.net>, 2024-2025
close all; clear;
% Reactor core parameters
Dc   = 3;           % Diameter of Core
Hc   = 11;          % Height of core [m]
Vc = pi/4*Dc^2*Hc;  % Volume of core [m³]
npts = 50;          % Cylinder discretization
capacity = 250e6;   % Thermal capacity [W]

k_graphite = 50.8; % [W/(m⋅K)­­­]
% https://mooseframework.inl.gov/bison/source/materials/GraphiteMatrixThermal.html

% Spherical fuel pellet parameters
Df = 4e-1;      % Diameter of fuel pellet [m]
As = pi*Df^2;   % Surface area of fuel pellet [m²]
Vf = 4/3*pi*(Df/2)^3;

% Flow area calculation
n_flow_points = 300;

% Real performance of reactor
N_real = 425e3;
As_real = N_real*As; % [m²]
q_real = capacity/As_real; % [W/m²]
e_gen_real = capacity/N_real; % [W/n]
T_max_real = e_gen_real*(Df/2)^2/(6*k_graphite);

% Max side length
Nmax = floorDiv(Dc,Df);   % Maximum number of side:to:side pellets

[Xf,Yf,Zf] = meshgrid(-Nmax/2:Nmax/2,-Nmax/2:Nmax/2,0:Hc/Df*30);

Xf = Xf.*Df;
Yf = Yf.*Df;
Zf = Zf.*Df + Df/2;

%% For basic packing ----------------------------------------------------

% Calculate internal points
inside_rad = sqrt(Xf(:).^2 + Yf(:).^2) <= Dc/2-Df/2;
inside_height = Zf(:) < Hc-(Df/2);
key = and(inside_rad,inside_height);

Xin_basic = Xf(key);
Yin_basic = Yf(key);
Zin_basic = Zf(key);

% Tally points inside reactor region
N_basic = nnz(key);
As_basic = N_basic*As;

%% For HCP packing ------------------------------------------------------

Xf_hcp = Xf;
Yf_hcp = Yf;
Zf_hcp = Zf;

% Applying transforms to correct:
Xf_hcp(:,:,2:2:end) = Xf(:,:,2:2:end) + Df/2;
Yf_hcp(:,:,2:2:end) = Yf(:,:,2:2:end) + Df/2;
%Zf_hcp(:,:,2:2:end) = Zf(:,:,2:2:end); %+ sqrt(Df^2-Df^2/2);

% Displace all sheets to hug each other
for i = 2:numel(unique(Zf_hcp))
    Zf_hcp(:,:,i:end) = Zf_hcp(:,:,i:end) - (Df - sqrt(6)/3*Df);
    fprintf("%i/%i HCP displacement calculation\n",i,numel(unique(Zf_hcp)));
end

% Calculate internal points
inside_rad = sqrt(Xf_hcp(:).^2 + Yf_hcp(:).^2) <= Dc/2-Df/2;
inside_height = Zf_hcp(:) < Hc-(Df/2);
key = and(inside_rad,inside_height);

Xin_hcp = Xf_hcp(key);
Yin_hcp = Yf_hcp(key);
Zin_hcp = Zf_hcp(key);

N_hcp = nnz(key);
As_hcp = N_hcp*As;

% For non-tiny pebble diameters, plot the geometry
if Df>=0.4
    %% Generate 3D plot of the packing problem ------------------------------
    [X,Y,Z] = cylinder(Dc/2,npts);
    Z = Z.*Hc; % Scale to correct height

    % Plot the geometry
    alpha = 0.12;

    % Plotting flat/box stacking --------------------------------------------
    figure(Name="HCP Stacking");
    surf(X,Y,Z,'FaceAlpha',alpha,'EdgeColor','none');
    hold on;
    patch(X(1,:),Y(1,:),Z(1,:),'b','facealpha',alpha,'edgecolor','none');
    patch(X(2,:),Y(2,:),Z(2,:),'b','facealpha',alpha,'edgecolor','none');

    for h = 1:numel(Xin_hcp)
        [i,j,k] = sphere(12);
        i = Df/2.*i + Xin_hcp(h);
        j = Df/2.*j + Yin_hcp(h);
        k = Df/2.*k + Zin_hcp(h);
        surface(i,j,k);
        fprintf("%i/%i sphere display calculation\n",h,numel(Xin_basic))
    end

    hold off;
    grid on;
    axis equal;
    xlabel('{\itx} [m]');
    ylabel('{\ity} [m]');
    zlabel('{\itz} [m]');

    % Plotting HCP stacking ------------------------------------------------
    figure(Name="Basic Stacking")
    surf(X,Y,Z,'FaceAlpha',alpha,'EdgeColor','none');
    hold on;
    patch(X(1,:),Y(1,:),Z(1,:),'b','facealpha',alpha,'edgecolor','none');
    patch(X(2,:),Y(2,:),Z(2,:),'b','facealpha',alpha,'edgecolor','none');

    for h = 1:numel(Xin_basic)
        [i,j,k] = sphere(12);
        i = Df/2.*i + Xin_basic(h);
        j = Df/2.*j + Yin_basic(h);
        k = Df/2.*k + Zin_basic(h);
        surface(i,j,k);
        fprintf("%i/%i sphere display calculation\n",h,numel(Xin_hcp))
    end

    hold off;
    grid on;
    axis equal;
    xlabel('{\itx} [m]');
    ylabel('{\ity} [m]');
    zlabel('{\itz} [m]');
end

%% Flow area calculations
% Range [0,√(6)*(2/3)*Df] [m] axially inside core

L = linspace(0,sqrt(6)*2/3*Df,n_flow_points);
n_per_layer_hcp = sum(Zin_hcp == Zin_hcp(1));
XS_area_total = pi/4*Dc^2;

bounds_hcp = Df.*[0;        % 1. MIDDLE of 1st sphere
         (sqrt(6)/3-0.5);   % 2. START of 2nd sphere
         (0.5);             % 3. END of 1st sphere
         (sqrt(6)/3);       % 4. MIDDLE of 2nd sphere
         (4*sqrt(6)-3)/6;   % 5. START of 3rd sphere
         (sqrt(6)/3+0.5);   % 6. END of 2nd sphere
         2*sqrt(6)/3];      % 7. MIDDLE of 3rd sphere

XS_area_hcp = @(L) XS_area_total - n_per_layer_hcp .* ((L<=bounds_hcp(3)).*pi.*((Df/2)^2-L.^2) ...
                    + and(L>=bounds_hcp(2), L<=bounds_hcp(6)) .* pi .* ((Df/2)^2-(L-bounds_hcp(4)).^2) ...
                    + (L>=bounds_hcp(5)) .* pi .* ((Df/2)^2-(L-bounds_hcp(7)).^2));

XS_flow_area_hcp = XS_area_hcp(L);

XS_area_avg_hcp = integral(XS_area_hcp,0,bounds_hcp(7))/bounds_hcp(7);

figure(Name="XS Flow Area Plot")
plot(L,XS_flow_area_hcp,'-b')
axis padded
grid on
xlabel("Axial Distance [m]")
ylabel("Cross-sectional Flow Area [m²]")
xticks([0,sqrt(6)/3,2*sqrt(6)/3]*Df)
xticklabels({'0.0 D_f','√(6)/3 D_f','2·√(6)/3 D_f'})
hold on
yline(XS_area_avg_hcp,'--r')
text(bounds_hcp(4),XS_area_avg_hcp+0.05,"Avg. "+XS_area_avg_hcp+" [m²]",...
    "VerticalAlignment","bottom","HorizontalAlignment","center","Color",'r')
hold off

%% Output and postprocessing --------------------------------------------

q_basic = capacity/As_basic;
dens_basic = N_basic/Vc;

q_hcp = capacity/As_hcp;
dens_hcp = N_hcp/Vc;

fprintf(" Calculated values using a core D/H (%3.2f [m])/(%3.2f [m]).\n Fuel diameter %3.2f [m].\n",Dc,Hc,Df);

fprintf("\n ======== Basic stacking (box/cube stacking) pattern ========\n")
fprintf("  -> Total number of contained pellets: %i [n]\n",N_basic);
fprintf("  -> Total surface area of all pellets: %3.2f [m²]\n",As_basic);
fprintf("  -> Calculated surface heat flux of the core: %3.2f [W/m²]\n",q_basic);
fprintf("  -> Average number of pellets per cubic meter: %3.2f [n/m³]\n",dens_basic);
fprintf("  -> Volume packing efficiency: %3.2f [%%]\n",N_basic*Vf/Vc*100);

fprintf("\n ======== HCP stacking (sheet-based stacking) pattern ========\n")
fprintf("  -> Total number of contained pellets: %i [n]\n",N_hcp);
fprintf("  -> Total surface area of all pellets: %3.2f [m²]\n",As_hcp);
fprintf("  -> Calculated surface heat flux of the core: %3.2f [W/m²]\n",q_hcp);
fprintf("  -> Average number of pellets per cubic meter: %3.2f [n/m³]\n",dens_hcp);
fprintf("  -> Volume packing efficiency: %3.2f [%%]\n",N_hcp*Vf/Vc*100);

fprintf("\n ======== Real analysis (no stacking calculation) ========\n")
fprintf("  -> Total claimed number of pellets: %i [n]\n",N_real)
fprintf("  -> Volume packing efficiency: %3.2f [%%]\n",N_real*Vf/Vc*100);
fprintf("  -> Total surface area of all pellets: %3.2f [m²]\n",As_real);
fprintf("  -> Calculated surface heat flux of the core: %6.2f [W/m²]\n",q_real)
fprintf("  -> ΔT to center of fuel pellet: %4.3f [°C]\n",T_max_real)

fprintf("\n ======== Flow cross-section calculations for HCP ========\n")
fprintf("  -> Center-to-center sphere spacing: %5.3f [m]\n",uniquetol(diff(unique(Zin_hcp))))
fprintf("  -> XS flow area for entire reactor core %5.3f [m²]\n",pi/4*Dc^2)
fprintf("  -> Maximum XS area for one pebble: %5.3f [m²]\n",pi/4*Df^2)
fprintf("  -> Maximum XS area for reactor core: %5.3f [m²]\n",max(XS_flow_area_hcp))
fprintf("  -> Minimum XS area for reactor core: %5.3f [m²]\n",min(XS_flow_area_hcp))
fprintf("  -> Average XS flow area for reactor core %5.3f [m²]\n",XS_area_avg_hcp)

fprintf("\n (Sphere packing analysis, Connor Moore, 2024-2025)")
fprintf("\n      <connor.moore@ontariotechu.net>\n")
fprintf(" (Calculations performed "+string(datetime)+")\n")