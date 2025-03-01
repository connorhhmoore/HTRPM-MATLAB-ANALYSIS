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
Df = 6e-2;      % Diameter of fuel pellet [m]
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

%% For BCC packing ------------------------------------------------------

Xf_bcc = Xf;
Yf_bcc = Yf;
Zf_bcc = Zf;

% Applying transforms to correct:
Xf_bcc(:,:,2:2:end) = Xf(:,:,2:2:end) + Df/2;
Yf_bcc(:,:,2:2:end) = Yf(:,:,2:2:end) + Df/2;
%Zf_bcc(:,:,2:2:end) = Zf(:,:,2:2:end); %+ sqrt(Df^2-Df^2/2);

% Displace all sheets to hug each other
for i = 2:numel(unique(Zf_bcc))
    Zf_bcc(:,:,i:end) = Zf_bcc(:,:,i:end) - (Df - sqrt(6)/3*Df);
    fprintf("%i/%i BCC displacement calculation\n",i,numel(unique(Zf_bcc)));
end

% Calculate internal points
inside_rad = sqrt(Xf_bcc(:).^2 + Yf_bcc(:).^2) <= Dc/2-Df/2;
inside_height = Zf_bcc(:) < Hc-(Df/2);
key = and(inside_rad,inside_height);

Xin_bcc = Xf_bcc(key);
Yin_bcc = Yf_bcc(key);
Zin_bcc = Zf_bcc(key);

N_bcc = nnz(key);
As_bcc = N_bcc*As;

% For non-tiny pebble diameters, plot the geometry
if Df>=0.4
    %% Generate 3D plot of the packing problem ------------------------------
    [X,Y,Z] = cylinder(Dc/2,npts);
    Z = Z.*Hc; % Scale to correct height

    % Plot the geometry
    alpha = 0.03;

    % Plotting BCC stacking --------------------------------------------
    figure(Name="BCC Stacking");
    surf(X,Y,Z,'FaceAlpha',alpha,'EdgeColor','none','FaceColor','b');
    hold on;
    patch(X(1,:),Y(1,:),Z(1,:),'b','facealpha',alpha,'edgecolor','none','FaceColor','b');
    patch(X(2,:),Y(2,:),Z(2,:),'b','facealpha',alpha,'edgecolor','none','FaceColor','b');

    for h = 1:numel(Xin_bcc)
        [i,j,k] = sphere(16);
        i = Df/2.*i + Xin_bcc(h);
        j = Df/2.*j + Yin_bcc(h);
        k = Df/2.*k + Zin_bcc(h);
        surface(i,j,k);
        fprintf("%i/%i sphere display calculation\n",h,numel(Xin_bcc))
    end

    hold off;
    grid on;
    axis equal;
    colormap([linspace(1,.5)',linspace(0,.5)',linspace(0,.5)'])
    c = colorbar('Ticks',[0,5.5,10.82],'TickLabels',[800,550,300]);
    c.Label.String='Approximate Temperature [°C]';
    xlabel('{\itx} [m]');
    ylabel('{\ity} [m]');
    zlabel('{\itz} [m]');
%{
    % Plotting flat/box stacking ------------------------------------------------
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
        fprintf("%i/%i sphere display calculation\n",h,numel(Xin_basic))
    end

    hold off;
    grid on;
    axis equal;
    xlabel('{\itx} [m]');
    ylabel('{\ity} [m]');
    zlabel('{\itz} [m]');
%}
end

%% Flow area calculations
% Range [0,√(6)*(2/3)*Df] [m] axially inside core

L = linspace(0,sqrt(6)*2/3*Df,n_flow_points);
n_per_layer_bcc = sum(Zin_bcc == Zin_bcc(1));
XS_area_total = pi/4*Dc^2;

bounds_bcc = Df.*[0;        % 1. MIDDLE of 1st sphere
         (sqrt(6)/3-0.5);   % 2. START of 2nd sphere
         (0.5);             % 3. END of 1st sphere
         (sqrt(6)/3);       % 4. MIDDLE of 2nd sphere
         (4*sqrt(6)-3)/6;   % 5. START of 3rd sphere
         (sqrt(6)/3+0.5);   % 6. END of 2nd sphere
         2*sqrt(6)/3];      % 7. MIDDLE of 3rd sphere

XS_area_bcc = @(L) XS_area_total - n_per_layer_bcc .* ((L<=bounds_bcc(3)).*pi.*((Df/2)^2-L.^2) ...
                    + and(L>=bounds_bcc(2), L<=bounds_bcc(6)) .* pi .* ((Df/2)^2-(L-bounds_bcc(4)).^2) ...
                    + (L>=bounds_bcc(5)) .* pi .* ((Df/2)^2-(L-bounds_bcc(7)).^2));

XS_flow_area_bcc = XS_area_bcc(L);

XS_area_avg_bcc = integral(XS_area_bcc,0,bounds_bcc(7))/bounds_bcc(7);

figure(Name="XS Flow Area Plot")
plot(L,XS_flow_area_bcc,'-b')
axis padded
grid on
xlabel("Axial Distance [m]")
ylabel("Cross-sectional Flow Area [m²]")
xticks([0,sqrt(6)/3,2*sqrt(6)/3]*Df)
xticklabels({'0.0 D_f','√(6)/3 D_f','2·√(6)/3 D_f'})
hold on
yline(XS_area_avg_bcc,'--r')
text(bounds_bcc(4),XS_area_avg_bcc+0.05,"Avg. "+XS_area_avg_bcc+" [m²]",...
    "VerticalAlignment","bottom","HorizontalAlignment","center","Color",'r')
hold off

%% Output and postprocessing --------------------------------------------

q_basic = capacity/As_basic;
dens_basic = N_basic/Vc;

q_bcc = capacity/As_bcc;
dens_bcc = N_bcc/Vc;

fprintf(" Calculated values using a core D/H (%3.2f [m])/(%3.2f [m]).\n Fuel diameter %3.2f [m].\n",Dc,Hc,Df);

fprintf("\n ======== Basic stacking (box/cube stacking) pattern ========\n")
fprintf("  -> Total number of contained pellets: %i [n]\n",N_basic);
fprintf("  -> Total surface area of all pellets: %3.2f [m²]\n",As_basic);
fprintf("  -> Calculated surface heat flux of the core: %3.2f [W/m²]\n",q_basic);
fprintf("  -> Average number of pellets per cubic meter: %3.2f [n/m³]\n",dens_basic);
fprintf("  -> Average number of pellets per sheet: %3.2f [n/sheet]\n",N_basic./(numel(unique(Zin_basic))));
fprintf("  -> Volume packing efficiency: %3.2f [%%]\n",N_basic*Vf/Vc*100);

fprintf("\n ======== BCC-style stacking (sheet-based stacking) pattern ========\n")
fprintf("  -> Total number of contained pellets: %i [n]\n",N_bcc);
fprintf("  -> Total surface area of all pellets: %3.2f [m²]\n",As_bcc);
fprintf("  -> Calculated surface heat flux of the core: %3.2f [W/m²]\n",q_bcc);
fprintf("  -> Average number of pellets per cubic meter: %3.2f [n/m³]\n",dens_bcc);
fprintf("  -> Average number of pellets per sheet: %3.2f [n/sheet]\n",N_bcc./(numel(unique(Zin_bcc))));
fprintf("  -> Volume packing efficiency: %3.2f [%%]\n",N_bcc*Vf/Vc*100);

fprintf("\n ======== Real analysis (no stacking calculation) ========\n")
fprintf("  -> Total claimed number of pellets: %i [n]\n",N_real)
fprintf("  -> Volume packing efficiency: %3.2f [%%]\n",N_real*Vf/Vc*100);
fprintf("  -> Total surface area of all pellets: %3.2f [m²]\n",As_real);
fprintf("  -> Calculated surface heat flux of the core: %6.2f [W/m²]\n",q_real)
fprintf("  -> ΔT to center of fuel pellet: %4.3f [°C]\n",T_max_real)

fprintf("\n ======== Flow cross-section calculations for BCC ========\n")
fprintf("  -> Center-to-center sphere spacing: %5.3f [m]\n",uniquetol(diff(unique(Zin_bcc))))
fprintf("  -> XS flow area for entire reactor core %5.3f [m²]\n",pi/4*Dc^2)
fprintf("  -> Maximum XS area for one pebble: %5.3f [m²]\n",pi/4*Df^2)
fprintf("  -> Maximum XS area for reactor core: %5.3f [m²]\n",max(XS_flow_area_bcc))
fprintf("  -> Minimum XS area for reactor core: %5.3f [m²]\n",min(XS_flow_area_bcc))
fprintf("  -> Average XS flow area for reactor core %5.3f [m²]\n",XS_area_avg_bcc)

fprintf("\n (Sphere packing analysis, Connor Moore, 2024-2025)")
fprintf("\n      <connor.moore@ontariotechu.net>\n")
fprintf(" (Calculations performed "+string(datetime)+")\n")