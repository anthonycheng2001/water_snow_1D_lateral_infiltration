function [phi,h,theta_w,theta_i,x,tdelt,T,Da,Ht,Pe,St] = WWCS_func_rev6(warm_temp, cold_temp, run_time)

%% Numerical Solutions for Warm Water Flowing Through Cold Snow, Dept-Independent Porosity
% This code simulates the lateral infiltration of water through porous ice.
% It takes the arguments of source (warm_temp), sink (cold_temp), and run
% time (run_time) to return the results of porosity (phi), water depth (h),
% water temperature (theta_w), ice temperature (theta_i), space (x), time
% step (tdelt), timescale (T), Darcy number (Da), heat tranfer timescales
% ratio (Ht), Peclet nuber (Pe), and Stefan number (St).

%% Prefactors, Groups, Constants

% scales
wt = warm_temp; % temperature upper bound [Celsius]
ct = cold_temp; % temperature lower bound [Celsius]
THETA = wt-ct; % temperature scale [Celsius]

% H = 100; % riley depth scale [m]
% X = 10000; % riley length scale [m]

H = 1; % jess depth scale [m]
X = 260; % jess length scale [m]

% (approximate orders of) physical constants
PHI_s = 0.5; % sphericity [unitless ratio]
d_p = 0.0004; % diameter of equivalent volume spheres [m]
k_0 = PHI_s^2*d_p^2/180; % permeability (Kozeny) constant [dimensionless]
% k_0 = 1e-12; % [m^2], just a check to compare with above permeability
rho = 900; % density [kg/m^3]
g = 9.81; % gravity [m/s^2]
mu = 1e-3; % viscosity [Pa s]
sigma = 5; % heat transfer coefficient [W/m^3/K] (extra areal dimension)
K = 0.5; % thermal conductivity [W/m/K]
c = 2000; % specific heat capacity [J/kg/K]
L = 3e5; % specific latent heat of fusion [J/kg]

% stefan number
St = L/c/THETA;

% timescale groups uncovered during nondimensionalization
T = mu*X^2/k_0/rho/g/H; % timescale [s]
Ht = sigma*mu*X^2/k_0/rho^2/g/c/H; % heat transfer timescales ratio
Da = rho*g*H*T*k_0/mu/X^2; % darcy number, ratio of porous velocity to nominal velocity
Pe = K*mu/k_0/rho^2/g/c/H; % peclet number, ratio of advection to diffusion (conduction)

% baseline porosity without water
phi_0 = 0.5;

% tolerance
tol = 1e-5;

% water minimum temperature
water_min_temp = 0;

%% Initialization

% step sizes
xx = linspace(0,1,101);
xdelt = xx(2)-xx(1);
tdelt = xdelt^2/Da/3;

% mesh
rt = run_time;
t = 0:tdelt:tdelt*rt;
x = 0:xdelt:1;

% create state variables
phi = phi_0*ones(length(t), length(x)); % all same initial value
h = tol*ones(length(t), length(x)); % small seeds to grow from
theta_w = water_min_temp*ones(length(t),length(x)); % no initial water temp
theta_w(1,1) = wt; % source water temp
theta_i = ct*ones(length(t), length(x)); % ice begins at sink temp
htw = water_min_temp*ones(size(theta_w)); % advected warm quantity
htw(1,1) = wt; % source

%% Initial, Boundary Conditions
% some are used, some commented out for different cases if desired

% baseline porosity everywhere (IC)
phi(1,:) = phi_0;

% water depth is full at the edge, none at other (IC+BC)
% h(:,1) = 1;
% h(:,end) = 0;

% flux boundary for h on cold side (better BC)
h(:,1) = 1;
h(:,end) = h(:,end-1);

% for value boundary conditions (less realistic)
% theta_w(:,1) = wt;
% theta_w(:,2:end) = 0;
% theta_i(:,1) = 0;
% theta_i(:,2:end) = ct;

% for flux boundary conditions (water on cold side, ice both sides)
theta_w(:,1) = wt;
theta_w(:,end) = theta_w(:,end-1);
theta_i(:,1) = theta_i(:,2);
theta_i(:,end) = theta_i(:,end-1);

% advected term, same BC as water temp
htw(:,1) = wt;
htw(:,end) = htw(:,end-1);

%% Fill Mesh

for j = 2:length(t)
    % porosity (no boundary conditions)
    phi(j,:) = phi(j-1,:)+Ht.*h(j-1,:).*tdelt.*((1-phi(j-1,:)).*theta_i(j-1,:)+...
    phi(j-1,:).*theta_w(j-1,:))./(-theta_i(j-1,:)+theta_w(j-1,:)+St);

    % crop
    phi(phi > 1-tol) = 1-tol;
    phi(phi < tol) = tol;

    % water depth
    h(j,2:end-1) = Da.*tdelt./2./phi(j-1,2:end-1)./xdelt^2.*(h(j-1,3:end).* ...
    (phi(j-1,3:end)).^3.*(h(j-1,3:end)-h(j-1,2:end-1))-h(j-1,1:end-2).* ...
    (phi(j-1,1:end-2)).^3.*(h(j-1,2:end-1)-h(j-1,1:end-2)))+h(j-1,2:end-1);

    % doctor boundary conditions
    h(j,1) = 1;
    h(j,end) = h(j,end-1);

    % crop
    h(h > 1) = 1;
    h(h < tol) = tol;

    % water temperature (if solving explicitly)
    % theta_w(j,2:end-1) = h(j-1,2:end-1).*theta_w(j-1,2:end-1)./h(j,2:end-1)+...
    %     Da.*tdelt./2./h(j,2:end-1)./xdelt^2.*((phi(j-1,3:end)).^3.*h(j-1,3:end)...
    %     .*theta_w(j-1,3:end).*(h(j-1,3:end)-h(j-1,2:end-1))-(phi(j-1,1:end-2)).^3.*...
    %     h(j-1,1:end-2).*theta_w(j-1,1:end-2).*(h(j-1,2:end-1)-h(j-1,1:end-2)))+...
    %     Pe.*tdelt./h(j,2:end-1)./xdelt^2.*(h(j-1,3:end).*theta_w(j-1,3:end)-...
    %     2*h(j-1,2:end-1).*theta_w(j-1,2:end-1)+h(j-1,1:end-2).*theta_w(j-1,1:end-2))-...
    %     Ht.*tdelt./h(j,2:end-1).*phi(j-1,2:end-1).*h(j-1,2:end-1).*theta_w(j-1,2:end-1);
    % theta_w(j,2:end-1) = h(j-1,2:end-1).*theta_w(j-1,2:end-1)./h(j,2:end-1)+...
    %     Da.*tdelt./2./h(j,2:end-1)./xdelt^2.*((phi(j-1,3:end)).^3.*h(j-1,3:end)...
    %     .*theta_w(j-1,3:end).*(h(j-1,3:end)-h(j-1,2:end-1))-(phi(j-1,1:end-2)).^3.*...
    %     h(j-1,1:end-2).*theta_w(j-1,1:end-2).*(h(j-1,2:end-1)-h(j-1,1:end-2)))+...
    %     Da.*tdelt.*(1-phi(j-1,2:end-1)).*theta_w(j-1,2:end-1)./2./h(j,2:end-1)./xdelt^2.*...
    %     (h(j-1,3:end).*(phi(j-1,3:end)).^3.*(h(j-1,3:end)-h(j-1,2:end-1))-...
    %     h(j-1,1:end-2).*(phi(j-1,1:end-2)).^3.*(h(j-1,2:end-1)-h(j-1,1:end-2)))+...
    %     Pe.*tdelt./h(j,2:end-1)./xdelt^2.*(h(j-1,3:end).*theta_w(j-1,3:end)-...
    %     2*h(j-1,2:end-1).*theta_w(j-1,2:end-1)+h(j-1,1:end-2).*theta_w(j-1,1:end-2))-...
    %     Ht.*tdelt./h(j,2:end-1).*phi(j-1,2:end-1).*h(j-1,2:end-1).*theta_w(j-1,2:end-1);
    
    % advected quantity
    htw(j,2:end-1) = htw(j-1,2:end-1)+Da.*tdelt./2./xdelt^2.*...
        ((phi(j-1,3:end)).^3.*htw(j-1,3:end).*(h(j-1,3:end)-h(j-1,2:end-1))...
        -(phi(j-1,1:end-2)).^3.*htw(j-1,1:end-2).*(h(j-1,2:end-1)-h(j-1,1:end-2)))+...
        Pe.*tdelt./xdelt^2.*(htw(j-1,3:end)-2.*htw(j-1,2:end-1)+htw(j-1,1:end-2))-...
        Ht.*tdelt.*phi(j-1,2:end-1).*htw(j-1,2:end-1);

    % BC
    htw(j,1) = wt;
    htw(j,end) = htw(j,end-1);

    % crop
    htw(htw>wt) = wt;
    htw(htw<water_min_temp) = water_min_temp;

    % solve water temp
    theta_w(j,:) = htw(j,:)./h(j,:);

    % BC check again
    theta_w(j,1) = wt;
    theta_w(j,end) = theta_w(j,end-1);

    % crop again
    theta_w(theta_w > wt) = wt;
    theta_w(theta_w < water_min_temp) = water_min_temp;

    % ice temperature
    theta_i(j,2:end-1) = theta_i(j-1,2:end-1)+Pe.*tdelt./xdelt^2.*(theta_i(j-1,3:end)...
        -2*theta_i(j-1,2:end-1)+theta_i(j-1,1:end-2))-Ht.*tdelt.*((1-phi(j-1,2:end-1)).*...
        h(j-1,2:end-1).*theta_i(j-1,2:end-1));

    % BC
    theta_i(j,1) = theta_i(j,2);
    theta_i(j,end) = theta_i(j,end-1);

    % crop
    theta_i(theta_i > 0) = 0;
    theta_i(theta_i < ct) = ct;

    % let user know what step, temperature is running
    disp(length(t) - j)
    disp(ct)
end

end