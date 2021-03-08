%%%%%%%%%%%%%%%%%%%%%%%%%%%   NP FINAL MODEL   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%    NESTOR PADRÓ    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%      PARAMETERS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.time = 0:365
points = 200; %total number of grids to analyze
%%%%plankton constants
%P = plankton concentration
param.Z = 200; %m         Total depthat
param.D = 5*8.64; %m2/s       Diffusivity rate
param.U = 0.04*24; %m/day      Sinking velocity
param.deltaZ = param.Z/points; %m     depth of the sections grid
param.z = ((param.deltaZ)/2) : param.deltaZ : (param.Z-((param.deltaZ)/2)) ;    %declare segments of depth each deltaZ
param.n = length(param.z); %No. of grid cells (cuantos segmentos hay)

%%%%light constants
%I = Light intensity
param.kbg = 0.045; % background turbidity-- (1 / m)
param.kp = 6e-12; % specific light attenuation of phytoplankton-- (m2 / cell)
param.l = 0.01*24; % specific loss rate-- (1 / day)
param.mumax = 0.04*24; % maximal specific production rate-- (1 / day)
param.I0 = 450*86400; % incident light intensity (initial)-- (micromol photons / m2 day)

param.H_I = 20*86400; %half-saturation constant of light-limited growth [μmol photons m-2 day-1]
param.H_N = 0.01; %Of nutrients

%%%%nutrients constants
%N = Nutrient concentration
param.alpha = 1*10^-9; % Nutrient content of a phytoplankton cell-- (mmol nutrient/cell)
param.eps = 0.005; % Nutrient recycling coefficient


%%%%%%%%%%%%%%%%%%%%%%       INITIAL CONDITIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.P0=ones(1, points)*1e06;

param.N0=ones(1, points)*100;


%%%%%%%%%%%%%%%%%%%%%%%%      FUNCTION CALL      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(param.z);
[t, C] = vertical(param);
N = C(:, 1:n);
P = C(:, n+1:2*n);
Plast=P(end,:);
plot(Plast, param.z,'r-');xlim([0 4e9]);ylim([0 200]);
axis ij
hold on


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      PLOTS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(param.z);
tiledlayout (3,1)
figure
nexttile
surface(t,param.z',N')
axis ij; shading interp; colorbar;
title('Nutrients     [mmol nutrient m-3]');

nexttile
surface(t,param.z',P')
axis ij; shading interp; colorbar;
ylabel('Depth (meters)'); title('Plankton         [cells m-3]');

%nexttile
% single light plot
%plot(lightcalc(param.z,C,param.deltaZ,param.kp,param.kbg,param.I0), param.z, 'Linewidth', 1.5)
%axis ij
%xlabel('Light (μmol photons m-2 day-1)'); title('Light');

%%%% LIGHT INTENSITY THROUGH TIME %%%%
lightvalues = zeros(length(param.time),n); %matrix of light intensity

for a = 1:length(param.time) % a for a loop in every time-step
        
lightvalues(a,:) = lightcalc(param.z,P(a,:),param.deltaZ,param.kp,param.kbg,param.I0); %vector of light intensity in each timestep "a" is assigned to a matrix

end
nexttile
hold on
surface(t,param.z,lightvalues')
axis ij; shading interp; colorbar;
xlabel('Time [days]'); title('Light intensity          [Mmol photons m-2 d-1]');
hold off

%%
%%%% SENSITIVITY %%%%%%
hnnew=[0.01 0.1 1 10 100] %different half-saturation constants to analyze
maxvalue = zeros;
depthblock = zeros;
lightsens = zeros(length(hnnew),n);
limNsens = zeros(length(hnnew),n); %sens to avoid overwriting the original limN/L
limLsens = zeros(length(hnnew),n);

for v = 1:length(hnnew) %how many times we will analyze (how many constants to analyze)
param.H_N = hnnew(v); %define our parameter with the constant to analyze at every loop
[t, C] = vertical(param); 
N = C(:, 1:n);
P = C(:, n+1:2*n);
Nlast = N(end,:);
Plast = P(end,:); %Plast will be the last vertical profile of plankton concentrations
[maxvalue(v), depthblock(v)]=max(Plast); %extract the max value and position of this last P, in every loop i

lightsens (v,:) = lightcalc(param.z,Plast,param.deltaZ,param.kp,param.kbg,param.I0);
limNsens (v,:) = Nlast ./ (param.H_N + Nlast);
limLsens (v,:) = lightsens(v,:) ./ (param.H_I + lightsens(v,:));
end
figure
tiledlayout(3,1)
nexttile
plot(hnnew,maxvalue,'LineWidth',1.5)
hold on
plot(hnnew,maxvalue,'*')
hold off
ylabel('Max P concentration');
nexttile
semilogx(hnnew,maxvalue,'LineWidth',1.5)
hold on
semilogx(hnnew,maxvalue,'*')
hold off
ylabel('Max P concentration');
nexttile
plot(hnnew,depthblock,'LineWidth',1.5)
hold on
plot(hnnew,depthblock,'*')
ylabel('Depth of Max P concentration');xlabel('half-saturation constant of nutrients');

%limiting factors for different constant analyzed
figure
plot (limNsens,param.z,'DisplayName','Limiting Nutrients');
legend('0.01','0.1','1','10','100','Location','Southeast');
hold on
plot (limLsens,param.z,'g','DisplayName','Limiting Light');xlabel('Limiting factor'); ylabel('Depth (meters)');
axis ij
title('Limiting factors for different Nutrient half-saturation constant [[μmol nutrient m-2 day-1]');




%%
%%%%%%% LIMITING FACTOR %%%%%%
limN = zeros(length(param.time),n); 
limL = zeros(length(param.time),n);
limitant = zeros(length(param.time),n);
for ti = 1:length(param.time)
    for j = 1:n
        %for each timestep "ti" we calculate the limitant of Nutrient and
        %light in each grid "j"
limN(ti,j) = (N(ti,j)/ (param.H_N + N(ti,j)));
limL(ti,j) = (lightvalues(ti,j) / (param.H_I + lightvalues(ti,j)));
% we compare them to see which is actually limitating
if limN(ti,j) < limL(ti,j)
    limitant(ti,j) = 0; %light is limiting
else
    limitant(ti,j) = 1; %nutrient is limiting
end

    end
end
figure;
hold on
surface(t,param.z,limitant')
axis ij; colorbar
xlabel('Time [days]'); ylabel('Depth [m]');title('Limiting factor');



%%
figure
%%%%%%PLOT OF DAYS COMPARISON (VERTICAL PROFILES)
%%%%%%  1 day
tiledlayout(3,3)
nexttile % PLANKTON
plot(P(1,:), param.z, 'Linewidth', 1.5); %xlim([0 7e6]);
axis ij; title('Plankton [cells/m3]');

nexttile % NUTRIENTS
plot(N(1,:), param.z, 'Linewidth', 1.5); %xlim([0 100]);
axis ij; xlabel('t = 0'); title('Nutrients [mmol nutrient/m3]');

nexttile % LIGHT
plot(lightcalc(param.z,C(1,:),param.deltaZ,param.kp,param.kbg,param.I0), param.z, 'Linewidth', 1.5)
axis ij; title('Light [μmol photons m-2 day-1]');

%%%%%% 80 days  (bloom?)
nexttile % PLANKTON
plot(P(80,:), param.z, 'Linewidth', 1.5); %xlim([0 7e6]);
axis ij; ylabel( 'Depth [m]'); 

nexttile % NUTRIENTS
plot(N(80,:), param.z, 'Linewidth', 1.5); %xlim([0 100]);
axis ij; xlabel('t = 80');

nexttile % LIGHT
plot(lightcalc(param.z,C(80,:),param.deltaZ,param.kp,param.kbg,param.I0), param.z, 'Linewidth', 1.5)
axis ij;

%%%%%% 300 days
nexttile % PLANKTON
plot(P(end,:), param.z, 'Linewidth', 1.5); %xlim([0 7e6]);
axis ij; 

nexttile % NUTRIENTS
plot(N(end,:), param.z, 'Linewidth', 1.5); %xlim([0 100]);
axis ij; xlabel('t = 300'); 

nexttile % LIGHT
plot(lightcalc(param.z,C(end,:),param.deltaZ,param.kp,param.kbg,param.I0), param.z, 'Linewidth', 1.5)
axis ij;
%

%%
%%%% DIFFERENT CELL GRID SIZES %%%%%%
points = 200;
n = length(param.z);
[t, C] = vertical(param);
N = C(:, 1:n);
P = C(:, n+1:2*n);
Plast=P(end,:);
figure
hold on
plot(Plast, param.z,'r-','DisplayName','200 grid cells');xlim([0 4e9]);ylim([0 200]);
hold on
axis ij


points = 300;
n = length(param.z);
[t, C] = vertical(param);
N = C(:, 1:n);
P = C(:, n+1:2*n);
Plast=P(end,:);
plot(Plast, param.z,'b.','DisplayName','300 grid cells');xlim([0 4e9]);ylim([0 200]);
axis ij
hold on

points = 1000;
n = length(param.z);
[t, C] = vertical(param);
N = C(:, 1:n);
P = C(:, n+1:2*n);
Plast=P(end,:);
plot(Plast, param.z,'g*','DisplayName','1000 grid cells');xlim([0 4e9]);ylim([0 200]);
xlabel('Phytoplankton [cells m-3]'); ylabel('Depth [m]');
legend('Location','Southeast'); 
title('Plankton concentrations at last time step');
axis ij
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, C] = vertical(param)

n = length(param.z); %no. of grid cells (cuantos segmentos hacia abajo se van midiendo)

[t,C] = ode45(@derivative, param.time,[param.N0 param.P0]);



function dCdt = derivative(t,C) %depends on time and y-concentrations

N = C(1 : n);     %concentration is one vector, half of nutrients
P = C(n+1 : 2*n); %half of plankton, %is it end? I also tried with 2*param.n

iaux = 1:n;
iauxd = 2:n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  FLUXES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% NUTRIENT FLUX %%%%%%%%%%
JNdiff(iauxd) = -param.D * (N(iauxd) - N(iauxd-1)) / param.deltaZ; 
JNdiff(1) = 0; %No flux of nutrients at surface
JNdiff(n+1) = -param.D * (100 - N(end))/ param.deltaZ; ; %or bottom


%%%%%%%% PHYTOPLANKTON FLUX %%%%%%%%
JPadv(iauxd) = param.U * P(iauxd-1);
JPadv(1) = 0; %No flux of nutrients at surface, surface of P is y(101)
JPadv(n+1) = 0; %or bottom y(201)

JPdiff(iauxd) = -param.D * (P(iauxd) - P(iauxd-1)) / param.deltaZ;
JPdiff(1) = 0; %No flux of nutrients at surface
JPdiff(n+1) = 0;

    %total
JPtot = JPadv + JPdiff;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  LIGHT EQUATIONS %%%%%%%%%%%%%%%%%%%%%%%%%

%I = lightcalc(param.z,C,param.deltaZ,param.kp,param.kbg,param.I0)*(max(0,(0.5-sin(2*pi*(t/180))))); %t/#days to change....0.5-sin because we want light goes from -0.5 to 0.5, and the max is to avoid negative light, just goes from 0 to 0.5
I = lightcalc(param.z,C,param.deltaZ,param.kp,param.kbg,param.I0); %constant light
mu = zeros;
mu(iaux) = param.mumax * min((N(iaux)'./(param.H_N+N(iaux)')) , (I(iaux)./(param.H_I+I(iaux))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  FINAL EQUATIONS %%%%%%%%%%%%%%%%%%%%%%%%%
%        Nutrients equation---------
uptake=zeros;  %just to allocate the chain and make the code faster
recycle=zeros;
uptake = param.alpha .* mu(iaux) .*P(iaux)';
recycle = param.eps .* param.alpha .* param.l .* P(iaux)';


dNdt=zeros;
dNdt = - uptake + recycle  - (JNdiff(2:(param.n+1))-JNdiff(1:param.n)) / param.deltaZ;

%        Phytoplankton--------
dPdt=zeros;
dPdt = (-(JPtot(2:(param.n+1))-JPtot(1:param.n)) / param.deltaZ) + mu .* P' - param.l .* P';


%dNdt = dNdt';
dCdt = [dNdt dPdt]';

end
end

