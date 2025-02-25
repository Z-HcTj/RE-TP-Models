%% Input data
clear
clc

% Key soil param change
load('Case_name.mat')
load('run_ii.mat')

% Project root directory
oper_filepath = fullfile(fileparts(fileparts(fileparts(fileparts(pwd)))), 'Main');
addpath(oper_filepath)

%% Domain 
% Location of points
x = [1];
z = [-0.005:-0.01:-0.995]';
% Number of grids (row and column)
nx = length(x);
nz = length(z);
% Grid size
xgs = 1 * ones(length(z),length(x));
zgs = 0.01 * ones(length(z),length(x));
% Distance between points
dxp = 1 * ones(length(z) + 2,length(x) + 1);
dzp = 0.01 * ones(length(z) + 1,length(x) + 2);
%% Perturbation
% Perturbation increment
dthew = -1 * 10 ^ -11;
dpw = -1 * 10 ^ -8;
%% Fluid properties
% Compressible coefficients
aw = 0;
bn = 1.189e-5;
% Reference density
rhow0 = 1000;
rhon0 = 1.29;
% Viscosity
vw = 1 * 10^-3;
vn = 1.8 * 10^-5;
% Specific storage
Ss = 1 * 10 ^ -5;
%% Gravity
% Gravity
g = 9.81;
% Gravity option 
gr = 1;
%% Precision
% Convergence threshold
ethew = 1 * 10 ^ -3;
epw = 1;
%% Soil
% Van Genuchten soil parameters
n = zeros(nz+2,nx+2);
a = zeros(nz+2,nx+2);
% Capillary parameter
pcm = zeros(nz+2,nx+2);
% Porosity
phi = zeros(nz+2,nx+2);
% Absolute permeability
k = zeros(nz+2,nx+2);
% Residual & saturated water / air content
thewr = zeros(nz+2,nx+2);
thews = zeros(nz+2,nx+2);
thenr = zeros(nz+2,nx+2);
thens = zeros(nz+2,nx+2);
% Value
n(:,:) = 4.4;
m = 1 - 1 ./ n;
a(:,:) = 2.434;
a(42:end,:) = key_param(ii);
pcm = rhow0 * g ./ a;
k(:,:) = 396 * 0.01 / 86400 / (rhow0 * g / vw);
phi(:,:) = 0.43;
thewr(:,:) = 0;
thenr(:,:) = 0;
thews(:,:) = phi - thenr;
thens(:,:) = phi - thewr;
%% Timestep & Simulation time
% Initial timestep size
dt = 1e-3;
% Timestep control parameters
dt_max = 10;
inu_ac = 4;
inu_limit = 10;
% Recond time points 
T_ran = [0:0.01:1]';
%% Initial & Boundary condition
thew = zeros(nz+2,nx+2);
pw = zeros(nz+2,nx+2);

% Boundary condition
pw_bd = 0 * rhow0 * g;

% Initial condition
pw_ini = [-1.8:0.01:-0.8]' * rhow0 * g;

% Combination
pw(:,1) = [pw_bd; pw_ini];
pw(:,2) = [pw_bd; pw_ini];
pw(:,3) = [pw_bd; pw_ini];
thew = pwtranthewVG(pw,phi,thewr,thenr,n,m,a,rhow0,g);
hw = pw / (rhow0 * g);

% Air breakthrough (1-open & 0-close)
FBAE_option = 1;

%% Flowrate Boundary
% Flowrate option (1-open & 0-close)
SStoFD = 0;
% Mass flux
q = -k(2,1) * (rhow0 * g / vw)  * 0.8 * rhow0; 
% Range
SSq_left = 1;
SSq_right = 1;
% Boundary water content change
UP_stage1_option = SStoFD;
UP_stage2_option = 0;
UP_stage3_option = 0;
pnae = (m(2,1) .^ (1 ./ n(2,1))) ./ a(2,1) * rhow0 * g;
thetanae = (1 ./ (1 + m(2,1))) .^ m(2,1) .* (phi(2,1) - thewr(2,1) - thenr(2,1)) + thewr(2,1);
epw_BM = 1e-3;

%% Boundary open/close options
% 1-close & 0-open
left_close = 1;
right_close = 1;
top_close = 0;
bottom_close = 0;
% Special setup for 2D heter. case
top_wclose_nopen = 0;

%% Depth of the wetting front when stop (RE model)
dep = 81;

%% Save
% Save input data
save(fullfile(oper_filepath,"Inputdata.mat"),"bottom_close","top_close","right_close","left_close"...
    ,"thew","pw","inu_ac","dt_max","dt","thens","thews","thenr","thewr"...
    ,"phi","k","pcm","m","a","n","epw","ethew","gr","g","Ss","vn","vw","rhon0"...
    ,"rhow0","bn","aw","dpw","dthew","dzp","dxp","zgs","xgs","nz","nx","z","x"...
    ,"T_ran","SStoFD","q","SSq_right","SSq_left","pnae","epw_BM"...
    ,"inu_limit","dep","thetanae","FBAE_option","top_wclose_nopen")
