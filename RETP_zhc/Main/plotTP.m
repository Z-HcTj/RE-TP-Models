%% Two-phase flow model
clear
clc

%% Input parameters
% Input data
id = load('Inputdata.mat');
% Input simulation time
tt = load('simulation time.mat');
% Output save
load('mydata.mat');

% Location of points
x = id.x;
z = id.z;
% Number of grids (row and column)
nx = id.nx;
nz = id.nz;
% Grid size
xgs = id.xgs;
zgs = id.zgs;
% Distance between points
dxp = id.dxp;
dzp = id.dzp;
% Tiny step 
dthew = id.dthew;
dpw = id.dpw;
% Compressible parameters
aw = id.aw;
bn = id.bn;
% Reference density
rhow0 = id.rhow0;
rhon0 = id.rhon0;
% Viscosity
vw = id.vw;
vn = id.vn;
% Specific storage
Ss = id.Ss;
% Gravity
g = id.g;
% Gracity option
gr = id.gr;
% Precision
ethew = id.ethew;
epw = id.epw;
% Soil parameters
n = id.n;
a = id.a;
m = id.m;
pcm = id.pcm;
k = id.k;
phi = id.phi;
thewr = id.thewr;
thenr = id.thenr;
thews = id.thews;
thens = id.thens;
% Initial timestep
dt = id.dt;
% Allocation
dt_max = id.dt_max;
inu_ac = id.inu_ac;
inu_limit = id.inu_limit;
% Initial simulation time
T = 0;
% Recond time
T_rec = tt.T_simulation * id.T_ran;
recond_index = 1;
% Intermediate valuables
pc = zeros(nz+2,nx+2);
pn = zeros(nz+2,nx+2);
seff = zeros(nz+2,nx+2);
b = zeros(nz*nx,1);
A = zeros(nz*nx,nz*nx);
% Initial and boundary conditions
% Pressure head boundary
pw = id.pw;
thew = id.thew;
% Reference valuables
thew0 = thew;
pw0 = pw;
pn0 = pn;
% Iteration valuables
pwn = pw;
thewn = thew;
pwnn = pw;
thewnn = thew;
inum = [];
% Boundary condition options
left_close = id.left_close;
right_close = id.right_close;
top_close = id.top_close;
bottom_close = id.bottom_close;
% Output recond valuables
thew_rec = zeros(nz+2,nx+2,length(T_rec));
pw_rec = zeros(nz+2,nx+2,length(T_rec));
pc_rec = zeros(nz+2,nx+2,length(T_rec));
pn_rec = zeros(nz+2,nx+2,length(T_rec));
qw_rec = zeros(nz+1,nx+2,length(T_rec));
qn_rec = zeros(nz+1,nx+2,length(T_rec));
un_rec = zeros(nz+1,nx+2,length(T_rec));
% Air entry effect
FBAE_option = id.FBAE_option;
pnae = id.pnae;
thewae = id.thetanae;

% Sink/Source term to flowrate boundary options (mass flux)
SStoFD = id.SStoFD;
if SStoFD == 1
    q = id.q;
    SSq = q / dzp(1,2);
    SSq_left = id.SSq_left;
    SSq_right = id.SSq_right;
end

% Change upper boundary to input boundary value
if SStoFD == 1
    UP_stage1_option = SStoFD;
    UP_stage2_option = 0;
    UP_stage3_option = 0;
    pwbmi = pw(2,SSq_left+1:SSq_right+1);
    pwbmm = -pw(2,SSq_left+1:SSq_right+1);
    epw_BM = id.epw_BM;
    curan = SSq_left+1:SSq_right+1;
    pw(1,curan) = 0;
end

% Fixed boundary condition
pn(pw >= 0) = pw(pw >= 0);    

% Bubbling time
airbub_time = [];
airclose_time = [];

% Mass balance valuables
dqw = 0;
qwn = 0;
qwdet = 0;
dqn = 0;
qnn = 0;
unn = 0;

% Ponding time recond
T_tp_pond = 0;

% Valuable for recording whether 
% a time point is incorrect or not
outtime = [];

% Operation time
tic;

% MVG model application
MVG_option = 0;
hc_mvg = id.hc_mvg;

%% Simulation


    % Time increasing
    T = T + dt;
    % Iteration times
    inu = 0;
    % Boundary Ponding setup 
    if SStoFD == 1
        if UP_stage2_option == 1
            pw_up2_ini = pw(1,curan);
        end
        if UP_stage3_option == 1
            pw_up3_ini = pw(1,curan);
        end
    end


        % Flowrate top boundary setup
        if SStoFD == 1
            % UP_stage1 -- no ponding occurs, update top boundary
            if UP_stage1_option == 1
                [thewtop, pwtop] = getupperboundaryBM(rhow0,pwbmi,pwbmm,pw(2,curan)...
                    ,thew(2,curan),q,dzp(1,curan),n(1,curan),m(1,curan),k(1,curan),vw,a(1,curan),...
                    g,gr,thewr(1,curan),thenr(1,curan),phi(1,curan),epw_BM);
                thew(1,curan) = thewtop;
                pw(1,curan) = pwtop;
                if pwtop > 0
                    UP_stage1_option = 0;
                    UP_stage2_option = 1;
                    pw(1,curan) = 0;
                    pw_up2_ini= pw(1,curan);
                    T_tp_pond = T;
                end
            end
            % UP_stage2 -- ponding occurs, update top boundary
            if UP_stage2_option == 1
                thew(1,curan) = phi(1,curan) - thenr(1,curan);
                seff2_up2 = (thew(2,curan) - thewr(2,curan)) ...
                    ./ (phi(2,curan) - thewr(2,curan) - thenr(2,curan));
                krw2_up2 = seff2_up2 .^ 0.5 .* (1 ...
                    - (1 - (seff2_up2 .^ (1 ./ m(2,curan)))) .^ m(2,curan)) .^ 2;
                K_up2 = (k(1,curan) * 1 / vw + k(2,curan) * krw2_up2 / vw) / 2;
                pw(1,curan) = pw_up2_ini;
                pw_dt_up2 = (-q / rhow0 ./ K_up2 .* dzp(1,curan) - pw(1,curan)...
                    + pw(2,curan)  - rhow0 * g * dzp(1,curan)) ...
                    / (1 + dzp(1,curan) ./ (rhow0 * g * dt * K_up2));
                pw(1,curan) = pw_up2_ini + pw_dt_up2;
                pc2_up2 = pcm(2,curan) .* (seff(2,curan) ...
                    .^ (-1 ./ m(2,curan)) -1) .^ (1 ./ n(2,curan));
                pn2_up2 = pw(2,curan) + pc2_up2;
                if pw(1,curan) < 0
                    UP_stage1_option = 1;
                    UP_stage2_option = 0;
                end
                if pn2_up2 >= (pnae + pw(1,curan))
                    UP_stage2_option = 0;
                    UP_stage3_option = 1;
                    pw_up3_ini = pw(1,curan);
                end
            end
            % UP_stage3 -- ponding occurs, air breakthrough
            if UP_stage3_option == 1
                thew(1,curan) = thewae(1,curan);
                seff1_up3 = (thew(1,curan) - thewr(1,curan)) ...
                    ./ (phi(1,curan) - thewr(1,curan) - thenr(1,curan));
                seff2_up3 = (thew(2,curan) - thewr(2,curan)) ...
                    ./ (phi(2,curan) - thewr(2,curan) - thenr(2,curan));
                krw1_up3 = seff1_up3 .^ 0.5 .* (1 ...
                    - (1 - (seff1_up3 .^ (1 ./ m(1,curan)))) .^ m(1,curan)) .^ 2;
                krw2_up3 = seff2_up3 .^ 0.5 .* (1 ...
                    - (1 - (seff2_up3 .^ (1 ./ m(2,curan)))) .^ m(2,curan)) .^ 2;
                K1_up3 = k(1,curan) * krw1_up3 / vw;
                K2_up3 = k(2,curan) * krw2_up3 / vw;
                K_up3 = (K1_up3 + K2_up3) / 2;
                pw(1,curan) = pw_up3_ini;
                pw_dt_up3 = (-q / rhow0 ./ K_up3 .* dzp(1,curan) - pw(1,curan) ...
                    + pw(2,curan)  - rhow0 * g * dzp(1,curan)) ...
                    / (1 + dzp(1,curan) ./ (rhow0 * g * dt * K_up3));
                pw(1,curan) = pw_up3_ini + pw_dt_up3;
            end
        end

        % Air breakthrough (not flowrate boundary)
        if FBAE_option == 1
            pc2_fb2 = pcm(2,2) .* (seff(2,2) ...
                .^ (-1 ./ m(2,2)) -1) .^ (1 ./ n(2,2));
            pn2_fb2 = pw(2,2) + pc2_fb2;
            if pn2_fb2 >= (pnae + pw(1,2))
                FBtopchange = 1;
                thew(1,2) = thewae;
                airbub_time = [airbub_time; toc];
            elseif pn2_fb2 < (pnae + pw(1,2))
                FBtopchange = 0;
                thew(1,2) = id.thew(1,2);
                airclose_time = [airclose_time; toc];
            end
        end

        % Iteration increasing
        inu = inu + 1;

        % Effective saturation & relative permeability
        seff = [0:0.01:1.01]';
        seffn = (thewn - thewr) ./ (phi - thewr - thenr);
        seff(seff > 1) = 1;
        seffn(seffn > 1) = 1;
        [~,krw] = TP_VGkrpc(MVG_option,seff,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);
        [~,~,krn] = TP_VGkrpc(MVG_option,seff,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);

        % Capillary pressure
        pc = TP_VGkrpc(MVG_option,seff,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);
        pcn = pcm .* (seffn .^ (-1 ./ m) -1) .^ (1 ./ n);
