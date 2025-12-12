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
masdqw = 0;
masdqn = 0;

% Ponding time recond
T_tp_pond = 0;

% Valuable for recording whether 
% a time point is incorrect or not
outtime = [];

% Operation time
tic;

% back parameter
inu_back = 0;
tp_outline = 0;

% MVG model application (not suitable for TP model)
MVG_option = 0;
hc_mvg = id.hc_mvg;

%% Simulation
while 1

    % Time increasing
    T = T + dt;
    % Iteration times
    inu = 0;
    % back paramater
    inu_back = 0;
    % Boundary Ponding setup 
    if SStoFD == 1
        if UP_stage2_option == 1
            pw_up2_ini = pw(1,curan);
        end
        if UP_stage3_option == 1
            pw_up3_ini = pw(1,curan);
        end
    end

    while 1
        % Flowrate top boundary setup
        if SStoFD == 1
            % UP_stage1 -- no ponding occurs, update top boundary
            if UP_stage1_option == 1
                [thewtop, pwtop] = getupperboundaryBM(MVG_option,hc_mvg,rhow0,pwbmi,pwbmm,pw(2,curan)...
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
        seff = (thew - thewr) ./ (phi - thewr - thenr);
        seffn = (thewn - thewr) ./ (phi - thewr - thenr);
        seff(seff > 1) = 1;
        seffn(seffn > 1) = 1;
        [~,krw] = TP_VGkrpc(MVG_option,seff,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);
        [~,~,krn] = TP_VGkrpc(MVG_option,seff,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);

        % Capillary pressure
        pc = TP_VGkrpc(MVG_option,seff,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);
        pcn = TP_VGkrpc(MVG_option,seffn,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);

        % Density
        rhow = rhow0 * exp(aw * (pw - pw0));
        rhon = rhon0 + bn * (pw + pc - pn0);
        rhown = rhow0 * exp(aw * (pwn - pw0));
        rhonn = rhon0 + bn * (pwn + pcn - pn0);

        % Permeability (intercell)
        Kw = k .* krw .* rhow / vw;
        Kn = k .* krn .* rhon / vn;
        % Arithmetic
        Kwxm = (Kw(2:end-1,1:end-1) + Kw(2:end-1,2:end)) / 2;
        Knxm = (Kn(2:end-1,1:end-1) + Kn(2:end-1,2:end)) / 2;
        Kwzm = (Kw(1:end-1,2:end-1) + Kw(2:end,2:end-1)) / 2;
        Knzm = (Kn(1:end-1,2:end-1) + Kn(2:end,2:end-1)) / 2;

        % Prepare valuables for matrix A dthew
        thew_dt = thew + dthew;
        seff_dt = (thew_dt - thewr) ./ (phi - thewr - thenr);
        seff_dt(seff_dt > 1) = 1;
        [~,krw_dt] = TP_VGkrpc(MVG_option,seff_dt,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);
        [~,~,krn_dt] = TP_VGkrpc(MVG_option,seff_dt,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);
        pc_dt = TP_VGkrpc(MVG_option,seff_dt,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);
        rhon_dt = rhon0 + bn * (pw + pc_dt - pn0);

        Kw_dt = k .* krw_dt .* rhow / vw;
        Kn_dt = k .* krn_dt .* rhon_dt / vn;
        Kwxm_dt_iaj = (Kw(2:end-1,1:end-1) + Kw_dt(2:end-1,2:end)) / 2;
        Kwxm_dt_ibj = (Kw_dt(2:end-1,1:end-1) + Kw(2:end-1,2:end)) / 2;
        Knxm_dt_iaj = (Kn(2:end-1,1:end-1) + Kn_dt(2:end-1,2:end)) / 2;
        Knxm_dt_ibj = (Kn_dt(2:end-1,1:end-1) + Kn(2:end-1,2:end)) / 2;
        Kwzm_dt_ija = (Kw(1:end-1,2:end-1) + Kw_dt(2:end,2:end-1)) / 2;
        Kwzm_dt_ijb = (Kw_dt(1:end-1,2:end-1) + Kw(2:end,2:end-1)) / 2;
        Knzm_dt_ija = (Kn(1:end-1,2:end-1) + Kn_dt(2:end,2:end-1)) / 2;
        Knzm_dt_ijb = (Kn_dt(1:end-1,2:end-1) + Kn(2:end,2:end-1)) / 2;

        % Prepare valuables for matrix A dpw
        pw_dp = pw + dpw;
        rhow_dp = rhow0 * exp(aw * (pw_dp - pw0));
        rhon_dp = rhon0 + bn * (pw_dp + pc - pn0);
        Kw_dp = k .* krw .* rhow_dp / vw;
        Kn_dp = k .* krn .* rhon_dp / vn;
        Kwxm_dp_iaj = (Kw(2:end-1,1:end-1) + Kw_dp(2:end-1,2:end)) / 2;
        Kwxm_dp_ibj = (Kw_dp(2:end-1,1:end-1) + Kw(2:end-1,2:end)) / 2;
        Knxm_dp_iaj = (Kn(2:end-1,1:end-1) + Kn_dp(2:end-1,2:end)) / 2;
        Knxm_dp_ibj = (Kn_dp(2:end-1,1:end-1) + Kn(2:end-1,2:end)) / 2;
        Kwzm_dp_ija = (Kw(1:end-1,2:end-1) + Kw_dp(2:end,2:end-1)) / 2;
        Kwzm_dp_ijb = (Kw_dp(1:end-1,2:end-1) + Kw(2:end,2:end-1)) / 2;
        Knzm_dp_ija = (Kn(1:end-1,2:end-1) + Kn_dp(2:end,2:end-1)) / 2;
        Knzm_dp_ijb = (Kn_dp(1:end-1,2:end-1) + Kn(2:end,2:end-1)) / 2;
        % Left close
        if left_close == 1
            Kwxm(:,1) = 0;
            Knxm(:,1) = 0;
            Kwxm_dt_iaj(:,1) = 0;
            Kwxm_dt_ibj(:,1) = 0;
            Knxm_dt_iaj(:,1) = 0;
            Knxm_dt_ibj(:,1) = 0;
            Kwxm_dp_iaj(:,1) = 0;
            Kwxm_dp_ibj(:,1) = 0;
            Knxm_dp_iaj(:,1) = 0;
            Knxm_dp_ibj(:,1) = 0;
        end
        % Right close
        if right_close == 1
            Kwxm(:,end) = 0;
            Knxm(:,end) = 0;
            Kwxm_dt_iaj(:,end) = 0;
            Kwxm_dt_ibj(:,end) = 0;
            Knxm_dt_iaj(:,end) = 0;
            Knxm_dt_ibj(:,end) = 0;
            Kwxm_dp_iaj(:,end) = 0;
            Kwxm_dp_ibj(:,end) = 0;
            Knxm_dp_iaj(:,end) = 0;
            Knxm_dp_ibj(:,end) = 0;
        end
        % Bottom
        if bottom_close == 1
            Kwzm(end,:) = 0;
            Knzm(end,:) = 0;
            Kwzm_dt_ija(end,:) = 0;
            Kwzm_dt_ijb(end,:) = 0;
            Knzm_dt_ija(end,:) = 0;
            Knzm_dt_ijb(end,:) = 0;
            Kwzm_dp_ija(end,:) = 0;
            Kwzm_dp_ijb(end,:) = 0;
            Knzm_dp_ija(end,:) = 0;
            Knzm_dp_ijb(end,:) = 0;
        end
        % Upper
        if top_close == 1
            Kwzm(1,:) = 0;
            Knzm(1,:) = 0;
            Kwzm_dt_ija(1,:) = 0;
            Kwzm_dt_ijb(1,:) = 0;
            Knzm_dt_ija(1,:) = 0;
            Knzm_dt_ijb(1,:) = 0;
            Kwzm_dp_ija(1,:) = 0;
            Kwzm_dp_ijb(1,:) = 0;
            Knzm_dp_ija(1,:) = 0;
            Knzm_dp_ijb(1,:) = 0;
        end

        % Special setup for 2D heter. case
        if id.top_wclose_nopen == 1
            Kwzm(1,:) = 0;
            Kwzm_dt_ija(1,:) = 0;
            Kwzm_dt_ijb(1,:) = 0;
            Kwzm_dp_ija(1,:) = 0;
            Kwzm_dp_ijb(1,:) = 0;
        end
        
        % Source/term to flowrate condition (for upper boundary)
        if SStoFD == 1
            if UP_stage1_option == 1
                Kwzm(1,SSq_left:SSq_right) = 0;
                Kwzm_dt_ija(1,SSq_left:SSq_right) = 0;
                Kwzm_dt_ijb(1,SSq_left:SSq_right) = 0;
                Kwzm_dp_ija(1,SSq_left:SSq_right) = 0;
                Kwzm_dp_ijb(1,SSq_left:SSq_right) = 0;
            end
            if UP_stage2_option == 1
                Knzm(1,SSq_left:SSq_right) = 0;
                Knzm_dt_ija(1,SSq_left:SSq_right) = 0;
                Knzm_dt_ijb(1,SSq_left:SSq_right) = 0;
                Knzm_dp_ija(1,SSq_left:SSq_right) = 0;
                Knzm_dp_ijb(1,SSq_left:SSq_right) = 0;
            end
        end

        % Saturated no air flow out
        maskkn = thew(1, 2:end-1) >= thews(1, 2:end-1);
        Knzm(1,maskkn) = 0;
        Knzm_dt_ija(1,maskkn) = 0;
        Knzm_dt_ijb(1,maskkn) = 0;
        Knzm_dp_ija(1,maskkn) = 0;
        Knzm_dp_ijb(1,maskkn) = 0;
            

        %% Vector b
        % Rw
        Rwp1 = (rhow(2:end-1,2:end-1) .* thew(2:end-1,2:end-1) - ...
            rhown(2:end-1,2:end-1) .* thewn(2:end-1,2:end-1)) / dt;
        Rwp2 = Ss ./ phi(2:end-1,2:end-1) .* thew(2:end-1,2:end-1) ...
            .* rhow(2:end-1,2:end-1) ...
            .* (pw(2:end-1,2:end-1) - pwn(2:end-1,2:end-1)) ...
            / (dt * g * rhow0);
        Rwp3_i = 1 ./ xgs .* (-Kwxm(:,2:end) ...
            .* (pw(2:end-1,3:end) - pw(2:end-1,2:end-1)) ./ dxp(2:end-1,2:end)...
            + Kwxm(:,1:end-1) ...
            .* (pw(2:end-1,2:end-1) - pw(2:end-1,1:end-2)) ./ dxp(2:end-1,1:end-1));
        Rwp3_j = 1 ./ zgs .* (-Kwzm(2:end,:) ...
            .* ((pw(3:end,2:end-1) - pw(2:end-1,2:end-1)) ./ dzp(2:end,2:end-1)...
            - (rhow(3:end,2:end-1) + rhow(2:end-1,2:end-1)) / 2 * g * gr) ...
            + Kwzm(1:end-1,:) ...
            .* ((pw(2:end-1,2:end-1) - pw(1:end-2,2:end-1)) ./ dzp(1:end-1,2:end-1) ...
            - (rhow(2:end-1,2:end-1) + rhow(1:end-2,2:end-1)) / 2 * g * gr));


        Rw = Rwp1 + Rwp2 + Rwp3_i + Rwp3_j;

        % Rn
        Rnp1 = (rhon(2:end-1,2:end-1) .* (phi(2:end-1,2:end-1) - thew(2:end-1,2:end-1)) ...
            - rhonn(2:end-1,2:end-1) .* (phi(2:end-1,2:end-1) - thewn(2:end-1,2:end-1))) ...
            / dt;
        Rnp3_i = 1 ./ xgs .* (-Knxm(:,2:end) ...
            .* (pw(2:end-1,3:end) - pw(2:end-1,2:end-1) + pc(2:end-1,3:end) - pc(2:end-1,2:end-1))...
            ./ dxp(2:end-1,2:end)...
            + Knxm(:,1:end-1) ...
            .* (pw(2:end-1,2:end-1) - pw(2:end-1,1:end-2) + pc(2:end-1,2:end-1) - pc(2:end-1,1:end-2)) ...
            ./ dxp(2:end-1,1:end-1));
        Rnp3_j = 1 ./ zgs .* (-Knzm(2:end,:) ...
            .* ((pw(3:end,2:end-1) - pw(2:end-1,2:end-1) + pc(3:end,2:end-1) - pc(2:end-1,2:end-1))...
            ./ dzp(2:end,2:end-1)...
            - (rhon(3:end,2:end-1) + rhon(2:end-1,2:end-1)) / 2 * g * gr) ...
            + Knzm(1:end-1,:) ...
            .* ((pw(2:end-1,2:end-1) - pw(1:end-2,2:end-1) + pc(2:end-1,2:end-1) - pc(1:end-2,2:end-1)) ...
            ./ dzp(1:end-1,2:end-1) ...
            - (rhon(2:end-1,2:end-1) + rhon(1:end-2,2:end-1)) / 2 * g * gr));
        Rn = Rnp1 + Rnp3_i +Rnp3_j;

        %% Matrix A dthew
        % dRwdthew i,j
        dRwdt_ij_p1 = (rhow(2:end-1,2:end-1) .* thew_dt(2:end-1,2:end-1) - ...
            rhown(2:end-1,2:end-1) .* thewn(2:end-1,2:end-1)) / dt;
        dRwdt_ij_p2 = Ss ./ phi(2:end-1,2:end-1) .* thew_dt(2:end-1,2:end-1) ...
            .* rhow(2:end-1,2:end-1) ...
            .* (pw(2:end-1,2:end-1) - pwn(2:end-1,2:end-1)) ...
            / (dt * g * rhow0);
        dRwdt_ij_p3 = 1 ./ xgs .* (-Kwxm_dt_ibj(:,2:end) ...
            .* (pw(2:end-1,3:end) - pw(2:end-1,2:end-1)) ./ dxp(2:end-1,2:end)...
            + Kwxm_dt_iaj(:,1:end-1) ...
            .* (pw(2:end-1,2:end-1) - pw(2:end-1,1:end-2)) ./ dxp(2:end-1,1:end-1))...
            + 1 ./ zgs .* (-Kwzm_dt_ijb(2:end,:) ...
            .* ((pw(3:end,2:end-1) - pw(2:end-1,2:end-1)) ./ dzp(2:end,2:end-1)...
            - (rhow(3:end,2:end-1) + rhow(2:end-1,2:end-1)) / 2 * g * gr) ...
            + Kwzm_dt_ija(1:end-1,:) ...
            .* ((pw(2:end-1,2:end-1) - pw(1:end-2,2:end-1)) ./ dzp(1:end-1,2:end-1) ...
            - (rhow(2:end-1,2:end-1) + rhow(1:end-2,2:end-1)) / 2 * g * gr));
        dRwdt_ij = ((dRwdt_ij_p1 + dRwdt_ij_p2 + dRwdt_ij_p3) - Rw) / dthew;

        % dRwdthew i-1,j
        dRwdt_ibj_p3 = 1 ./ xgs .* (-Kwxm(:,2:end) ...
            .* (pw(2:end-1,3:end) - pw(2:end-1,2:end-1)) ./ dxp(2:end-1,2:end)...
            + Kwxm_dt_ibj(:,1:end-1) ...
            .* (pw(2:end-1,2:end-1) - pw(2:end-1,1:end-2)) ./ dxp(2:end-1,1:end-1));
        dRwdt_ibj = (dRwdt_ibj_p3 - Rwp3_i) / dthew;

        % dRwdthew i+1,j
        dRwdt_iaj_p3 = 1 ./ xgs .* (-Kwxm_dt_iaj(:,2:end) ...
            .* (pw(2:end-1,3:end) - pw(2:end-1,2:end-1)) ./ dxp(2:end-1,2:end)...
            + Kwxm(:,1:end-1) ...
            .* (pw(2:end-1,2:end-1) - pw(2:end-1,1:end-2)) ./ dxp(2:end-1,1:end-1));
        dRwdt_iaj = (dRwdt_iaj_p3 - Rwp3_i) / dthew;

        % dRwdthew i,j-1
        dRwdt_ijb_p3 = 1 ./ zgs .* (-Kwzm(2:end,:) ...
            .* ((pw(3:end,2:end-1) - pw(2:end-1,2:end-1)) ./ dzp(2:end,2:end-1)...
            - (rhow(3:end,2:end-1) + rhow(2:end-1,2:end-1)) / 2 * g * gr) ...
            + Kwzm_dt_ijb(1:end-1,:) ...
            .* ((pw(2:end-1,2:end-1) - pw(1:end-2,2:end-1)) ./ dzp(1:end-1,2:end-1) ...
            - (rhow(2:end-1,2:end-1) + rhow(1:end-2,2:end-1)) / 2 * g * gr));
        dRwdt_ijb = (dRwdt_ijb_p3 - Rwp3_j) / dthew;

        % dRwdthew i,j+1
        dRwdt_ija_p3 = 1 ./ zgs .* (-Kwzm_dt_ija(2:end,:) ...
            .* ((pw(3:end,2:end-1) - pw(2:end-1,2:end-1)) ./ dzp(2:end,2:end-1)...
            - (rhow(3:end,2:end-1) + rhow(2:end-1,2:end-1)) / 2 * g * gr) ...
            + Kwzm(1:end-1,:) ...
            .* ((pw(2:end-1,2:end-1) - pw(1:end-2,2:end-1)) ./ dzp(1:end-1,2:end-1) ...
            - (rhow(2:end-1,2:end-1) + rhow(1:end-2,2:end-1)) / 2 * g * gr));
        dRwdt_ija = (dRwdt_ija_p3 - Rwp3_j) / dthew;

        % dRndthew i,j
        dRndt_ij_p1 = (rhon_dt(2:end-1,2:end-1) .* (phi(2:end-1,2:end-1) - thew_dt(2:end-1,2:end-1)) ...
            - rhonn(2:end-1,2:end-1) .* (phi(2:end-1,2:end-1) - thewn(2:end-1,2:end-1))) ...
            / dt;
        dRndt_ij_p3 = 1 ./ xgs .* (-Knxm_dt_ibj(:,2:end) ...
            .* (pw(2:end-1,3:end) - pw(2:end-1,2:end-1) + pc(2:end-1,3:end) - pc_dt(2:end-1,2:end-1))...
            ./ dxp(2:end-1,2:end)...
            + Knxm_dt_iaj(:,1:end-1) ...
            .* (pw(2:end-1,2:end-1) - pw(2:end-1,1:end-2) + pc_dt(2:end-1,2:end-1) - pc(2:end-1,1:end-2)) ...
            ./ dxp(2:end-1,1:end-1)) ...
            + 1 ./ zgs .* (-Knzm_dt_ijb(2:end,:) ...
            .* ((pw(3:end,2:end-1) - pw(2:end-1,2:end-1) + pc(3:end,2:end-1) - pc_dt(2:end-1,2:end-1))...
            ./ dzp(2:end,2:end-1)...
            - (rhon(3:end,2:end-1) + rhon_dt(2:end-1,2:end-1)) / 2 * g * gr) ...
            + Knzm_dt_ija(1:end-1,:) ...
            .* ((pw(2:end-1,2:end-1) - pw(1:end-2,2:end-1) + pc_dt(2:end-1,2:end-1) - pc(1:end-2,2:end-1)) ...
            ./ dzp(1:end-1,2:end-1) ...
            - (rhon_dt(2:end-1,2:end-1) + rhon(1:end-2,2:end-1)) / 2 * g * gr));
        dRndt_ij = ((dRndt_ij_p1 + dRndt_ij_p3) - Rn) / dthew;

        % dRndthew i-1,j
        dRndt_ibj_p3 = 1 ./ xgs .* (-Knxm(:,2:end) ...
            .* (pw(2:end-1,3:end) - pw(2:end-1,2:end-1) + pc(2:end-1,3:end) - pc(2:end-1,2:end-1))...
            ./ dxp(2:end-1,2:end)...
            + Knxm_dt_ibj(:,1:end-1) ...
            .* (pw(2:end-1,2:end-1) - pw(2:end-1,1:end-2) + pc(2:end-1,2:end-1) - pc_dt(2:end-1,1:end-2)) ...
            ./ dxp(2:end-1,1:end-1));
        dRndt_ibj = (dRndt_ibj_p3 - Rnp3_i) / dthew;

        % dRndthew i+1,j
        dRndt_iaj_p3 = 1 ./ xgs .* (-Knxm_dt_iaj(:,2:end) ...
            .* (pw(2:end-1,3:end) - pw(2:end-1,2:end-1) + pc_dt(2:end-1,3:end) - pc(2:end-1,2:end-1))...
            ./ dxp(2:end-1,2:end)...
            + Knxm(:,1:end-1) ...
            .* (pw(2:end-1,2:end-1) - pw(2:end-1,1:end-2) + pc(2:end-1,2:end-1) - pc(2:end-1,1:end-2)) ...
            ./ dxp(2:end-1,1:end-1));
        dRndt_iaj = (dRndt_iaj_p3 - Rnp3_i) / dthew;

        % dRndthew i,j-1
        dRndt_ijb_p3 = 1 ./ zgs .* (-Knzm(2:end,:) ...
            .* ((pw(3:end,2:end-1) - pw(2:end-1,2:end-1) + pc(3:end,2:end-1) - pc(2:end-1,2:end-1))...
            ./ dzp(2:end,2:end-1)...
            - (rhon(3:end,2:end-1) + rhon(2:end-1,2:end-1)) / 2 * g * gr) ...
            + Knzm_dt_ijb(1:end-1,:) ...
            .* ((pw(2:end-1,2:end-1) - pw(1:end-2,2:end-1) + pc(2:end-1,2:end-1) - pc_dt(1:end-2,2:end-1)) ...
            ./ dzp(1:end-1,2:end-1) ...
            - (rhon(2:end-1,2:end-1) + rhon_dt(1:end-2,2:end-1)) / 2 * g * gr));
        dRndt_ijb = (dRndt_ijb_p3 - Rnp3_j) / dthew;

        % dRndthew i,j+1
        dRndt_ija_p3 = 1 ./ zgs .* (-Knzm_dt_ija(2:end,:) ...
            .* ((pw(3:end,2:end-1) - pw(2:end-1,2:end-1) + pc_dt(3:end,2:end-1) - pc(2:end-1,2:end-1))...
            ./ dzp(2:end,2:end-1)...
            - (rhon_dt(3:end,2:end-1) + rhon(2:end-1,2:end-1)) / 2 * g * gr) ...
            + Knzm(1:end-1,:) ...
            .* ((pw(2:end-1,2:end-1) - pw(1:end-2,2:end-1) + pc(2:end-1,2:end-1) - pc(1:end-2,2:end-1)) ...
            ./ dzp(1:end-1,2:end-1) ...
            - (rhon(2:end-1,2:end-1) + rhon(1:end-2,2:end-1)) / 2 * g * gr));
        dRndt_ija = (dRndt_ija_p3 - Rnp3_j) / dthew;

        %% Matrix A dpw
        % dRwdpw i,j
        dRwdpw_ij_p1 = (rhow_dp(2:end-1,2:end-1) .* thew(2:end-1,2:end-1) - ...
            rhown(2:end-1,2:end-1) .* thewn(2:end-1,2:end-1)) / dt;
        dRwdpw_ij_p2 = Ss ./ phi(2:end-1,2:end-1) .* thew(2:end-1,2:end-1) ...
            .* rhow_dp(2:end-1,2:end-1) ...
            .* (pw_dp(2:end-1,2:end-1) - pwn(2:end-1,2:end-1)) ...
            / (dt * g * rhow0);
        dRwdpw_ij_p3 = 1 ./ xgs .* (-Kwxm_dp_ibj(:,2:end) ...
            .* (pw(2:end-1,3:end) - pw_dp(2:end-1,2:end-1)) ./ dxp(2:end-1,2:end)...
            + Kwxm_dp_iaj(:,1:end-1) ...
            .* (pw_dp(2:end-1,2:end-1) - pw(2:end-1,1:end-2)) ./ dxp(2:end-1,1:end-1))...
            + 1 ./ zgs .* (-Kwzm_dp_ijb(2:end,:) ...
            .* ((pw(3:end,2:end-1) - pw_dp(2:end-1,2:end-1)) ./ dzp(2:end,2:end-1)...
            - (rhow(3:end,2:end-1) + rhow_dp(2:end-1,2:end-1)) / 2 * g * gr) ...
            + Kwzm_dp_ija(1:end-1,:) ...
            .* ((pw_dp(2:end-1,2:end-1) - pw(1:end-2,2:end-1)) ./ dzp(1:end-1,2:end-1) ...
            - (rhow_dp(2:end-1,2:end-1) + rhow(1:end-2,2:end-1)) / 2 * g * gr));
        dRwdpw_ij = (dRwdpw_ij_p1 + dRwdpw_ij_p2 + dRwdpw_ij_p3 - Rw) / dpw;

        % dRwdpw i-1,j
        dRwdpw_ibj_p3 = 1 ./ xgs .* (-Kwxm(:,2:end) ...
            .* (pw(2:end-1,3:end) - pw(2:end-1,2:end-1)) ./ dxp(2:end-1,2:end)...
            + Kwxm_dp_ibj(:,1:end-1) ...
            .* (pw(2:end-1,2:end-1) - pw_dp(2:end-1,1:end-2)) ./ dxp(2:end-1,1:end-1));
        dRwdpw_ibj = (dRwdpw_ibj_p3 - Rwp3_i) / dpw;

        % dRwdpw i+1,j
        dRwdpw_iaj_p3 = 1 ./ xgs .* (-Kwxm_dp_iaj(:,2:end) ...
            .* (pw_dp(2:end-1,3:end) - pw(2:end-1,2:end-1)) ./ dxp(2:end-1,2:end)...
            + Kwxm(:,1:end-1) ...
            .* (pw(2:end-1,2:end-1) - pw(2:end-1,1:end-2)) ./ dxp(2:end-1,1:end-1));
        dRwdpw_iaj = (dRwdpw_iaj_p3 - Rwp3_i) / dpw;

        % dRwdpw i,j-1
        dRwdpw_ijb_p3 = 1 ./ zgs .* (-Kwzm(2:end,:) ...
            .* ((pw(3:end,2:end-1) - pw(2:end-1,2:end-1)) ./ dzp(2:end,2:end-1)...
            - (rhow(3:end,2:end-1) + rhow(2:end-1,2:end-1)) / 2 * g * gr) ...
            + Kwzm_dp_ijb(1:end-1,:) ...
            .* ((pw(2:end-1,2:end-1) - pw_dp(1:end-2,2:end-1)) ./ dzp(1:end-1,2:end-1) ...
            - (rhow(2:end-1,2:end-1) + rhow_dp(1:end-2,2:end-1)) / 2 * g * gr));
        dRwdpw_ijb = (dRwdpw_ijb_p3 - Rwp3_j) / dpw;

        % dRwdpw i,j+1
        dRwdpw_ija_p3 = 1 ./ zgs .* (-Kwzm_dp_ija(2:end,:) ...
            .* ((pw_dp(3:end,2:end-1) - pw(2:end-1,2:end-1)) ./ dzp(2:end,2:end-1)...
            - (rhow_dp(3:end,2:end-1) + rhow(2:end-1,2:end-1)) / 2 * g * gr) ...
            + Kwzm(1:end-1,:) ...
            .* ((pw(2:end-1,2:end-1) - pw(1:end-2,2:end-1)) ./ dzp(1:end-1,2:end-1) ...
            - (rhow(2:end-1,2:end-1) + rhow(1:end-2,2:end-1)) / 2 * g * gr));
        dRwdpw_ija = (dRwdpw_ija_p3 - Rwp3_j) / dpw;

        % dRndpw i,j
        dRndpw_ij_p1 = (rhon_dp(2:end-1,2:end-1) .* (phi(2:end-1,2:end-1) - thew(2:end-1,2:end-1)) ...
            - rhonn(2:end-1,2:end-1) .* (phi(2:end-1,2:end-1) - thewn(2:end-1,2:end-1))) ...
            / dt;
        dRndpw_ij_p3 = 1 ./ xgs .* (-Knxm_dp_ibj(:,2:end) ...
            .* (pw(2:end-1,3:end) - pw_dp(2:end-1,2:end-1) + pc(2:end-1,3:end) - pc(2:end-1,2:end-1))...
            ./ dxp(2:end-1,2:end)...
            + Knxm_dp_iaj(:,1:end-1) ...
            .* (pw_dp(2:end-1,2:end-1) - pw(2:end-1,1:end-2) + pc(2:end-1,2:end-1) - pc(2:end-1,1:end-2)) ...
            ./ dxp(2:end-1,1:end-1))...
            + 1 ./ zgs .* (-Knzm_dp_ijb(2:end,:) ...
            .* ((pw(3:end,2:end-1) - pw_dp(2:end-1,2:end-1) + pc(3:end,2:end-1) - pc(2:end-1,2:end-1))...
            ./ dzp(2:end,2:end-1)...
            - (rhon(3:end,2:end-1) + rhon_dp(2:end-1,2:end-1)) / 2 * g * gr) ...
            + Knzm_dp_ija(1:end-1,:) ...
            .* ((pw_dp(2:end-1,2:end-1) - pw(1:end-2,2:end-1) + pc(2:end-1,2:end-1) - pc(1:end-2,2:end-1)) ...
            ./ dzp(1:end-1,2:end-1) ...
            - (rhon_dp(2:end-1,2:end-1) + rhon(1:end-2,2:end-1)) / 2 * g * gr));
        dRndpw_ij = (dRndpw_ij_p1 + dRndpw_ij_p3 - Rn) / dpw;

        % dRndpw i-1,j
        dRndpw_ibj_p3 = 1 ./ xgs .* (-Knxm(:,2:end) ...
            .* (pw(2:end-1,3:end) - pw(2:end-1,2:end-1) + pc(2:end-1,3:end) - pc(2:end-1,2:end-1))...
            ./ dxp(2:end-1,2:end)...
            + Knxm_dp_ibj(:,1:end-1) ...
            .* (pw(2:end-1,2:end-1) - pw_dp(2:end-1,1:end-2) + pc(2:end-1,2:end-1) - pc(2:end-1,1:end-2)) ...
            ./ dxp(2:end-1,1:end-1));
        dRndpw_ibj = (dRndpw_ibj_p3 - Rnp3_i) / dpw;

        % dRndpw i+1,j
        dRndpw_iaj_p3 = 1 ./ xgs .* (-Knxm_dp_iaj(:,2:end) ...
            .* (pw_dp(2:end-1,3:end) - pw(2:end-1,2:end-1) + pc(2:end-1,3:end) - pc(2:end-1,2:end-1))...
            ./ dxp(2:end-1,2:end)...
            + Knxm(:,1:end-1) ...
            .* (pw(2:end-1,2:end-1) - pw(2:end-1,1:end-2) + pc(2:end-1,2:end-1) - pc(2:end-1,1:end-2)) ...
            ./ dxp(2:end-1,1:end-1));
        dRndpw_iaj = (dRndpw_iaj_p3 - Rnp3_i) / dpw;

        % dRndpw i,j-1
        dRndpw_ijb_p3 = 1 ./ zgs .* (-Knzm(2:end,:) ...
            .* ((pw(3:end,2:end-1) - pw(2:end-1,2:end-1) + pc(3:end,2:end-1) - pc(2:end-1,2:end-1))...
            ./ dzp(2:end,2:end-1)...
            - (rhon(3:end,2:end-1) + rhon(2:end-1,2:end-1)) / 2 * g * gr) ...
            + Knzm_dp_ijb(1:end-1,:) ...
            .* ((pw(2:end-1,2:end-1) - pw_dp(1:end-2,2:end-1) + pc(2:end-1,2:end-1) - pc(1:end-2,2:end-1)) ...
            ./ dzp(1:end-1,2:end-1) ...
            - (rhon(2:end-1,2:end-1) + rhon_dp(1:end-2,2:end-1)) / 2 * g * gr));
        dRndpw_ijb = (dRndpw_ijb_p3 - Rnp3_j) / dpw;

        % dRndpw i,j+1
        dRndpw_ija_p3 = 1 ./ zgs .* (-Knzm_dp_ija(2:end,:) ...
            .* ((pw_dp(3:end,2:end-1) - pw(2:end-1,2:end-1) + pc(3:end,2:end-1) - pc(2:end-1,2:end-1))...
            ./ dzp(2:end,2:end-1)...
            - (rhon_dp(3:end,2:end-1) + rhon(2:end-1,2:end-1)) / 2 * g * gr) ...
            + Knzm(1:end-1,:) ...
            .* ((pw(2:end-1,2:end-1) - pw(1:end-2,2:end-1) + pc(2:end-1,2:end-1) - pc(1:end-2,2:end-1)) ...
            ./ dzp(1:end-1,2:end-1) ...
            - (rhon(2:end-1,2:end-1) + rhon(1:end-2,2:end-1)) / 2 * g * gr));
        dRndpw_ija = (dRndpw_ija_p3 - Rnp3_j) / dpw;

        %% Assemble matrix A and vector b
        % Source/term to flowrate condition (for upper boundary)
        if SStoFD == 1
            if UP_stage1_option == 1
                Rw(1,SSq_left:SSq_right) = Rw(1,SSq_left:SSq_right) + SSq;
            end
        end
        % Vector b
        b = [Rw(:);Rn(:)];

        % Process A_Rwthew
        % Last row remove
        dRwdt_ija(end,:) = 0;
        dRwdt_ija = dRwdt_ija(:);
        dRwdt_ija = dRwdt_ija(1:end-1);
        % First row remove
        dRwdt_ijb(1,:) = 0;
        dRwdt_ijb = dRwdt_ijb(:);
        dRwdt_ijb = dRwdt_ijb(2:end);
        % Last column remove
        dRwdt_iaj = dRwdt_iaj(:,1:end-1);
        dRwdt_iaj = dRwdt_iaj(:);
        % First column remove
        dRwdt_ibj = dRwdt_ibj(:,2:end);
        dRwdt_ibj = dRwdt_ibj(:);
        % Assemble matrix A_Rwthew
        A_Rwthew = diag(dRwdt_ij(:))...
            + diag(dRwdt_ija,1)...
            + diag(dRwdt_ijb,-1)...
            + diag(dRwdt_iaj,nz)...
            + diag(dRwdt_ibj,-nz);

        % Process A_Rwpw
        % Last row remove
        dRwdpw_ija(end,:) = 0;
        dRwdpw_ija = dRwdpw_ija(:);
        dRwdpw_ija = dRwdpw_ija(1:end-1);
        % First row remove
        dRwdpw_ijb(1,:) = 0;
        dRwdpw_ijb = dRwdpw_ijb(:);
        dRwdpw_ijb = dRwdpw_ijb(2:end);
        % Last column remove
        dRwdpw_iaj = dRwdpw_iaj(:,1:end-1);
        dRwdpw_iaj = dRwdpw_iaj(:);
        % First column remove
        dRwdpw_ibj = dRwdpw_ibj(:,2:end);
        dRwdpw_ibj = dRwdpw_ibj(:);
        % Assemble matrix A_Rwpw
        A_Rwpw = diag(dRwdpw_ij(:))...
            + diag(dRwdpw_ija,1)...
            + diag(dRwdpw_ijb,-1)...
            + diag(dRwdpw_iaj,nz)...
            + diag(dRwdpw_ibj,-nz);

        % Process A_Rnthew
        % Last row remove
        dRndt_ija(end,:) = 0;
        dRndt_ija = dRndt_ija(:);
        dRndt_ija = dRndt_ija(1:end-1);
        % First row remove
        dRndt_ijb(1,:) = 0;
        dRndt_ijb = dRndt_ijb(:);
        dRndt_ijb = dRndt_ijb(2:end);
        % Last column remove
        dRndt_iaj = dRndt_iaj(:,1:end-1);
        dRndt_iaj = dRndt_iaj(:);
        % First column remove
        dRndt_ibj = dRndt_ibj(:,2:end);
        dRndt_ibj = dRndt_ibj(:);
        % Assemble matrix A_Rnthew
        A_Rnthew = diag(dRndt_ij(:))...
            + diag(dRndt_ija,1)...
            + diag(dRndt_ijb,-1)...
            + diag(dRndt_iaj,nz)...
            + diag(dRndt_ibj,-nz);

        % Process A_Rnpw
        % Last row remove
        dRndpw_ija(end,:) = 0;
        dRndpw_ija = dRndpw_ija(:);
        dRndpw_ija = dRndpw_ija(1:end-1);
        % First row remove
        dRndpw_ijb(1,:) = 0;
        dRndpw_ijb = dRndpw_ijb(:);
        dRndpw_ijb = dRndpw_ijb(2:end);
        % Last column remove
        dRndpw_iaj = dRndpw_iaj(:,1:end-1);
        dRndpw_iaj = dRndpw_iaj(:);
        % First column remove
        dRndpw_ibj = dRndpw_ibj(:,2:end);
        dRndpw_ibj = dRndpw_ibj(:);
        % Assemble matrix A_Rnpw
        A_Rnpw = diag(dRndpw_ij(:))...
            + diag(dRndpw_ija,1)...
            + diag(dRndpw_ijb,-1)...
            + diag(dRndpw_iaj,nz)...
            + diag(dRndpw_ibj,-nz);

        % Assemble matrix A
        A = [A_Rwthew A_Rwpw; A_Rnthew A_Rnpw];
        A = sparse(A);

        % NR solution and iteration
        % Del solution
        thew_pw_del = A \ -b;
        thew_del = thew_pw_del(1:nx*nz);
        thew_cri = max(abs(thew_del));
        thew_del = reshape(thew_del,nz,nx);
        pw_del = thew_pw_del(nx*nz+1:end);
        pw_cri = max(abs(pw_del));
        pw_del = reshape(pw_del,nz,nx);

        % Add
        thew(2:end-1,2:end-1) = thew(2:end-1,2:end-1) + thew_del;
        pw(2:end-1,2:end-1) = pw(2:end-1,2:end-1) + pw_del;

        % Criterion
        if thew_cri <= ethew && pw_cri <= epw 
            break
        end
        % Adjust timestep for convergence
        if inu >= inu_limit
            thew = thewn;
            pw = pwn;
            T = T - dt;
            dt = 0.8 * dt;
            inu_back = 1;
            break
        end
    end

    if inu_back == 0
        % Mass balance
        % Parameters
        seff_mas = (thew - thewr) ./ (phi - thewr - thenr);
        seff_mas(seff_mas > 1) = 1;

        pc_mas = TP_VGkrpc(MVG_option,seff_mas,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);
        [~,krw_mas] = TP_VGkrpc(MVG_option,seff_mas,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);
        [~,~,krn_mas] = TP_VGkrpc(MVG_option,seff_mas,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);

        rhow_mas = rhow0 * exp(aw * (pw - pw0));
        rhon_mas = rhon0 + bn * (pw + pc_mas - pn0);

        Kw_mas = k .* krw_mas .* rhow_mas / vw;
        Kn_mas = k .* krn_mas .* rhon_mas / vn;
        Kwzm_mas = (Kw_mas(1:end-1,:) + Kw_mas(2:end,:)) / 2;
        Knzm_mas = (Kn_mas(1:end-1,:) + Kn_mas(2:end,:)) / 2;

        Knzm_topclose = find(thew(1,:) >= thews(1,:));
        Knzm_mas(1,Knzm_topclose) = 0;

        if bottom_close == 1
            Kwzm_mas(end,:) = 0;
            Knzm_mas(end,:) = 0;
        end
        if id.top_wclose_nopen == 1
            Kwzm(1,:) = 0;
            Kwzm_dt_ija(1,:) = 0;
            Kwzm_dt_ijb(1,:) = 0;
            Kwzm_dp_ija(1,:) = 0;
            Kwzm_dp_ijb(1,:) = 0;
        end

        % Calculate
        qw = -Kwzm_mas .* ((pw(2:end,:) - pw(1:end-1,:)) ./ dzp(:,:) - (rhow_mas(1:end-1,:) + rhow_mas(2:end,:)) / 2 * g * gr);
        dqw = dqw + sum(qw(1,id.SSq_left+1:id.SSq_right+1) * dt) - ...
            sum(qw(end,2:end-1) * dt);
        dmqw = (thew(2:end-1,2:end-1)) .* zgs .* rhow_mas(2:end-1,2:end-1) ...
            - (id.thew(2:end-1,2:end-1)) .* zgs .* rhow_mas(2:end-1,2:end-1);
        masdqw = (sum(sum(dmqw)) - dqw);

        % 计算应有的积水量值(单位rhow*m/s)
        if SStoFD == 1
            qwdet = qwdet + (id.q + qw(1)) * dt;
        end

        qn = -Knzm_mas .* ((pw(2:end,:) - pw(1:end-1,:)) ./ dzp(:,:) ...
            + (pc_mas(2:end,:) - pc_mas(1:end-1,:)) ./ dzp(:,:)...
            - (rhon_mas(1:end-1,:) + rhon_mas(2:end,:)) / 2 * g * gr);
        dqn = dqn + sum(qn(1,id.SSq_left+1:id.SSq_right+1) * dt) - ...
            sum(qn(end,2:end-1) * dt);
        dmqn = (phi(2:end-1,2:end-1) - thew(2:end-1,2:end-1)) .* zgs .* rhon_mas(2:end-1,2:end-1) ...
            - (phi(2:end-1,2:end-1) - id.thew(2:end-1,2:end-1)) .* zgs .* id.rhon0;
        masdqn = (sum(sum(dmqn)) - dqn);
        un = qn ./ ((rhon_mas(1:end-1,:) + rhon_mas(2:end,:)) / 2);

        % Recond output
        if T >= T_rec(recond_index)
            if T < T_rec(end)
                if T>= T_rec(recond_index + 1)
                    outtime = [outtime recond_index];
                end
            end
            thew_rec(:,:,recond_index) = thew ...
                - ((T - T_rec(recond_index)) / dt) .* (thew - thewn);
            pw_rec(:,:,recond_index) = pw...
                - ((T - T_rec(recond_index)) / dt) .* (pw - pwn);
            seff_pc = (thew - thewr) ...
                ./ (phi - thewr - thenr);
            seff_pcn = (thewn - thewr) ...
                ./ (phi - thewr - thenr);
            seff_pc(seff_pc >= 1) = 1;
            seff_pcn(seff_pcn >= 1) = 1;
            pc_pc = TP_VGkrpc(MVG_option,seff_pc,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);
            pc_pcn = TP_VGkrpc(MVG_option,seff_pcn,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g);
            pc_rec(:,:,recond_index) = pc_pc ...
                - ((T - T_rec(recond_index)) / dt) .* (pc_pc - pc_pcn);
            pn_rec(:,:,recond_index) = pw_rec(:,:,recond_index) ...
                + pc_rec(:,:,recond_index);
            qw_rec(:,:,recond_index) = qw ...
                - ((T - T_rec(recond_index)) / dt) .* (qw - qwn);
            qn_rec(:,:,recond_index) = qn ...
                - ((T - T_rec(recond_index)) / dt) .* (qn - qnn);
            un_rec(:,:,recond_index) = un ...
                - ((T - T_rec(recond_index)) / dt) .* (un - unn);
            recond_index = recond_index + 1;
            if T >= T_rec(end)
                break
            end
        end

    end
    % Update
    unn = un;
    qwn = qw;
    qnn = qn;
    thewn = thew;
    pwn = pw;

    % Recond iteration times
    inum = [inum; inu];

    % Timestep adjust
    if inu <= inu_ac
        dt = dt * 1.2;
    end
    if dt >= dt_max
        dt = dt_max;
    end

    % not convergence
    if dt < id.dt * 0.1
        tp_outline = 1;
        break
    end
%     if sum(inum == 10) > 20
%         T = 0;
%         break
%     end   

end

% Simulation end ponding pw (if exist)
if SStoFD == 1
    if UP_stage2_option == 1
        pw_up2_ini = pw(1,curan);
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

    end
    if UP_stage3_option == 1
        pw_up3_ini = pw(1,curan);
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

t_operTP = toc;

% Save to ''mydata.mat''
ds.thewr = thewr;
ds.thenr = thenr;
ds.phi = phi;
ds.tp_outtime = outtime;
ds.tp_nt = inum;
ds.tp_oper = t_operTP;
ds.tp_thew = thew_rec;
ds.tp_pw = pw_rec;
ds.tp_pn = pn_rec;
ds.tp_pc = pc_rec;
ds.tp_uw = qw_rec / rhow0;
ds.tp_un = un_rec;
ds.tp_wmass = abs(masdqw / dqw);
ds.tp_nmass = abs(masdqn / dqn);
ds.T_tp_pond = T_tp_pond;
ds.tp_outline = tp_outline;
save('mydata.mat', 'ds');