% 1D infiltration problem from Celia et al.(1992)
% 采用特殊的相对渗透率-有效饱和度关系，见'Kwtvcb.m'与'Kntvcb.m'
clear
clc

% Domain (zooming only)
% Location of points
x = [1];
z = [-0.005:-0.01:-1.065]';
% Number of grids (row and column)
nx = length(x);
nz = length(z);
% Grid size
xgs = 1 * ones(length(z),length(x));
zgs = 0.01 * ones(length(z),length(x));
% Distance between points
dxp = 1 * ones(length(z) + 2,length(x) + 1);
dzp = 0.01 * ones(length(z) + 1,length(x) + 2);
% Time step size
dt = 0.1;
% Tiny step for der
dthew = -1 * 10 ^ -5;
dpw = -1 * 10 ^ -3;
% Compressible parameters
aw = 0;
bn = 1.189e-5;
% Reference density
rhow0 = 1000;
rhon0 = 1.29;
% Viscosity
vw = 1 * 10^-3;
vn = 1.83 * 10^-5;
% Soil parameters
n = zeros(nz+2,nx+2);
a = zeros(nz+2,nx+2);
% Specific storage
Ss = 1 * 10 ^ -5;
% Porosity
phi = zeros(nz+2,nx+2);
% Gravity
g = 9.81;
% Gracity option
gr = 1;
% Absolute permeability
k = zeros(nz+2,nx+2);
% Residual & Saturated
thewr = zeros(nz+2,nx+2);
thews = zeros(nz+2,nx+2);
thenr = zeros(nz+2,nx+2);
thens = zeros(nz+2,nx+2);
% Set homosoil
n(:,:) = 2.2;
a(:,:) = 4.4;
m = 1 - 1 ./ n;
pcm = rhow0 * g ./ a;
Aw = 18130 * 0.01 / 3600;
Bw = 6.07;
Aa = 3.6e-5;
Ba = -2.4;
Kas = 2800 * 0.01 / 3600;
phi(:,:) = 0.37;
thewr(:,:) = 0.0265;
thenr(:,:) = 0.058;
thews(:,:) = phi - thenr;
thens(:,:) = phi - thewr;
% Precision
ethew = 1 * 10 ^ -3;
epw = 1;
% Timestep max value & adjust criterion
dt_max = 10;
inu_ac = 10;
% Total simulation time
T = 0;
% Recond time
T_rec = 60 * [ 3 6 12 18 24];
recond_index = 1;
% Create output variables
thew = zeros(nz+2,nx+2);
pw = zeros(nz+2,nx+2);
% Intermediate valuables
pc = zeros(nz+2,nx+2);
pn = zeros(nz+2,nx+2);
seff = zeros(nz+2,nx+2);
b = zeros(nz*nx,1);
A = zeros(nz*nx,nz*nx);
% Initial and boundary conditions
% Pressure head boundary
hw = [-1.075:0.01:0.005]';
pw(:,1) = hw * rhow0 * g;
pw(:,2) = hw * rhow0 * g;
pw(:,3) = hw * rhow0 * g;

% Van Genuchten transformation
thew = (1 + abs(a .* pw ./ (rhow0 * g)).^n).^-m .* ...
    (phi - thewr - thenr) + thewr;
% Reference
thew0 = thew;
pw0 = pw;
pn0 = pn;
% Iteration valuables
pwn = pw;
thewn = thew;
pwnn = pw;
thewnn = thew;
inum = [];
% Output recond valuables
thew_rec = zeros(nz+2,nx+2,length(T_rec));
pw_rec = zeros(nz+2,nx+2,length(T_rec));
pc_rec = zeros(nz+2,nx+2,length(T_rec));
pn_rec = zeros(nz+2,nx+2,length(T_rec));

% Boundary condition options
left_close = 1;
right_close = 1;
top_close = 0;
bottom_close = 1;

% Sink/Source term to flowrate boundary options
SStoFD = 1;
if SStoFD == 1
    q = -20 * 0.01 / 3600 * rhow0;
    SSq = q / dzp(1,2);
    SSq_left = 1;
    SSq_right = 1;
end

% Change upper boundary to input boundary value
UP_stage1_option = 1;
UP_stage2_option = 0;
UP_stage3_option = 0;
pwbmi = pw(1,2);
pwbmm = -pw(1,2);
pnae = rhow0 * g * 14 * 0.01;
thewae = (1 + abs(a .* pnae ./ rhow0 / g).^n).^-m .* ...
    (phi - thewr - thenr) + thewr;
epw_BM = 1e-1;
curan = 2;

% Operation time
tic;

% Simulation
while 1

    % Time increasing
    T = T + dt;
    % Iteration times
    inu = 0;

    if SStoFD == 1
        if UP_stage2_option == 1
            pw_up2_ini = pw(1,curan);
        end
        if UP_stage3_option == 1
            pw_up3_ini = pw(1,curan);
        end
    end

    while 1
        
        % Get upper boundary water content 
        if UP_stage1_option == 1
            [thewtop, pwtop] = getupperboundaryBM(rhow0,pwbmi,pwbmm,pw(2,curan)...
                ,thew(2,curan),q,dzp(1,curan),n(1,curan),m(1,curan),k(1,curan),vw,a(1,curan),...
                g,gr,thewr(1,curan),thenr(1,curan),phi(1,curan),epw_BM,Aw,Bw);
            thew(1,curan) = thewtop;
            pw(1,curan) = pwtop;
            if pwtop > 0
                UP_stage1_option = 0;
                UP_stage2_option = 1;
                SStoFD = 0;
                T_ponding = T;
                pw(1,curan) = 0;
                pw_up2_ini = pw(1,curan);
            end
        end
        if UP_stage2_option == 1
            thew(1,curan) = phi(1,curan) - thenr(1,curan);
            K1_up2 = (Kwtvcb(thew(1,curan),Aw,Bw,g) + Kwtvcb(thew(2,curan),Aw,Bw,g)) / rhow0 / 2;
            K_up2 = K1_up2;
            pw_dtup2 = (-q / rhow0 ./ K_up2 .* dzp(1,curan) - pw(1,curan) ...
                + pw(2,curan)  - rhow0 * g * dzp(1,curan)) ...
                / (1 + dzp(1,curan) ./ (rhow0 * g * dt * K_up2));
            pw(1,curan) = pw_up2_ini + pw_dtup2;
            pc2_up2 = pcm(2,curan) .* (seff(2,curan) ...
                .^ (-1 ./ m(2,curan)) -1) .^ (1 ./ n(2,curan));
            pn2_up2 = pw(2,curan) + pc2_up2;
            if pw(1,curan) < 0
                UP_stage2_option = 0;
                UP_stage1_option = 1;
            end
            if pn2_up2 >= (pnae + pw(1,curan)) && pw(1,curan) > 0
                UP_stage2_option = 0;
                UP_stage3_option = 1;
                T_crit = T;
                pw_up3_ini = pw(1,curan);
            end
        end
        if UP_stage3_option == 1
            thew(1,curan) = thewae(1,curan);
            K1_up3 = Kwtvcb(thew(1,curan),Aw,Bw,g) / rhow0;
            K2_up3 = Kwtvcb(thew(2,curan),Aw,Bw,g) / rhow0;
            K_up3 = (K1_up3 + K2_up3) / 2;
            pw_dtup3 = (-q / rhow0 ./ K_up3 .* dzp(1,curan) - pw(1,curan) ...
                + pw(2,curan)  - rhow0 * g * dzp(1,curan)) ...
                / (1 + dzp(1,curan) ./ (rhow0 * g * dt * K_up3));
            pw(1,curan) = pw_up3_ini + pw_dtup3;
        end  

        % Iteration increasing
        inu = inu + 1;

        % Effective saturation & relative permeability
        seff = (thew - thewr) ./ (phi - thewr - thenr);
        seffn = (thewn - thewr) ./ (phi - thewr - thenr);
        seff(seff > 1) = 1;
        seffn(seffn > 1) = 1;

        % Capillary pressure
        pc = pcm .* (seff .^ (-1 ./ m) -1) .^ (1 ./ n);
        pcn = pcm .* (seffn .^ (-1 ./ m) -1) .^ (1 ./ n);

        % Density
        rhow = rhow0 * exp(aw * (pw - pw0));
        rhon = rhon0 + bn * (pw + pc - pn0);
        rhown = rhow0 * exp(aw * (pwn - pw0));
        rhonn = rhon0 + bn * (pwn + pcn - pn0);

        % Permeability (intercell)
        Kw = Kwtvcb(thew,Aw,Bw,g);
        Kn = Kntvcb(Kas,Aa,Ba,pc/(rhow0*g),g,rhon0,rhow0);
        % Arithmetic
        Kwxm = (Kw(2:end-1,1:end-1) + Kw(2:end-1,2:end)) / 2;
        Knxm = (Kn(2:end-1,1:end-1) + Kn(2:end-1,2:end)) / 2;
        Kwzm = (Kw(1:end-1,2:end-1) + Kw(2:end,2:end-1)) / 2;
        Knzm = (Kn(1:end-1,2:end-1) + Kn(2:end,2:end-1)) / 2;

        % Prepare valuables for matrix A dthew
        thew_dt = thew + dthew;
        seff_dt = (thew_dt - thewr) ./ (phi - thewr - thenr);
        pc_dt = pcm .* (seff_dt .^ (-1 ./ m) -1) .^ (1 ./ n);
        rhon_dt = rhon0 + bn * (pw + pc_dt - pn0);

        Kw_dt = Kwtvcb(thew_dt,Aw,Bw,g);
        Kn_dt = Kntvcb(Kas,Aa,Ba,pc_dt/(rhow0*g),g,rhon0,rhow0);
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
        Kw_dp = Kwtvcb(thew,Aw,Bw,g);
        Kn_dp = Kntvcb(Kas,Aa,Ba,pc/(rhow0*g),g,rhon0,rhow0);
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
        
        % Source/term to flowrate condition (for upper boundary)
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
        if UP_stage1_option == 1
            Rw(1,SSq_left:SSq_right) = Rw(1,SSq_left:SSq_right) + SSq;
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

        %% NR solution and iteration
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
        if inu >= 15
            thew = thewn;
            pw = pwn;
            T = T - dt;
            dt = 0.8 * dt;
            break
        end
    end

    % Recond output
    if T >= T_rec(recond_index)
        thew_rec(:,:,recond_index) = thew ...
            - ((T - T_rec(recond_index)) / dt) .* (thew - thewn);
        pw_rec(:,:,recond_index) = pw...
            - ((T - T_rec(recond_index)) / dt) .* (pw - pwn);
        seff_pc = (thew_rec(:,:,recond_index) - thewr) ...
            ./ (phi - thewr - thenr);
        pc_rec(:,:,recond_index) = pcm .* (seff_pc .^ (-1 ./ m) -1) .^ (1 ./ n);
        pn_rec(:,:,recond_index) = pw_rec(:,:,recond_index) ...
            + pc_rec(:,:,recond_index);
        fprintf('Simulation Time = %d\n',T_rec(recond_index));
        recond_index = recond_index + 1;
        if T >= T_rec(end)
            break
        end
    end

    % Update
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

end

% Recond operation time
fprintf('Operation Time = %d\n',toc);

%% Plot
figure
pl = [2:61];
zp = [99.5:-1:40.5]';
T_n = {'T = 180 s' 'T = 360 s' 'T = 720 s' 'T = 1080 s' 'T = 1440 s'};
hold on
for ii = 1:5
    plot(zp,thew_rec(pl,2,ii),'DisplayName',T_n{ii})
end
hold off
legend
xlim([0 100])
view(90,-90)

