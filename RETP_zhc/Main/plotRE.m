 %% 2D Richards-based model
clear
clc

%% Input parameters
id = load('Inputdata.mat');
% Simulation time
tt = load('simulation time.mat');
% Output save
load('mydata.mat');

% Domain
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
% Initial timestep
dt = id.dt;
% Time step criterion
dt_max = id.dt_max;
dt_ac = id.inu_ac;
dt_limit = id.inu_limit;
% Small perturbation
delta = id.dpw / (id.rhow0 * id.g);
% Gracity option
gr = id.gr;
% Specific storage
Ss = id.Ss;
% Soil parameters
n = 1.5;
a = id.a;
theta_r = id.thewr;
theta_s = id.thews;
phi = id.phi;
Ks = id.k * (id.rhow0 * id.g / id.vw);
m = id.m;
% Converge parameter
e = id.epw / (id.rhow0 * id.g);
% Initial simulation time
T = 0;
% Recond at set time points
T_rec = tt.T_simulation * id.T_ran;
recond_index = 1;
theta_recond = zeros(nz+2,nx+2,length(T_rec));
pw_recond = zeros(nz+2,nx+2,length(T_rec));
uw_recond = zeros(nz+1,nx+2,length(T_rec));
% Create output variables
H = zeros(nz+2,nx+2);
Hn = zeros(nz+2,nx+2);

% Initial and boundary conditions
% Initial condition (including pressure head boundary)
H = id.pw / (id.rhow0 * id.g);
H = [-1.01:0.01:0]';
Hn = H;
% Flowrate boundary
SStoFD = id.SStoFD;
if SStoFD == 1
    q = id.q / id.rhow0;
    SSq = q ./ dzp(1,2);
    qrxb = id.SSq_left;
    qrxf = id.SSq_right;
    UP_stage1_option = SStoFD;
    UP_stage2_option = 0;
    pwbmi = H(2,qrxb+1:qrxf+1) * (id.rhow0 * id.g);
    pwbmm = -H(2,qrxb+1:qrxf+1) * (id.rhow0 * id.g);
    epw_BM = id.epw_BM;
    curan = qrxb+1:qrxf+1;
    uwdet = 0;
end
% Ponding time (if exist)
T_re_pond = 0;
% Recond iteration times
ntt = [];
% Boundary options 
BD_leftclose = id.left_close;
BD_rightclose = id.right_close;
BD_topclose = id.top_close;
BD_bottomclose = id.bottom_close;

% Mass balance valuables
duw = 0;
uwn = 0;

% Valuable for recording whether 
% a time point is incorrect or not
outtime = [];

% MVG model application
MVG_option = 1;
hc_mvg = id.hc_mvg;

%% simulation
% Operation time
tic;


    % Simulation time
    T = T + dt;
    % Iteration times for every NR iteration
    nt = 0;
    % Boundary Ponding setup
    if SStoFD == 1
        if UP_stage2_option == 1
            H_up2_ini = H(1,curan);
        end
    end


        nt = nt + 1;
        % Flowrate top boundary setup
        if SStoFD == 1
            % UP_stage1 -- no ponding occurs, update boundary pressure head
            if UP_stage1_option == 1
                thew = (1 + abs(a.* H).^n).^-m .* (theta_s - theta_r) + theta_r;
                [thewtop, pwtop] = getupperboundaryBM(id.rhow0,pwbmi,pwbmm,H(2,curan) * (id.rhow0 * id.g)...
                    ,thew(2,curan),q*id.rhow0,dzp(1,curan),n(1,curan),m(1,curan),id.k(1,curan),id.vw,a(1,curan),...
                    id.g,gr,theta_r(1,curan),0,phi(1,curan),epw_BM);
                H(1,curan) = pwtop / (id.rhow0 * id.g);
                if H(1,curan) > 0
                    UP_stage1_option = 0;
                    UP_stage2_option = 1;
                    H(1,curan) = 0;
                    H_up2_ini = H(1,curan);
                    T_re_pond = T;
                end
            end
            % UP_stage2 -- ponding occurs, update boundary pressure head
            if UP_stage2_option == 1
                S_up2 = (1 + abs(a.*H).^n).^-m;
                S_up2(H > 0) = 1;
                K_up2 = Ks .* S_up2.^0.5 .* (1-(1-S_up2.^(1./m)).^m).^2;
                K_up2_mid = (K_up2(1,2:end-1) + K_up2(2,2:end-1)) / 2;
                H(1,curan) = H_up2_ini;
                H_dt_up2 = (-q - K_up2_mid * (H(1,curan) - H(2,curan)) / dzp(1,curan) - K_up2_mid) ./ ...
                    (1 / dt + K_up2_mid / dzp(1,curan));
                H(1,curan) = H_up2_ini + H_dt_up2;

                % if H_bc < 0, no ponding occurs
                if H(1,curan) < 0
                    UP_stage1_option = 1;
                    UP_stage2_option = 0;
                end
            end
        end

        % Der H
        H_delta = H + delta;

        % Effective saturation
        S = (1 + abs(a.*H).^n).^-m;
        S_delta = (1 + abs(a.*H_delta).^n).^-m;
        Sn = (1 + abs(a.*Hn).^n).^-m;
        S(H > 0) = 1;
        S_delta(H > 0) = 1;
        Sn(H > 0) = 1;

        % Water Content
        theta = RE_VGthew(MVG_option,theta_r,theta_s,S,hc_mvg,a,n,m);
        theta_delta = RE_VGthew(MVG_option,theta_r,theta_s,S_delta,hc_mvg,a,n,m);
        theta_n = RE_VGthew(MVG_option,theta_r,theta_s,Sn,hc_mvg,a,n,m);

        % Variable K
        Kr = RE_VGKr(MVG_option,S,theta_r,theta_s,hc_mvg,a,n,m);
        K = Ks .* RE_VGKr(MVG_option,S,theta_r,theta_s,hc_mvg,a,n,m);
        K_delta = Ks .* RE_VGKr(MVG_option,S,theta_r,theta_s,hc_mvg,a,n,m);

        figure
        plot(H,Kr)
        figure
        plot(H,theta)
        figure
        plot(S,Kr)