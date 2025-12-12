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
n = id.n;
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
masduw = 0;

% back parameter
ntt_back = 0;
re_outline = 0;

% Valuable for recording whether 
% a time point is incorrect or not
outtime = [];

% MVG model application
MVG_option = id.MVG_option;
hc_mvg = id.hc_mvg;

%% simulation
% Operation time
tic;

while 1
    % Simulation time
    T = T + dt;
    % Iteration times for every NR iteration
    nt = 0;
    % back parameter
    ntt_back = 0;
    % Boundary Ponding setup
    if SStoFD == 1
        if UP_stage2_option == 1
            H_up2_ini = H(1,curan);
        end
    end

    while 1
        nt = nt + 1;
        % Flowrate top boundary setup
        if SStoFD == 1
            % UP_stage1 -- no ponding occurs, update boundary pressure head
            if UP_stage1_option == 1
                thew = (1 + abs(a.* H).^n).^-m .* (theta_s - theta_r) + theta_r;
                [thewtop, pwtop] = getupperboundaryBM(MVG_option,hc_mvg,id.rhow0,pwbmi,pwbmm,H(2,curan) * (id.rhow0 * id.g)...
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
        K = Ks .* RE_VGKr(MVG_option,S,theta_r,theta_s,hc_mvg,a,n,m);
        K_delta = Ks .* RE_VGKr(MVG_option,S_delta,theta_r,theta_s,hc_mvg,a,n,m);

        % Inter-grid K
        K_mx = (K(2:end-1,1:end-1) + K(2:end-1,2:end)) / 2;
        K_mxf = (K(2:end-1,1:end-1) + K_delta(2:end-1,2:end)) / 2;
        K_mxb = (K_delta(2:end-1,1:end-1) + K(2:end-1,2:end)) / 2;
        K_mz = (K(1:end-1,2:end-1) + K(2:end,2:end-1)) / 2;
        K_mzf = (K(1:end-1,2:end-1) + K_delta(2:end,2:end-1)) / 2;
        K_mzb = (K_delta(1:end-1,2:end-1) + K(2:end,2:end-1)) / 2;

        % Boundary close or not
        if BD_leftclose == 1
            K_mx(:,1) = 0;
            K_mxf(:,1) = 0;
            K_mxb(:,1) = 0;
        end
        if BD_rightclose == 1
            K_mx(:,end) = 0;
            K_mxf(:,end) = 0;
            K_mxb(:,end) = 0;
        end
        if BD_topclose == 1
            K_mz(1,:) = 0;
            K_mzf(1,:) = 0;
            K_mzb(1,:) = 0;
        end
        if BD_bottomclose == 1
            K_mz(end,:) = 0;
            K_mzf(end,:) = 0;
            K_mzb(end,:) = 0;
        end

        if id.top_wclose_nopen == 1
            K_mz(1,:) = 0;
            K_mzf(1,:) = 0;
            K_mzb(1,:) = 0;
        end

        % Flowrate boundary setup
        if SStoFD == 1
            if UP_stage1_option == 1
                K_mz(1,qrxb:qrxf) = 0;
                K_mzf(1,qrxb:qrxf) = 0;
                K_mzb(1,qrxb:qrxf) = 0;
            end
        end

        % Calculate F
        % Storage term
        St = Ss * theta(2:end-1,2:end-1) ./ phi(2:end-1,2:end-1) ...
            .* (H(2:end-1,2:end-1) - Hn(2:end-1,2:end-1)) /dt;
        St_delta = Ss * theta_delta(2:end-1,2:end-1) ./ phi(2:end-1,2:end-1) ...
            .* (H_delta(2:end-1,2:end-1) - Hn(2:end-1,2:end-1)) /dt;
        % Time term
        Tt = (theta(2:end-1,2:end-1) - theta_n(2:end-1,2:end-1))...
            / dt;
        Tt_delta = (theta_delta(2:end-1,2:end-1) - theta_n(2:end-1,2:end-1))...
            / dt;
        % Tran term 
        Trt = - 1 ./ xgs .* (K_mx(:,2:end) .* (H(2:end-1,3:end) - H(2:end-1,2:end-1)) ./ dxp(2:end-1,2:end) ...
            - K_mx(:,1:end-1) .* (H(2:end-1,2:end-1) - H(2:end-1,1:end-2)) ./ dxp(2:end-1,1:end-1)) ...
            -1 ./ zgs .* (K_mz(2:end,:) .* ((H(3:end,2:end-1) - H(2:end-1,2:end-1)) ./ dzp(2:end,2:end-1) - gr) ...
            - K_mz(1:end-1,:) .* ((H(2:end-1,2:end-1) - H(1:end-2,2:end-1)) ./ dzp(1:end-1,2:end-1) - gr));
        Trt_delta = - 1 ./ xgs .* (K_mxb(:,2:end) .* (H(2:end-1,3:end) - H_delta(2:end-1,2:end-1)) ./ dxp(2:end-1,2:end) ...
            - K_mxf(:,1:end-1) .* (H_delta(2:end-1,2:end-1) - H(2:end-1,1:end-2)) ./ dxp(2:end-1,1:end-1)) ...
            -1 ./ zgs .* (K_mzb(2:end,:) .* ((H(3:end,2:end-1) - H_delta(2:end-1,2:end-1)) ./ dzp(2:end,2:end-1) - gr) ...
            - K_mzf(1:end-1,:) .* ((H_delta(2:end-1,2:end-1) - H(1:end-2,2:end-1)) ./ dzp(1:end-1,2:end-1) - gr));

        Trt_dxip1 = - 1 ./ xgs .* (K_mxf(:,2:end) .* (H_delta(2:end-1,3:end) - H(2:end-1,2:end-1)) ./ dxp(2:end-1,2:end) ...
            - K_mx(:,1:end-1) .* (H(2:end-1,2:end-1) - H(2:end-1,1:end-2)) ./ dxp(2:end-1,1:end-1)) ...
            -1 ./ zgs .* (K_mz(2:end,:) .* ((H(3:end,2:end-1) - H(2:end-1,2:end-1)) ./ dzp(2:end,2:end-1) - gr) ...
            - K_mz(1:end-1,:) .* ((H(2:end-1,2:end-1) - H(1:end-2,2:end-1)) ./ dzp(1:end-1,2:end-1) - gr));
        Trt_dxid1 = - 1 ./ xgs .* (K_mx(:,2:end) .* (H(2:end-1,3:end) - H(2:end-1,2:end-1)) ./ dxp(2:end-1,2:end) ...
            - K_mxb(:,1:end-1) .* (H(2:end-1,2:end-1) - H_delta(2:end-1,1:end-2)) ./ dxp(2:end-1,1:end-1))...
            -1 ./ zgs .* (K_mz(2:end,:) .* ((H(3:end,2:end-1) - H(2:end-1,2:end-1)) ./ dzp(2:end,2:end-1) - gr) ...
            - K_mz(1:end-1,:) .* ((H(2:end-1,2:end-1) - H(1:end-2,2:end-1)) ./ dzp(1:end-1,2:end-1) - gr));    

        Trt_dzjp1 = - 1 ./ xgs .* (K_mx(:,2:end) .* (H(2:end-1,3:end) - H(2:end-1,2:end-1)) ./ dxp(2:end-1,2:end) ...
            - K_mx(:,1:end-1) .* (H(2:end-1,2:end-1) - H(2:end-1,1:end-2)) ./ dxp(2:end-1,1:end-1)) ...
            -1 ./ zgs .* (K_mzf(2:end,:) .* ((H_delta(3:end,2:end-1) - H(2:end-1,2:end-1)) ./ dzp(2:end,2:end-1) - gr) ...
            - K_mz(1:end-1,:) .* ((H(2:end-1,2:end-1) - H(1:end-2,2:end-1)) ./ dzp(1:end-1,2:end-1) - gr));
        Trt_dzjd1 = - 1 ./ xgs .* (K_mx(:,2:end) .* (H(2:end-1,3:end) - H(2:end-1,2:end-1)) ./ dxp(2:end-1,2:end) ...
            - K_mx(:,1:end-1) .* (H(2:end-1,2:end-1) - H(2:end-1,1:end-2)) ./ dxp(2:end-1,1:end-1)) ...
            -1 ./ zgs .* (K_mz(2:end,:) .* ((H(3:end,2:end-1) - H(2:end-1,2:end-1)) ./ dzp(2:end,2:end-1) - gr) ...
            - K_mzb(1:end-1,:) .* ((H(2:end-1,2:end-1) - H_delta(1:end-2,2:end-1)) ./ dzp(1:end-1,2:end-1) - gr));  

        % F
        F = St + Tt + Trt;
        F_delta = St_delta + Tt_delta + Trt_delta;
        F_der = (F_delta - F) / delta;
        F_der = F_der(:);
        F_der_jp1 = (Trt_dzjp1 - Trt) / delta;
        F_der_jd1 = (Trt_dzjd1 - Trt) / delta;
        F_der_ip1 = (Trt_dxip1 - Trt) / delta;
        F_der_id1 = (Trt_dxid1 - Trt) / delta;
        % Last row remove
        F_der_jp1(end,:) = 0;
        F_der_jp1 = F_der_jp1(:);
        F_der_jp1 = F_der_jp1(1:end-1);
        % First row remove
        F_der_jd1(1,:) = 0;
        F_der_jd1 = F_der_jd1(:);
        F_der_jd1 = F_der_jd1(2:end);
        % Last cloumn remove
        F_der_ip1 = F_der_ip1(:,1:end-1);
        F_der_ip1 = F_der_ip1(:);
        % First column remove
        F_der_id1 = F_der_id1(:,2:end);
        F_der_id1 = F_der_id1(:);

        % Array b (flowrate boundary setup)
        if SStoFD == 1
            if UP_stage1_option == 1
                F(1,qrxb:qrxf) = F(1,qrxb:qrxf) + SSq;
            end
        end
        b = F(:);

        % Matrix A
        A = diag(F_der) ...
            + diag(F_der_jp1,1) ...
            + diag(F_der_jd1,-1) ...
            + diag(F_der_ip1,nz) ...
            + diag(F_der_id1,-nz);
        A = sparse(A);

        % calculate deltaH
        delta_h = A\(-b);
        delta_h = reshape(delta_h,nz,nx);
        H(2:end-1,2:end-1) = H(2:end-1,2:end-1) + delta_h;

        % break condition
        abs_delta_h = abs(delta_h);
        if max(max(abs_delta_h)) < e
            break
        end

        % Timestep
        if nt >= dt_limit
            H = Hn;
            T = T - dt;
            dt = 0.8 * dt;
            ntt_back = 1;
            break
        end
    end

    if ntt_back == 0
        % Flowrate boundary setup
        % In case of unsaturated boundary, the top head value is supplemented by
        if SStoFD == 1
            if UP_stage1_option == 1
                thew = (1 + abs(a.* H).^n).^-m .* (theta_s - theta_r) + theta_r;
                [thewtop, pwtop] = getupperboundaryBM(MVG_option,hc_mvg,id.rhow0,pwbmi,pwbmm,H(2,curan) * (id.rhow0 * id.g)...
                    ,thew(2,curan),q*id.rhow0,dzp(1,curan),n(1,curan),m(1,curan),id.k(1,curan),id.vw,a(1,curan),...
                    id.g,gr,theta_r(1,curan),0,phi(1,curan),epw_BM);
                H(1,curan) = pwtop / (id.rhow0 * id.g);
            end
        end

        % Water content
        S= (1 + abs(a.*H).^n).^-m;
        S(H > 0) = 1;
        theta = RE_VGthew(MVG_option,theta_r,theta_s,S,hc_mvg,a,n,m);

        % Mass balance
        % Parameter
        K_mas = Ks .* RE_VGKr(MVG_option,S,theta_r,theta_s,hc_mvg,a,n,m);
        K_mz_mas = (K_mas(1:end-1,:) + K_mas(2:end,:)) / 2;
        if BD_bottomclose == 1
            K_mz_mas(end,:) = 0;
        end

        % Calculate
        uw = -K_mz_mas .* ((H(2:end,:) - H(1:end-1,:)) ./ dzp(:,:) - 1 * gr);
        duw = duw + sum(uw(1,id.SSq_left+1:id.SSq_right+1) * dt) - sum(uw(end,2:end-1) * dt);
        dmvw = (theta(2:end-1,2:end-1)) .* zgs ...
            - (id.thew(2:end-1,2:end-1)) .* zgs;
        masduw = (sum(sum(dmvw)) - duw);

        % If flowrate boundary
        if SStoFD == 1
            uwdet = uwdet + (uw(1) + q) * dt;
        end

        % Recond theta_change_sum
        if T >= T_rec(recond_index)
            if T < T_rec(end)
                if T>= T_rec(recond_index + 1)
                    outtime = [outtime recond_index];
                end
            end
            S_n_rec = (1 + abs(a .* Hn) .^ n) .^ -m;
            S_n_rec(H >= 0) = 1;
            theta_n = RE_VGthew(MVG_option,theta_r,theta_s,S_n_rec,hc_mvg,a,n,m);
            theta_recond(:,:,recond_index) = theta ...
                - ((T - T_rec(recond_index)) / dt) .* (theta - theta_n);
            pw_recond(:,:,recond_index) = H * 1000 * 9.81 ...
                - ((T - T_rec(recond_index)) / dt) .* 1000 * 9.81 * (H - Hn);
            uw_recond(:,:,recond_index) = uw ...
                - ((T - T_rec(recond_index)) / dt) .* (uw - uwn);
            recond_index = recond_index + 1;
            if T >= T_rec(end)
                break
            end
        end
    end

    % Timestep
    if nt <= dt_ac
        dt = 1.2 * dt;
        if dt >= dt_max
            dt = dt_max;
        end
    end

        % not convergence
    if dt < id.dt * 0.1
        re_outline = 1;
        break
    end
%     if sum(ntt == 10) > 20
%         re_outline = 1;
%         break
%     end   

    % Update
    uwn = uw;
    Hn = H;
    ntt = [ntt; nt];

end

% Recond Operation time
t_oper = toc;

%% Save to ''mydata.mat''
ds.re_outtime = outtime;
ds.re_nt = ntt;
ds.trec = id.T_ran; 
ds.tsim = T;
ds.re_timeoper = t_oper;
ds.re_thew = theta_recond;
ds.re_pw = pw_recond;
ds.re_uw = uw_recond;
ds.re_wmass = abs(masduw / duw);
ds.T_re_pond = T_re_pond;
ds.re_outline = re_outline;
save('mydata.mat', 'ds');