% 2D Richards equation
% Pressure head boundary & Flowrate boundary
% Model parameters

clear
clc

% Domain (zooming only)
% Location of points
x = [2.5:5:497.5]' * 0.01;
z = [-2.5:-5:-297.5]' * 0.01;
% Number of grids (row and column)
nx = length(x);
nz = length(z);
% Grid size
xgs = 0.05 * ones(length(z),length(x));
zgs = 0.05 * ones(length(z),length(x));
% Distance between points
dxp = 0.05 * ones(length(z) + 2,length(x) + 1);
dzp = 0.05 * ones(length(z) + 1,length(x) + 2);
% Time step size
dt = 0.001;
% Tiny step for der
delta = -1 * 10 ^ -5;
% Soil parameters
n = zeros(nz+2,nx+2);
a = zeros(nz+2,nx+2);
% Specific storage
Ss = 1 * 10 ^ -5;
% Porosity
phi = zeros(nz+2,nx+2);
% Gracity option
gr = 1;
% Saturated conductivity
Ks = zeros(nz+2,nx+2);
% Residual water content
theta_r = zeros(nz+2,nx+2);
% Saturated water content
theta_s = zeros(nz+2,nx+2);
% Set homosoil
% Part 1
% Clay
n(1:21,:) = 1.3954;
a(1:21,:) = 1.04;
theta_r(1:21,:) = 0.1060;
theta_s(1:21,:) = 0.4686;
phi(1:21,:) = 0.4686;
Ks(1:21,:) = 13.1 * 0.01 / 86400;
% Sand
n(1:21,32:71) = 2.239;
a(1:21,32:71) = 2.8;
theta_r(1:21,32:71) = 0.0286;
theta_s(1:21,32:71) = 0.3658;
phi(1:21,32:71) = 0.3658;
Ks(1:21,32:71) = 541 * 0.01 / 86400;

% Part 2
% Sand
n(22:41,:) = 2.239;
a(22:41,:) = 2.8;
theta_r(22:41,:) = 0.0286;
theta_s(22:41,:) = 0.3658;
phi(22:41,:) = 0.3658;
Ks(22:41,:) = 541 * 0.01 / 86400;
% Clay
n(22:41,32:71) = 1.3954;
a(22:41,32:71) = 1.04;
theta_r(22:41,32:71) = 0.1060;
theta_s(22:41,32:71) = 0.4686;
phi(22:41,32:71) = 0.4686;
Ks(22:41,32:71) = 13.1 * 0.01 / 86400;

% Part 3
% Clay
n(42:62,:) = 1.3954;
a(42:62,:) = 1.04;
theta_r(42:62,:) = 0.1060;
theta_s(42:62,:) = 0.4686;
phi(42:62,:) = 0.4686;
Ks(42:62,:) = 13.1 * 0.01 / 86400;
% Sand
n(42:62,32:71) = 2.239;
a(42:62,32:71) = 2.8;
theta_r(42:62,32:71) = 0.0286;
theta_s(42:62,32:71) = 0.3658;
phi(42:62,32:71) = 0.3658;
Ks(42:62,32:71) = 541 * 0.01 / 86400;

m = 1 - 1 ./ n;
% Flowrate
SStoFD = 1;
q = -5 * 0.01 / 86400;
SSq = q ./ dzp(1,2);
qrxb = 41;
qrxf = 60;
% Precision
e = 1 * 10 ^ -3;
% Total simulation time
T = 0;
% Recond time
T_rec = 86400 * 12.5;
recond_index = 1;
% Create output variables
H = zeros(nz+2,nx+2);
Hn = zeros(nz+2,nx+2);

% Initial and boundary conditions
% Saturation boundary (top)
% theta_top = zeros(1,nx+2);
% theta_top(:) = 0.3;
% H(1,:) = -(((theta_top - theta_r) / (phi - theta_r)) .^ (-1/m) - 1) .^ (1/n) / a;
% Hn(1,:) = H(1,:);
% 
% % Initial conditions
% theta_ini = 0.1;
% Hini = -(((theta_ini - theta_r) / (phi - theta_r)) .^ (-1/m) - 1) .^ (1/n) / a;
% H(2:end,:) = Hini;
% Hn(2:end,:) = Hini;

H(:,:) = -500;
Hn = H;

% Output data
theta_recond = zeros(nz+2,nx+2,length(T_rec));
% Recond iteration times
ntt = [];

tic

% simulation
tic;
while 1
    
    % Simulation time
    T = T + dt;
    % Iteration times
    nt = 0;

    while 1

        nt = nt + 1;
        % Der H
        H_delta = H + delta;
        % Effective saturation
        S = (1 + abs(a.*H).^n).^-m;
        S_delta = (1 + abs(a.*H_delta).^n).^-m;
        Sn = (1 + abs(a.*Hn).^n).^-m;
        [cho_pHx, cho_pHy] = find(H >= 0);
        [cho_pHdelx, cho_pHdely] = find(H_delta >= 0);
        [cho_pHnx, cho_pHny] = find(Hn >= 0);
        S(cho_pHx, cho_pHy) = 1;
        S_delta(cho_pHx, cho_pHy) = 1;
        Sn(cho_pHx, cho_pHy) = 1;

        % Water Content
        theta = theta_r + (phi - theta_r) .* S;
        theta_delta = theta_r + (phi - theta_r) .* S_delta;
        theta_n = theta_r + (phi - theta_r) .* Sn;

        % Variable K
        K = Ks .* S.^0.5 .* (1-(1-S.^(1./m)).^m).^2;
        K_delta = Ks .* S_delta.^0.5 .* (1-(1-S_delta.^(1./m)).^m).^2;

        % Inter-grid K
        K_mx = (K(2:end-1,1:end-1) + K(2:end-1,2:end)) / 2;
        K_mxf = (K(2:end-1,1:end-1) + K_delta(2:end-1,2:end)) / 2;
        K_mxb = (K_delta(2:end-1,1:end-1) + K(2:end-1,2:end)) / 2;
        K_mz = (K(1:end-1,2:end-1) + K(2:end,2:end-1)) / 2;
        K_mzf = (K(1:end-1,2:end-1) + K_delta(2:end,2:end-1)) / 2;
        K_mzb = (K_delta(1:end-1,2:end-1) + K(2:end,2:end-1)) / 2;

        % Adjust
        K_mx(:,1) = 0;
        K_mxf(:,1) = 0;
        K_mxb(:,1) = 0;
        K_mx(:,end) = 0;
        K_mxf(:,end) = 0;
        K_mxb(:,end) = 0;
        K_mz(1,:) = 0;
        K_mzf(1,:) = 0;
        K_mzb(1,:) = 0;
        K_mz(end,:) = 0;
        K_mzf(end,:) = 0;
        K_mzb(end,:) = 0;

        if SStoFD == 1
            K_mz(1,qrxb:qrxf) = 0;
            K_mzf(1,qrxb:qrxf) = 0;
            K_mzb(1,qrxb:qrxf) = 0;
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

        % Array b
        % Flowrate boundary
        if SStoFD == 1
            F(1,qrxb:qrxf) = F(1,qrxb:qrxf) + SSq;
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
        if nt >= 10
            H = Hn;
            T = T - dt;
            dt = 0.8 * dt;
            break
        end
    end

    % Water content
    S= (1 + abs(a.*H).^n).^-m;
    theta = theta_r + (theta_s - theta_r) .* S;

    % Recond theta_change_sum
    if T >= T_rec(recond_index)
        theta_n = theta_r + (theta_s - theta_r) ...
            .* (1 + abs(a .* Hn) .^ n) .^ -m;
        theta_recond(:,:,recond_index) = theta ...
            - ((T - T_rec(recond_index)) / dt) .* (theta - theta_n);
        fprintf("Simulation time = %4.2f", T_rec(recond_index))
        recond_index = recond_index + 1;
        if T >= T_rec(end)
            break
        end
    end

    if nt <= 3
        dt = 1.2 * dt;
        if dt >= 3000
            dt = 3000;
        end
    end

    % update
    Hn = H;
    ntt = [ntt; nt];

end

% Marching
% figure
% x = -[0.5:1:99.5];
% wc1 = [0.0825 ;0.0825 ;0.0825 ;0.165 ;0.165 ;0.165 ;0.2475 ;0.2475 ;0.2475];
% z1 = [74.58 ;61.32; 39.02; 75.20 ;62.13 ; 40.04; 77.12 ;64.54 ;42.80] - 100;
% 
% plot(x,theta_recond(2:101,2,1),x,theta_recond(2:101,2,2),z1,wc1,'*')
% view(90,-90)
figure
v = [-400 -1 ]
contour(x,z,H(2:end-1,2:end-1),v)


toc


