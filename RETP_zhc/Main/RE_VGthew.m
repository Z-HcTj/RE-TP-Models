function [theta] = RE_VGthew(MVG_option,theta_r,theta_s,S,hc_mvg,a,n,m)

% From saturation to water content (VG/MVG)

if MVG_option == 0
    theta = theta_r + (theta_s - theta_r) .* S;
end

if MVG_option == 1
    theta_s_mvg = theta_r + (theta_s - theta_r) .* (1 + abs(a .* hc_mvg) .^ n) .^ m;
    S_ori = (theta_s_mvg - theta_r) ./ (theta_s - theta_r) .* S;
    theta = theta_r + (theta_s_mvg - theta_r) .* S;
    theta(S_ori>=1) = theta_s(S_ori>=1);
end

