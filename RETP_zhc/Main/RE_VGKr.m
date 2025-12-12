function [Kr] = RE_VGKr(MVG_option,S,theta_r,theta_s,hc_mvg,a,n,m)

% From saturation to relative permeability (VG/MVG)

if MVG_option == 0
    Kr = S.^0.5 .* (1-(1-S.^(1./m)).^m).^2;
end

if MVG_option == 1
    theta_s_mvg = theta_r + (theta_s - theta_r) .* (1 + abs(a .* hc_mvg) .^ n) .^ m;
    S_ori = (theta_s_mvg - theta_r) ./ (theta_s - theta_r) .* S;

    F_1 = (1 - ((theta_s - theta_r) ./ (theta_s_mvg - theta_r)) .^ (1 ./ m)) .^ m;
    F = (1 - S .^ (1 ./ m)) .^ m;

    Kr = S .^ 0.5 .* ((1 - F) ./ (1 - F_1)) .^ 2;
    Kr(S_ori >= 1) = 1;
end

end

