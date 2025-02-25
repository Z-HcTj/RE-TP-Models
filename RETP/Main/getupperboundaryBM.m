function [thewtop, pwtop, qout] = getupperboundaryBM(rhow0,pw1,pw1_ex,pw2,thew2,qw,dz,n,m,k,vw,a,g,gr,thewr,thenr,phi,epw)
% 获取每次迭代中的顶部含水量，作为迭代的输入条件
% 以上一次迭代输出的水头值和含水量（第二个网格）作为输入，求取此次迭代过程中应该具备边界含水量
% 二分法求解(仅限于顶部气压为0的情况)

% 初始值
while 1
    if pw1 < 0
        thew1 = pwtranthewVG(pw1,phi,thewr,thenr,n,m,a,rhow0,g);
    else
        thew1 = phi - thenr;
    end
    seff1 = (thew1 - thewr) ./ (phi - thewr - thenr);
    seff2 = (thew2 - thewr) ./ (phi - thewr - thenr);
    seff1(seff1 > 1) = 1;
    seff2(seff2 > 1) = 1;
    krw1 = seff1 .^ 0.5 .* (1 - (1 - (seff1 .^ (1 ./ m))) .^ m) .^ 2;
    krw2 = seff2 .^ 0.5 .* (1 - (1 - (seff2 .^ (1 ./ m))) .^ m) .^ 2;
    K1 = rhow0 * krw1 .* k / vw;
    K2 = rhow0 * krw2 .* k / vw;
    Km = (K1 + K2) / 2;
    F = qw - Km .* ((pw2 - pw1) ./ dz - rhow0 * g * gr);
    qout1 = - Km .* ((pw2 - pw1) ./ dz - rhow0 * g * gr);


    % 极限值
    if pw1_ex < 0
        thew1_ex = (1 + abs(a .* pw1_ex ./ rhow0 / g).^n).^-m .* ...
            (phi - thewr - thenr) + thewr;
    else
        thew1_ex = phi - thenr;
    end
    seff1_ex = (thew1_ex - thewr) ./ (phi - thewr - thenr);
    seff1_ex(seff1_ex > 1) = 1;
    krw1_ex = seff1_ex .^ 0.5 .* (1 - (1 - (seff1_ex .^ (1 ./ m))) .^ m) .^ 2;
    K1_ex = rhow0 * krw1_ex .* k / vw;
    Km_ex = (K1_ex + K2) / 2;
    F_ex = qw - Km_ex .* ((pw2 - pw1_ex) ./ dz - rhow0 * g * gr);
    qout2 = - Km_ex .* ((pw2 - pw1_ex) ./ dz - rhow0 * g * gr);


    % 中间值
    pw1_mid = (pw1 + pw1_ex) / 2;
    if pw1_mid < 0
        thew1_mid = (1 + abs(a .* pw1_mid ./ rhow0 / g).^n).^-m .* ...
            (phi - thewr - thenr) + thewr;
    else
        thew1_mid = phi - thenr;
    end
    seff1_mid = (thew1_mid - thewr) ./ (phi - thewr - thenr);
    seff1_mid(seff1_mid > 1) = 1;
    krw1_mid = seff1_mid .^ 0.5 .* (1 - (1 - (seff1_mid .^ (1 ./ m))) .^ m) .^ 2;
    K1_mid = rhow0 * krw1_mid .* k / vw;
    Km_mid = (K1_mid + K2) / 2;
    F_mid = qw - Km_mid .* ((pw2 - pw1_mid) ./ dz - rhow0 * g * gr);

    % 二分法
    for ii = 1:length(F)
        if F(ii) * F_mid(ii) < 0
            pw1_ex(ii) = pw1_mid(ii);
        elseif F_mid(ii) * F_ex(ii) < 0
            pw1(ii) = pw1_mid(ii);
        end
    end

    if abs(max(abs(pw1)) - max(abs(pw1_ex))) < epw
        break
    end

end

pwtop = (pw1 + pw1_ex) / 2;
if pwtop < 0
    thewtop = (1 + abs(a .* pwtop ./ rhow0 / g).^n).^-m .* ...
        (phi - thewr - thenr) + thewr;
else
    thewtop = phi - thenr;
end

 
