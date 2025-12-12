function [thew_tran] = pwtranthewVG(MVG_option,hc_mvg,pw_inp,phi,thewr,thenr,n,m,a,rhow0,g)
% Water phase pressure converted to water content (air pressure is equal to zero)
if MVG_option == 0
    thew_tran = (1 + abs(a .* pw_inp ./ rhow0 / g).^n).^-m ...
        .* (phi - thewr - thenr) + thewr;
    thew_sat = phi - thewr - thenr;
    thew_tran(pw_inp >= 0) = thew_sat(pw_inp >= 0);
end

if MVG_option == 1
    thews_mvg = thewr + (phi - thewr - thenr) .* (1 + abs(a .* hc_mvg) .^ n) .^ m;
    seff = (1 + abs(a .* pw_inp ./ rhow0 / g).^n).^-m;
    seff_ori = (thews_mvg - thewr - thenr) ./ (phi - thewr - thenr) .* seff;
    thew_tran = thewr + thenr + (thews_mvg - thewr - thenr) .* seff;
    thew_sat = phi - thewr - thenr;
    thew_tran(seff_ori>=1) = thew_sat(seff_ori>=1);
end

end