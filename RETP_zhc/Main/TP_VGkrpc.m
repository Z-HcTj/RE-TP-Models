function [pc,krw,krn] = TP_VGkrpc(MVG_option,seff,n,m,a,pcm,hc_mvg,phi,thewr,thenr,rhow0,g)

if MVG_option == 0
    pc = pcm .* (seff .^ (-1 ./ m) -1) .^ (1 ./ n);
    krw = seff .^ 0.5 .* (1 - (1 - (seff .^ (1 ./ m))) .^ m) .^ 2;
    krn = (1 - seff) .^ 0.5 .* ((1 - seff .^ (1 ./ m)) .^ m) .^ 2;
end

if MVG_option == 1
    thews_mvg = thewr + (phi - thewr - thenr) .* (1 + abs(a .* hc_mvg) .^ n) .^ m;
    seff_aft = (phi - thewr - thenr) ./ (thews_mvg - thewr - thenr) .* seff;
    seff_aft_n = (phi - thewr - thenr) ./ (thews_mvg - thewr - thenr) .* (1 - seff);

    F_1 = (1 - ((phi - thewr - thenr) ./ (thews_mvg - thewr - thenr)) .^ (1 ./ m)) .^ m;
    F = (1 - seff_aft .^ (1 ./ m)) .^ m;
    F_n = (1 - seff_aft_n .^ (1 ./ m)) .^ m;

    krw = seff .^ 0.5 .* ((1 - F) ./ (1 - F_1)) .^ 2;
    krn = (1 - seff) .^ 0.5 .* (1 - F_n) .^ 2;
    pc = pcm .* (seff_aft .^ (-1 ./ m) -1) .^ (1 ./ n);

    krw(seff >= 1) = 1;
    krn(seff >= 1) = 0;
    pc(seff >= 1) = rhow0 * g * abs(hc_mvg);

end

end

