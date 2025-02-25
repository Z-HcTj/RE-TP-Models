function [thew_tran] = pwtranthewVG(pw_inp,phi,thewr,thenr,n,m,a,rhow0,g)
% Water phase pressure converted to water content (air pressure is equal to zero)

thew_tran = (1 + abs(a .* pw_inp ./ rhow0 / g).^n).^-m ...
    .* (phi - thewr - thenr) + thewr;
thew_sat = phi - thewr - thenr;
thew_tran(pw_inp >= 0) = thew_sat(pw_inp >= 0);

end