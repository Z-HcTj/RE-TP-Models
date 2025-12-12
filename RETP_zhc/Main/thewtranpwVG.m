function [pw_tran] = thewtranpwVG(thew_inp,phi,thewr,thenr,n,m,a,rhow0,g)
% Water content converted to water phase pressure (air pressure is equal to zero)

pw_tran = -(((thew_inp - thewr) ./ (phi - thewr - thenr)) ...
    .^ (-1 ./ m) - 1) .^ (1 ./ n) ./ a * rhow0 * g;

end

