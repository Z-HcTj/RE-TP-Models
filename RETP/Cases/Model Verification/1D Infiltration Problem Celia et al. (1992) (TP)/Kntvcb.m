function Kn = Kntvcb(Kas,Aa,Ba,hc,g,rhon0,rhow0)

Kn = ((Kas * Aa) ./ (Aa + (hc * 100) .^ Ba)) * rhon0 / (g * rhow0);

end

