function [out] = r2rgb(in,max_in,min_in,lvl)

l = round(lvl/2);
s = linspace(0,1,l)';
z = zeros(l,1);

rgb = [ z s flipud(s); s(2:end) flipud(s(1:end-1)) z(1:end-1)];
l_eff = size(rgb,1);



in_norm = (in - min_in)/(max_in-min_in);


if isnan(round(in_norm*l_eff)) == 1
    out = [0 0 0];
elseif (in_norm*l_eff) == 0
    out = [0 0 1];
else
    
    index = min([round(in_norm*l_eff),l_eff]);
    out = rgb(index,:);
end

out = 255*out;
