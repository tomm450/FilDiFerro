function [out, nonzero ] = TH01( in,th,plt )

if nargin == 2
    plt = 0;
end


out = 255*( in > th);
out = uint8(out);

nonzero = (sum(sum(out)))/255;

if plt == 1
    figure; imshow(out); 
    title(sprintf('TH = %d ; elementi diveri da 0 -> %d',th,nonzero))

end

