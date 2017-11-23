function [mask] = mirino(TL,BR,dim,thickness)
%MIRINO Summary of this function goes here
%   Detailed explanation goes here
if nargin == 3
    thickness = 0;
end

mask = zeros(dim);

for r = 1:dim(1)
%    for c  = 1:dim(2)
        if r < TL(1) - 2*thickness
            % do nothing
        elseif r >= TL(1)-2*thickness && r <= TL(1)
        %elseif r == TL(1)
            mask(r,TL(2)-thickness:BR(2)+thickness) = 255;
        elseif r > TL(1) && r < BR(1)
            mask(r,TL(2)-thickness:TL(2)+thickness) = 255;
            mask(r,BR(2)-thickness:BR(2)+thickness) = 255;
        
        elseif r >= BR(1) && r <= BR(1)+2*thickness    
        %elseif r == BR(1)
            mask(r,TL(2)-thickness:BR(2)+thickness) = 255;
        elseif r > BR(1) +2*thickness
            % do nothing
        end
end

% mask_rgb(:,:,1) = mask;
% mask_rgb(:,:,2) = zeros(size(mask));
% mask_rgb(:,:,3) = zeros(size(mask));

mask = uint8(mask);
            
           
end

