function [i_new] = CS_fast(itemp,PLT)
  % SCRIPT INTESO PER RGB O TENSORI DI RGB
  if nargin == 1
      PLT = 0;
  end
  
  % find the min. value of pixel
  if size(itemp,3) == 1
         rmin = min(min(itemp)); 
  elseif size(itemp,3) > 1
         PLT = 0;
         rmin = min(min(min(itemp))); 
  end
  
  if size(itemp,3) == 1
         rmax = max(max(itemp)); 
  elseif size(itemp,3) > 1
         PLT = 0;
         rmax = max(max(max(itemp))); 
  end     

%   %% UFF
%   m = 255/(rmax - rmin);  % find the slope of line joining point (0,255) to (rmin,rmax)
%   c = 255 - m*rmax;       % find the intercept of the straight line with the axis
%   i_new = m.*itemp + c;   % transform the image according to new slope
%   
%   i_new = uint8(i_new);
%   %% TOM
    if rmin >=0 
        itemp = double(uint8(itemp) - uint8(rmin));
    else
        itemp = itemp - rmin;
    end
    
       
   rmax2 = rmax-rmin;
   
   i_new = itemp.*(255/double(rmax2));
   i_new = uint8(i_new);
  
  
  
  if size(i_new) == size(itemp)
      % ok
%       size(i_new)
%       size(itemp)
  else 
      error('cambio size...')
  end
  
  
  if PLT == 1
      figure
      subplot(1,2,1)
      imshow(itemp)
      title('W/O CS')% display original image
      subplot(1,2,2)
      imshow(i_new);   % display transformed image
      title('W/ CS')
      
      figure
      subplot(1,2,1)
      imhist(itemp)
      title('W/O CS')% display original image
      subplot(1,2,2)
      imhist(i_new);   % display transformed image
      title('W/ CS')
  end