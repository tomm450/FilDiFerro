function [b_f,diff] = filt2plot(a, FILT, PLT, string)

% convoluzione
b_f_temp = filter2(FILT,a,'same');
% riporto su 8 bit
b_f_temp = CS_fast(b_f_temp,0);
%converto a immagine
b_f = uint8(b_f_temp);

if nargin == 2
    PLT = 0;
    string = '';
elseif nargin == 3
    string = '';
end

if nargout == 2
  
    b_f_double = double(b_f);
    a_double = double(a);
    diff_double = b_f_double - a_double;
    
    diff_temp = CS_fast(diff_double,0);
    diff = uint8(diff_temp);
   
end


if PLT >= 1
    
    figure
    imshow(b_f)
    title(sprintf('Immagine Filtra %s',string))
    %subplot(1,3,3)
    if nargout == 2
        figure
        imshow(diff)
        title(sprintf('Differenza Filtra %s',string))
    end
end

end