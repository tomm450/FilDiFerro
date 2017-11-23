function [ lumib,tempf3] = shear_correction( tempf2,fract,opt )
j = 1; % ADATTAMENTO SINTASSI, NON RISCRIVO TUTTO...

y_camp = ceil(size(tempf2,1)/fract);
lumi(:,j) = mean(tempf2,2); 
    
% if opt2 == 'm'
%     lumiup = max(lumi(1:y_camp,j));
%     lumidwn = mean(lumi(end-y_camp,j));
% else
    
    lumiup = mean(lumi(1:y_camp,j));
    lumidwn = mean(lumi(end-y_camp,j));
    
% end
    [lumimin,imin] = min([lumiup,lumidwn]);
    if imin == 1
        lumimax = lumidwn;
        k = linspace(1,lumimin/lumimax,size(lumi,1));
    else
        lumimax = lumiup;
        k = linspace(lumimin/lumimax,1,size(lumi,1));
    end
%    MATLAB 2016       
%    tempf3 = double(tempf2).*k';

%   MATLAB 2015
    tempf3 = zeros(size(tempf2));
    
    for w = 1: size(tempf2,2)
       tempf3(:,w) = double(tempf2(:,w)).*k';
    end
        
    tempf3 = uint8(tempf3);
    
    if opt == 'CS'
        tempf3 = CS_fast(tempf3,0);
    elseif opt == 'HE'
        tempf3 = histeq(tempf3);
    elseif opt == 'no'
        % niente
    else 
        error('opt sbagliato!')
    end
    
    
    
    lumib(:,j) = mean(tempf3,2);


end

