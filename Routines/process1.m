function [outNM,outM] = process1(temp,k,PLT,j)
% prende in ingresso un crop dell'img di partenza e la sviluppa
% outNM{1} = out;
% outNM{2} = str;
% outNM{3} = index;
% outNM{4} = outpt;

%%
if nargin == 1
    k = 1.32;
    PLT = 0;
end

if nargin == 2
    PLT = 0;
end

%% MEDIA 
med = 1/9*ones(3);
[temp_med] = filt2plot(temp, med, 0, 'med');
%temp_med = temp;
%% KW5
%temp_k = Kuwahara(temp_med,9); 
% temp_k = uint8(temp_k);
temp_k = temp_med;
%% CS
temp_cs = CS_fast(temp_k);
% temp_cs = temp_k;
%% SHEAR CORRECTION
[ ~,tempf3] = shear_correction( temp_cs,3,'HE');
% tempf3 = temp_cs;
%% CURVE

lumi_t = mean(temp,    2);
lumi_m = mean(temp_med,2);
lumi_k = mean(temp_k,  2);
lumi_cs = mean(temp_cs,2);
lumif3 = mean(tempf3,  2);


%% NON MEDIATA
out  = uint8(255*(tempf3 > k*mean(lumif3)));

% marco i punti trovati
index = [];

str = 0; % numero banda
b = 0; % trovato bianco?
for w = 1:size(out,1)-1
       
   if mean(out(w,:)) > 0.3*255
        if b == 0
            ws = w;
        end
        b = 1;      
    else
        b = 0;
    end
        
    if (mean(out(w+1,:)) > 0.3*255) == (mean(out(w,:)) > 0.3*255) % non cambio
        
    else
        if b == 1
            we = w;
            index = [index,round(0.5*(ws+we))];
            str = str+1;
        end
    end
    
    if w == size(out,1)-1
        if b == 1
            %fprintf('filo a contatto con profilo\n')
            we = w+1;
            index = [index,round(0.5*(ws+we))];
            str = str+1;
        end
    end
            
end

outpt = zeros(size(out));
outpt(index,round(size(out,2)/2)) = 255;

outNM{1} = out;
outNM{2} = str;
outNM{3} = index;
outNM{4} = outpt;
%% MEDIATA
% size(tempf3)
% size((lumif3 > k*mean(lumif3)))

out_m = uint8(255*(ones(size(tempf3,1),1).*(lumif3 > k*mean(lumif3))));
out_m = repmat(out_m,1,size(tempf3,2));
% marco i punti trovati
indexM = [];

strM = 0; % numero banda
b = 0; % trovato bianco?
for w = 1:size(out_m,1)-1
       
   if out_m(w,1) > 0
        if b == 0
            ws = w;
        end
        b = 1;      
    else
        b = 0;
    end
        
    if out_m(w+1,1) == out_m(w,1) % non cambio
        
    else
        if b == 1
            we = w;
            indexM = [indexM,round(0.5*(ws+we))];
            strM = strM+1;
        end
        
        
    end
end
    
outPT = zeros(size(out_m));
outPT(indexM,round(size(out_m,2)/2)) = 255;

outM{1} = out_m;
outM{2} = strM;
outM{3} = indexM;
outM{4} = outPT;
%% PLOT
if PLT == 1
    np = 6;
    spl = 1;
    fig_name = figure;
    set(fig_name,'Position',[50 50 1324 728])
    
    subplot(1,np,spl)
    imshow(temp);
    title(sprintf('j = %d',j))
    
    spl = spl+1;
    subplot(1,np,spl)
    imshow(temp_med);
    title('Media 3x3')
    
%    spl = spl+1;
%    subplot(1,np,spl)
%    imshow(temp_k);
%    title('Kw5')
    
    spl = spl+1;
    subplot(1,np,spl)
    imshow(temp_cs);
    title('CS')
    
    spl = spl+1;
    subplot(1,np,spl)
    imshow(tempf3);
    title('Shear')
    
    spl = spl+1;
    subplot(1,np,spl);
    plot(lumi_t,[size(lumi_t,1):-1:1],'b');
    hold on
    %plot(lumi_m,[size(lumi_t,1):-1:1],'r');
    %plot(lumi_k,[size(lumi_t,1):-1:1],'r');
    %plot(lumi_cs,[size(lumi_t,1):-1:1],'c');
    plot(lumif3,[size(lumi_t,1):-1:1],'r');
    grid on
    legend('normale','PP')%,'Location','NorthEastOutside')
    
    spl = spl+1;
    subplot(1,np,spl);
    imshow(out_m);
    
end
dummy = 0;
end