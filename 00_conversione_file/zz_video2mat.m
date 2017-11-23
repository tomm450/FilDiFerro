% %%% DA LANCIARE SOLO LA PRIMA VOLTA

clc
clear all
close all

%%
folder_in  = './IN';
folder_out = './'; 

addpath ./Routines
addpath(folder_in);
addpath(folder_out);

video_n = 'MVI_9911';
ext   = 'MOV';

video_name = strcat(video_n,'.',ext); 

fps = 50;

PROC = 0; % 0 = imresize; 1 = binned

if PROC == 0
    h_out = 640;
elseif PROC == 1
    r = 4;
else 
    error('PROC ~= 0 | 1')
end

%%%%
IN = importfile(video_name);

if PROC == 0
    scale_f = h_out/size(IN,1);
end

screen_to_load = size(IN,4);

% probe to dimension
if PROC == 0
    probe = imresize(rgb2gray(uint8(IN(:,:,:,1))),scale_f);
elseif PROC == 1
    probe = binning(rgb2gray(uint8(IN(:,:,:,1))),0,r );
end
size(rgb2gray(uint8(IN(:,:,:,1))),2)
h_out = size(probe,2)
pause(3)
sdeng

T = (zeros([size(probe,1),size(probe,2),screen_to_load]));
res = double(T);

clear probe

for k = 1:screen_to_load

    [dummy] = process_bar(k,screen_to_load,'FASE 1 - scansione filmato');
        
    t0 = (IN(:,:,:,k));
    
    t1 = rgb2gray(t0);
    
    if PROC == 0
        t2 = imresize(t1,scale_f);
    else
        t2 = binning( t1,0,r );
    end
    %t = uint8(rgb2gray(t1));
    
    T(:,:,k) = t2;
    
    
end

Tim = sum(T,3)./k;
Tim = uint8(Tim);




%%%
for k = 1:screen_to_load

    [dummy] = process_bar(k,screen_to_load,'FASE 2 - sottraggo media');
       
    res(:,:,k) = (double(T(:,:,k)) - double(Tim));

end

res = CS_fast(res);

clc
disp('-------------------->|  FASE 3 - salvo')
PAR.T = T;
PAR.Tim = Tim;
PAR.res = res;

save(strcat(folder_out,'/',video_n,'_frame_fps_',num2str(fps),'_hout_',num2str(h_out),'_ext_',ext,'.mat'),'PAR');
imwrite(Tim,strcat(folder_out,'/',video_n,'_media_fps_',num2str(fps),'_hout_',num2str(h_out),'_ext_',ext,'.png'));


clc
disp('-------------------->|  FASE 4 - DONE!')
clear all




