% %%% DA LANCIARE SOLO LA PRIMA VOLTA

clc
clear all
close all

%%
folder_in  = './IN/MVI_9911_bmp'; 
folder_out = '.'; %'SHED50FPS'; %'./SHED50FPS';

addpath ./Routines
addpath(folder_in);
addpath(folder_out);

% video_n = 'MVI_9921';
% INstr = 'vlcsnap-00';
% ext   = 'png';

video_n = 'shed2';
INstr = 'image-';
ext   = 'bmp';

screen_to_load = 115;
offset         = 0;

%video_name = strcat(video_n,'.',ext); 

fps = 50;
h_out = 640;


%%
% probe to dimension
%probe =  imread(strcat(folder_in,'/',video_n,'/',INstr,sprintf('%3.3d',offset+1),'.',ext));   
probe =  imread(strcat(folder_in,'/',INstr,sprintf('%3.3d',offset+1),'.',ext));   
probe = probe(100:end,1:1000,:);
h_in  = size(probe,1);
scale_f = h_out/h_in;
%probe = imresize(probe,scale_f);
probe = binning(probe,0,0);
T = (zeros(size(probe,1),size(probe,2),screen_to_load));
res = double(T);

clear probe


for k = 1:screen_to_load

    [dummy] = process_bar(k,screen_to_load,strcat('FASE 1 - scansione filmato')); %: IMG: ',folder_in,'/',INstr,sprintf('%3.3d',offset+k),'.',ext));
    %IN =  imread(strcat(folder_in,'/',video_n,'/',INstr,sprintf('%3.3d',offset+k),'.',ext));    
    IN =  imread(strcat(folder_in,'/',INstr,sprintf('%3.3d',offset+k),'.',ext));
    IN = IN(100:end,1:1000,:);
    
    t1 = rgb2gray(IN);
    t2 = binning(t1,0,0);
    %t2 = imresize(t1,scale_f);

    %t = uint8(rgb2gray(t1));
    
    T(:,:,k) = t2;
    
end

Tim = mean(T,3);
Tim = uint8(Tim);

imwrite(Tim,strcat(folder_out,'/',video_n,'_media_fps_',num2str(fps),'_hout_',num2str(h_out),'_ext_',ext,'.png'));

%%%
for k = 1:screen_to_load

    [dummy] = process_bar(k,screen_to_load,'FASE 2 - sottraggo media');
       
    res(:,:,k) = (double(T(:,:,k)) - double(Tim));
    
end

%res = CS_fast(res);

clc
disp('-------------------->|  FASE 3 - salvo')
PAR.T = T;
PAR.Tim = Tim;
PAR.res = res;

save(strcat(folder_out,'/',video_n,'_frame_fps_',num2str(fps),'_hout_',num2str(h_out),'_ext_',ext,'.mat'),'PAR');

clc
disp('-------------------->|  FASE 4 - DONE!')
clear all