%% CREO VIDEO DA FIGURE

nome_out = 'Shedding24fpsMask(12pfs).avi'; % compreso percorso


dir2read = '.';          % percorso figure da assemblare
str_p = 'figure';        % prefisso figure da caricare
str_s = '_fps_50.png';   % suffisso
frame2read = 120;        % numero figure da leggere

framerate_out = 12;
%% fine input

images = cell(frame2read,1);

for k = 1:frame2read % assumo che inizio a leggere sempre dal primo fotogramma
    images{k} = imread(strcat(dir2read,'/',str_p,num2str(k),str_s));
end

 writerObj = VideoWriter(nome_out);
 writerObj.FrameRate = framerate_out;  
 
 open(writerObj);
 
 for u=1:length(images)
     
     % convert the image to a frame
     frame = im2frame(images{u});
     writeVideo(writerObj, frame);
 end

 % close the writer object
 close(writerObj);