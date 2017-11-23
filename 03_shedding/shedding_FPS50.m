%%
clc
clear 
close all

% Variabili Video
addpath ../Routines
addpath ../MAT

load('MVI_9911_frame_fps_50_hout_640_ext_png.mat');

% scelta sequenza da processare
seq = PAR.T;
%seq = PAR.res;

fps = 50;

% incidenza 

LE = [102 340];
TE = [379 508];

d_x = TE(1) -LE(1); 
d_y = TE(2) -LE(2);
alpha = 180/pi*(atan(d_y/d_x));

c = 0.15; %[m]
c_maestra = abs(c*sin(alpha*pi/180));
v_media = 1.5; %[ms-1]

fprintf('Acquisizione a %d FPS, alpha = %2.2f deg\n',fps,alpha)

k_start = 1;
k_end   = size(seq,3);

FLAG_cycle = 1;
  FLAG_visible = 'off';
  FLAG_save = 1;
FLAG_lumi = 0;

KC = 0;  % luminosit� campionamento
KS = 0;  % luminosit� subcampionamento

t = 1/fps*[1:k_end-k_start+1];
freq = (fps/max(size(t)))*(0:(max(size(t))-1));

toll_freq    = (freq(2)-freq(1))/2;
toll_sthrual = (toll_freq)*c_maestra/v_media;

figure(1); 
imshow(histeq(PAR.Tim));
hold on
plot(LE(1),LE(2),'ro','LineWidth',7)
plot(TE(1),TE(2),'ro','LineWidth',7)
title(sprintf('histeq(Media acquisizione); alpha = %2.2f deg\n',alpha))
%% PARAMETRI FINESTRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crop #1 (per ragioni grafiche)
l = 415; r = 715; 
u = 200; d = 580;
% crop #2 (per analisi relativi a crop #1)
U = 220; D = 320; 
L = 120; R = 240; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = round(0.5*(l+r));    
M = round(0.5*(L+R)-L);

% Crop ( e superCrop)
vis   = uint8(KC + seq(u:d,l:r,k_start:k_end));
crop = uint8(KS + vis(U:D,L:R,:));

% Media nella terza dimensione
mean_vis  = mean(vis,3);
mean_crop = mean(crop,3);

% ARRAY 3d formato dalla media
mean_vis   = repmat(mean_vis,1,1,size(vis,3));
mean_crop  = repmat(mean_crop,1,1,size(vis,3));

% POST PROCESS GENERALE non posso fare modifiche al singolo fotogramma
% in un analisi complessiva 
seqPP  = visteq(seq);      
visPP  = visteq(vis);
cropPP = CS_fast(crop);

% pp togliendo media
vis_ppo  = CS_fast(vis-uint8(mean_vis));
crio_ppo = CS_fast(crop-uint8(mean_crop));

[mask_vis]  = mirino([u-1,l-1],[d+1,r+1],size(seq(:,:,1)),4);
[mask_crop] = mirino([U-1,L-1],[D+1,R+1],size(vis(:,:,1)),4);

%% ESTRAGGO VETTORI LUMINOSITA'
% non normalizzati %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% valore medio casella di crop
L_crop = permute(mean(mean(crop,2),1),[3 1 2]);
% valore medio casella di crop + PP
L_crop_pp = permute(mean(mean(cropPP,2),1),[3 1 2]);
% valore medio casella di (crop - media_crop )+ PP
L_crop_ppo = permute(mean(mean(crio_ppo,2),1),[3 1 2]);

% valori colonna mezzeria a crop + PP
LINEA = permute(mean(cropPP(:,M,:),1),[3 1 2]);

%% normalizzo e pulisco da media segnale

L_crop_n = (L_crop-mean(L_crop))./max(L_crop-mean(L_crop));
L_crop_ppn = (L_crop_pp-mean(L_crop_pp))./max(L_crop_pp-mean(L_crop_pp));
L_crop_ppno = (L_crop_ppo-mean(L_crop_ppo))./max(L_crop_ppo-mean(L_crop_ppo));

LINEAn = (LINEA-mean(LINEA))./max(LINEA-mean(LINEA));


%% Analisi Segnali

nfig = 10000;
spl  = 0;
[s1] = Ps(L_crop_n,freq,nfig,spl);

[s2] = Ps(L_crop_ppn,freq,nfig,spl);

[s3] = Ps(L_crop_ppno,freq,nfig,spl);

[s4] = Ps(LINEAn,freq,nfig,spl);

S_mat = [s1,s2,s3,s4];

FFF = figure(nfig); 
if spl == 1
    subplot(1,2,1);
end
legend('Crop','Crop+PP','(Crop-Offset)+PP','Linea','Location','best')
if spl == 1
    subplot(1,2,2);
    legend('Crop','Crop+PP','(Crop-Offset)+PP','Linea','Location','best')
end
set(FFF,'Position',[50 50 800 300])

C = [(max(s1)-mean(s1)),(max(s2)-mean(s2)),(max(s3)-mean(s3)),(max(s4)-mean(s4))];              
[dummy,imax] = max(C);

[Fmax,iFREQ] = max(S_mat(:,imax));
FREQ = freq(iFREQ);

sthrual = 0.2;                            % valore teorici
sthrual_calc = FREQ*c/v_media;            % valore calcolato per corpo aerodinamico
sthrual_maestro = FREQ*c_maestra/v_media; % valore calcolato per sezione maestra 


%% PLOT SEQUENZA
if FLAG_cycle == 1
    for k = k_start:k_end
        
        [dummy] = process_bar(k,k_end-k_start+1,'PLOT');
        
        if strcmp(FLAG_visible,'on') == 1
            pause(0.25)
        end
        
         % impostazione posizione,dimensione e visibilit� 
        hFig = figure;
        set(hFig, 'Position', [0 0 1280 800],'visible',FLAG_visible)
                
        subplot(2,3,1);
        imshow(uint8(seqPP(:,:,k))+mask_vis);
        title(strcat('t = ',num2str(k*1/fps),'[s] Seq + PP'));
        
        subplot(2,3,2);
        imshow(uint8(0.8*visPP(:,:,k))+mask_crop);
        
        subplot(2,3,3);
        imshow(uint8(CS_fast(crio_ppo(:,:,k))))
        title('Crop + PP')
        
        subplot(2,3,[4 5 6])
        plot(t(1:k),L_crop_ppno(1:k),'b','LineWIdth',2)
        hold on
        plot(t(k),L_crop_ppno(k),'ro','LineWIdth',2)
        grid on
        axis([(t(1)-0.01) (t(end)+0.01) -1.05 1.05])
        xlabel('Tempo [s]')
        ylabel('Luminosità segnale')
        
        if FLAG_save == 1
            saveas(hFig,...
            strcat('./50fps/figure',num2str(k),'_fps_',num2str(fps),'.png'));
        else
            if k == 30
                pause(2)
                disp('...pulisco plot...')
                close all
                % facciamo riposare la scheda grafica..
            end
        end
    end
end

%% PLOT LUMINOSITA'
if FLAG_lumi == 1
  
    % PLOT andamenti luminosit�
    figure; plot(t,L_crop); hold on; plot(t,L_crop_pp); plot(t,L_crop_ppo) ; title('Crop');
    legend('sandard','PP','PP-offset')
    ylabel('Luminosit�')
    xlabel('Tempo acquisizione')
    
    figure; plot(t,LINEA); title('Linea')
    ylabel('Luminosit�')
    xlabel('Tempo acquisizione')
    
    figure;
    plot(t,L_crop_n)
    plot(t,L_crop_ppn)
    plot(t,L_crop_ppno)
    plot(t,LINEAn)
    title('Grandezze normalizzate');
    grid on
    legend('Crop','Crop+PP','(Crop-Offset)+PP','Linea')
    ylabel('Luminosit�')
    xlabel('Tempo acquisizione')
    

      
end
%%
fprintf('Frequenza registrata .. = %6.3f +- %2.2f [Hz] \n',FREQ,toll_freq);
fprintf('Sthrual teorico ....... = %6.3f \n',sthrual);
fprintf('Sthrual sez maestra ... = %6.3f +- %2.3f\n',sthrual_maestro,toll_sthrual);









