clc
clear 
close all
%% PREPROCESS FOTO CAMPO DI MOTO

% importo file video avvio galleria
load('MVI_9904_frame_fps_24_hout_640_ext_png.mat');
% importo curve di correzione normalizzate 
load('botte.mat');


a = PAR.T;
% media galleria spenta
b = uint8( mean(a(:,:,1:100),3));

% caratteristiche f
LE = [155 353]; % leading edge
TE = [465 464]; % trailing edge

m = [346 196];  % coordinate punto luminoso preso come riferimento
lumi_ref = 111; % luminosit� 

d = (TE - LE); 
l = norm(d);   % lunghezza di riferimento
ref = round(LE+0.3*(TE - LE)); % coordinate asse rotazione
                

figure; 
imshow(histeq(uint8(b))); title('Alone')


core = double(b);
core(325:490,130:510) = zeros(size(core(325:490,130:510)));
figure; 
imshow(histeq(uint8(core))); title('histeq(Patch_alone)')

% normalizzo patch
core = 1.1*core./lumi_ref;

% importo foto da correggere
c = imread('video.bmp');

% caratteristich� foto
LE_v = [313 273];
TE_v = [589 262];
d_v = TE_v - LE_v;
l_v = norm(d_v);

a_v = atan2(d_v(2),d_v(1));
m_v = [467 111];
p_v = round(LE_v + 0.3*d_v);

ref_v = m_v - p_v;

%% CORREZIONE ALONE
% riscalo 
d = rgb2gray(imresize(c,(l/l_v)));

% figure
% imshow(d)

% il punto di riferimento ha cambiato coordinate
m_v = [557 132];
lumi_ref_v = 171;


d2    = [zeros(64,size(d,2));d];
core2 = [zeros(size(core,1),211),core];

d3 = [d2,zeros(size(d2,1),325)];
core3 = [core2;zeros(6,size(core2,2))];

% figure;imshow(d3)
% figure;imshow(core3)

% Sottraggo alone
mega = double(d3) - double(core3)*lumi_ref_v;
figure;
imshow(histeq(uint8(mega)));
% crop 

return
mega2 = uint8(mega(100:end,1:1000));

% ricampio punti notevoli e correggo botte
m_v_mega = [557 97];
LE = [373 289];
TE = [703 277];
d = TE-LE;
p = round(LE + 0.3*d);

ref_v = m_v_mega - p;
[dritto] = correction_botte(mega2,cup_export,cdwn_export,m_v_mega,ref_v);

best = uint8(dritto(50:580,10:end));
figure; imshow(best);title('Risultato finale')
imwrite(best,'pulita.bmp'); % <- immagine di partenza per calcolare 
                              %    il campo di moto
