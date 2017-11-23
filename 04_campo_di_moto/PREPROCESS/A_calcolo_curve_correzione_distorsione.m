close all
clc
clear 

% CALCOLO CURVE CORREZIONE DISTORSIONE A BOTTE
% carico fot galleria spenta
f = imread('foto.JPG');
f = rgb2gray(f);

% caratteristiche f
LE_f = [2303 2029]; % leading edge
TE_f = [3753 1981]; % trailing edge
d_f = TE_f-LE_f;    
a_f = atan2(d_f(2),d_f(1));  % incidenza
m_f = [3112 1185];           % coordinate punto luminoso preso come riferimento
p_f = round(LE_f + 0.3*d_f); % coordinate asse rotazione
ref_f = m_f-p_f;             % lunghezza di riferimento



% in f si nota la distorsione dell'immagine, un filtraggio in y rende
% ancora meglio visibile 
d_y = [1 1 1; 0 0 0 ; -1 -1 -1];
k1 = filt2plot(f,[1 1 1; 0 0 0; -1 -1 -1],0,'');


% campiono punti su curve
% stesso inizio e stessa fine in X
pt_up = [1195 610; 2233 577; 3280 562; 4690 575; 5184 577];
cup = round(spline(pt_up(:,1),pt_up(:,2),[pt_up(1,1):pt_up(end,1)]'));
cup = [ [pt_up(1,1):pt_up(end,1)]', cup];

pt_dwn = [1195 3413; 1815 3427; 2883 3433; 4292 3403; 5184 3384];
cdwn = round(spline(pt_dwn(:,1),pt_dwn(:,2),[pt_dwn(1,1):pt_dwn(end,1)]'));
cdwn = [ [pt_dwn(1,1):pt_dwn(end,1)]', cdwn];

% normalizzo curve in modo da poter correggere anche immagine di diversa
% dimensione
cdwn_export = [(cdwn(:,1) - m_f(1)),(cdwn(:,2) - m_f(2))]; cdwn_export = cdwn_export./norm(ref_f);

cup_export = [(cup(:,1) - m_f(1)),(cup(:,2) - m_f(2))]; cup_export = cup_export./norm(ref_f);

save('botte.mat','cdwn_export','cup_export');

figure
imshow(histeq(k1))
hold on
plot(cup(:,1),cup(:,2),'g','LineWidth',4)
plot(cdwn(:,1),cdwn(:,2),'g','LineWidth',4)
title('Curve distorsione obiettivo')
