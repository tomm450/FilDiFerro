clc
close all
clear all

%%% DIMENSIONI RICIRCOLI
% file a 24fps
load('MVI_9906_frame_fps_24_hout_640_ext_png.mat');

% fotogramma di interesse
n = 39; % fotogramma migliore


% dopo aver plottato il frame di interesse so che 
LE = [102 340];
TE = [379 508];


d_x = TE(1) -LE(1); 
d_y = TE(2) -LE(2);
alpha = 180/pi*(atan(d_y/d_x));
c = 0.15; %[m]

% crop per plot
l = 460; r = 815; 
u = 390; d = 580;
crop =  (uint8(PAR.T(u:d,l:r,n))); % FRAME SPECIFICO
K = c/sqrt(d_x^2 + d_y^2); 


%% post process
med = filt2plot(CS_fast(crop),1/9*ones(3),0);
kw5 = Kuwahara((med),5);

z = kw5;
zz = (z>1.35*(mean(mean(z))));
zz2 = [zeros(u,size(zz,2));zz;zeros(size(PAR.T,1)-d-1,size(zz,2))];
zz2 = [zeros(size(zz2,1),l),zz2,zeros(size(zz2,1),size(PAR.T,2)-r-1)];

rgb(:,:,1) = uint8(zz2*255);
rgb(:,:,2) = CS_fast(PAR.Tim);
rgb(:,:,3) = zeros(size(rgb,1),size(rgb,2));




% punti ricircolo

% a mano inserisco i punti di ricircolo
pt = [461 772 705; % X
      508 456 406];
fprintf('Ricircoli:\n')
Dx = K*(max(pt(1,:)) - min(pt(1,:)))
Dy = K*(max(pt(2,:)) - min(pt(2,:)))

% plot
figure;
imshow(uint8(rgb));
hold on 
plot(pt(1,1),pt(2,1),'bo','LineWidth',3)
plot(pt(1,2),pt(2,2),'co','LineWidth',3)
plot(pt(1,3),pt(2,3),'go','LineWidth',3)

legend(sprintf('[%d %d]',pt(1,1),pt(2,1)),sprintf('[%d %d]',pt(1,2),pt(2,2)),sprintf('[%d %d]',pt(1,3),pt(2,3)))
% plot([461,772],[590 590],'w','LineWidth',3)
% plot([820,820],[406 508],'w','LineWidth',3)
title(sprintf('Dimensioni ricircoli: Dx = %1.4f ,Dy = %1.4f [m] ',Dx,Dy))


%%
n = 12; % fotogramma migliore�

% crop
l = 675; r = 837; 
u = 386; d = 490;
crop =  (uint8(PAR.T(u:d,l:r,n))); % FRAME SPECIFICO
K = c/sqrt(d_x^2 + d_y^2); 


% pp
med = filt2plot(CS_fast(crop),1/9*ones(3),0);
kw5 = Kuwahara((med),5);


z = kw5;
zz = (z>1.35*(mean(mean(z))));
zz2 = [zeros(u,size(zz,2));zz;zeros(size(PAR.T,1)-d-1,size(zz,2))];
zz2 = [zeros(size(zz2,1),l),zz2,zeros(size(zz2,1),size(PAR.T,2)-r-1)];

rgb(:,:,1) = uint8(zz2*255);
rgb(:,:,2) = CS_fast(PAR.Tim);
rgb(:,:,3) = zeros(size(rgb,1),size(rgb,2));
clear pt


a = (sum(zz2,1) > 0); [err,pt(1,1)] = max(a);

a2 = a; a2(1:pt(1,1)+1) = 1; [err,pt(1,2)] = min(a2);

a = (sum(zz2,2) > 0); [err,pt(2,1)] = max(a);

a2 = a; a2(1:pt(2,1)+1) = 1; [err,pt(2,2)] = min(a2);

fprintf('Vortice:\n')
Dxv = K*(max(pt(1,:)) - min(pt(1,:)))
Dyv = K*(max(pt(2,:)) - min(pt(2,:)))

dv = mean([Dxv,Dyv])

toll = 0.5*K
% plot
figure;
imshow(uint8(rgb));
hold on 
plot(pt(1,1),pt(2,1),'c+','LineWidth',2)
plot(pt(1,2),pt(2,2),'g+','LineWidth',2)

plot(mean(pt(1,:)),mean(pt(2,:)),'go','LineWidth',2)

legend(sprintf('[%d %d]',pt(1,1),pt(2,1)),sprintf('[%d %d]',pt(1,2),pt(2,2)),...
    sprintf('Asse vortice [%d %d]',mean(pt(1,:)),mean(pt(2,:))));
title(sprintf('Dimensioni vortice: d_{medio} = %1.4f [m]',dv))


fprintf('\nSALVA A MANO\n');