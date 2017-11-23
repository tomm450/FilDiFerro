%% STALL
clc
clear all
close all

%%

disp('ESEGUIRE SEZIONE PER SEZIONE (Ctrl+Enter)...')
disp('L accuratezza del risultato dipende dalla scelta dei fotogrammi...')

addpath ./Routines

folder_in = './IN/STALL/';
ext = '.png';

%img_str = ['a1';'a2';'a3'];
img_str = ['a2'];
n_caso = 1;

% PP e indentificazioni punti (a mano...)
DO_pp = 1;

if DO_pp == 1
    for caso = 1:n_caso
        a = imread(strcat(folder_in,img_str(caso,:),ext));
        a_or = rgb2gray(a);
        a = a_or+200;
        
        a_cs = CS_fast(a);
        %a_he = histeq(a);
        
        kcs = uint8(Kuwahara(a_cs,9));
        %khe = uint8(Kuwahara(a_he,9));
        
        se = strel('diamond',1);
        se = se.Neighborhood;
        
        k = imerode(kcs,se);
        diff = kcs - k;
        
        difft = uint8(255*(diff > 15));
        
        figure; imshow(difft); title(img_str(caso,:));
        
        for k = 1:2:11
            b = uint8(k.*a_or); 
            figure;
            imshow(b); 
            title(sprintf('CASO %d; k = %s',caso,num2str(k)));
        end
    end
end

%% MODIFICA (A MANO) MATRICE IN BASE A PLOT
% CASO n -> (x,y)
% 1  LE = 108 698  - TE = 1339 895
% 2  LE = 110 687  - TE = 1338 907
% 3  LE = 112 682  - TE = 1336 925

%    [x1   x2   x3 ;  y1  y2  y3 ];
LE = [108  110  112;  698 687 682];
TE = [1339 1338 1336; 895 907 925];

d = TE - LE;
c = sqrt(d(1,:).^2 + d(2,:).^2);

disp('angoli rilevati [deg]')
a = 180/pi.*(atan(d(2,:)./d(1,:)))


figure(7); 
hold on; 
plot([110 1338],[687 907],'o-','LineWidth',4)
plot([110],[687],'rx-','LineWidth',6);
plot(1338,907,'gx-','LineWidth',6);
legend('Corda',...
    sprintf('Leading Edge = [%d %d]',[110],[687]),...
    sprintf('Trailing Edge = [%d %d]',1338,907));
title(sprintf('Visualizzazione ad alpha = %2.2f [deg]',a(2)))

fprintf('SALVA A MANO');


