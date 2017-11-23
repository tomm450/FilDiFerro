function [O,index,N,X,mat] = Calcolo_strisce(load_str,K,n,TL,BR,th_lim,ff_x)
clc
close all

addpath ../Routines/
addpath ../IN/

PLT_pre = 0;


a = imread(load_str);

if size(a,3) == 3
    a = rgb2gray(a);
end
%figure; imshow(k5); title('KW(5)')
% figure(1)
% imshow(a)
% title('Immagine di Partenza')

% D_y
d_y = 1/3*[ -1 -1 -1; 0 0 0; 1 1 1];

[b_y1] = filt2plot(a, d_y         ,PLT_pre, '(dy_1)');
[b_y2] = filt2plot(a, flipud(d_y), PLT_pre, '(dy_2)');

[b_x1] = filt2plot(a, d_y'         ,PLT_pre, '(dx_1)');
[b_x2] = filt2plot(a, flipud(d_y)' ,PLT_pre, '(dy_1)');


b_y = uint8(sqrt( double(b_y1).^2 + double(b_y2).^2 ));
b_x = uint8(sqrt( double(b_x1).^2 + double(b_x2).^2 ));



%b_ycs = CS_fast(b_y,1);
[b_ycs] = CS_fast(b_y, 0); %, speed, b_l, b_h)
[b_xcs] = CS_fast(b_x, 0);
%b_yhe = histeq(b_y);


if PLT_pre == 1
    figure
    imshow(b_y)
    title('module dy')
    figure
    imshow(b_x)
    title('module dx')
    figure
    imshow(b_ycs)
    title('module dy CS')
    figure
    imshow(b_xcs)
    title('module dx CS')
end

b_xy = uint8(sqrt( double(b_xcs).^2 + double(b_ycs).^2 ));
figure
imshow(b_xy)
title('grad')

% figure
% b_yhe = histeq(b_y);
% imshow(b_yhe)
% title('module dy HE')

%% IDENTIFICO PROFILO
% coordinate punto Top-Left e Bottom-Right rettangolo che contiene profilo
%TL = [240 230];
%BR = [320 585];

% limite soglia tc profilo chiuso
%th_lim = 43;

[mat,~,img_out] = profile_id(CS_fast(b_xy,0),TL,BR,1,th_lim);


a = a-img_out;

imshow(img_out)

% imwrite(img_out,'airfoil_cm.bmp');
% imwrite(a,'maschera.bmp');
%% KUWUHARA
k5 = Kuwahara(a,5); k5 = uint8(k5);
%figure; imshow(k5); title('Kw(5)')
%% CONTRAST STRETCH
[k_cs] = CS_fast(k5,0);
%figure; imshow(k_cs); title('Kw(5) + CS')
% imwrite(k_cs,'starting_point_CS.png')

%% INIZIO
%b_f = a;      % scelgo da quale immagine partire
%b_f = k_he;
b_f = k_cs;

% divisione orizzontale
dy = round(size(b_f,1)/2);
divisX = floor(linspace(1,size(a,2),n+1));

%lumi = zeros(size(b_f,1),n);
ev = zeros(size(b_f)); 
em = ev;
pt = ev;
%fili = [];
N1 = []; N2 = [];


for j = 1: n
    
    % filo di fumo impatta il profilo e si divide
    % a seconda della "fetta di immagine dovr� cercare
    % un numero differente di fili
    
    if j <= round(ff_x(1,1)*n)
        fup  = ff_x(2,1);
        fdwn = ff_x(3,1);
    elseif j>round(ff_x(1,1)*n) && j <= round(ff_x(1,2)*n)
        fup  = ff_x(2,2);
        fdwn = ff_x(3,2);
    else
        fup  = ff_x(2,3);
        fdwn = ff_x(3,3);
    end
    
    [~] = process_bar(j,n,''); 
    
    temp = b_f(:,divisX(j):divisX(j+1),:);
    left = num2str(divisX(j)); right = num2str(divisX(j+1));
    X(j) = mean([divisX(j),divisX(j+1)]);
    
    % DIVISIONE VERTICALE sopra/sotto profilo 
    if divisX(j) >= mat(1,1) && divisX(j+1) <= mat(end,1)
        
        % sono in zona profilo, divido immagine
        [~,imin1] = min(abs(mat(:,1) - divisX(j)));
        [~,imin2] = min(abs(mat(:,1) - divisX(j+1)));
        
        
        d1 = floor(min([mat(imin1,2),mat(imin2,2)]));         
        u2 = ceil(max([mat(imin1,3),mat(imin2,3)]));
        
    else
        
        d1 = dy;
        u2 = dy +1;
        
    end
    % divido striscia in parte alta e bassa
    temp1 = temp(1:d1,:,:);
    temp2 = temp(u2:end,:,:);
    plot_res = 0;
  
%     if j == 5; plot_res = 1; end
%     
%     if j == 25; plot_res = 1; end
%     
%     if j == 45; plot_res = 1; end  

    % parte alta
    % processo secondo input
    [ou1,outM1] = process1(temp1,K,plot_res,j);
    [ N1(j)] = contatore( mean(outM1{1},2) );
    
    k_vect = [1:0.05:2]; kk = 1;
    
    
    if N1(j) ~= fup % non ho trovato numero corretto fili; 
                    % itero soglia fino a trovarli
        while N1(j) ~= fup
            [~,outM1] = process1(temp1,k_vect(kk),plot_res);
            [ Np(kk)] = contatore( mean(outM1{1},2) );
            kk = kk +1;
            
            if kk == 21
                break
            end
        end
        [~,imin] = min(abs(Np - fup));
        
        [ou1,outM1] = process1(temp1,k_vect(imin),plot_res);
        [ N1(j)] = contatore( mean(outM1{1},2) );
    end

    % parte bassa
    % processo secondo input
    [ou2,outM2] = process1(temp2,K,plot_res,j);
    [ N2(j)] = contatore( mean(outM2{1},2) );
    k_vect = [1:0.05:2]; kk = 1;
    
    
    if N2(j) ~= fdwn % non ho trovato numero corretto fili; 
                     % itero soglia fino a trovarli
                     
        while N1(j) ~= fdwn
            [~,outM2] = process1(temp2,k_vect(kk),plot_res);
            [ Np(kk)] = contatore( mean(outM2{1},2) );
            kk = kk +1;
            
            if kk == 21
                break
            end
        end
        [~,imin] = min(abs(Np - fdwn));
        
        [ou2,outM2] = process1(temp2,k_vect(imin),plot_res);
        [ N2(j)] = contatore( mean(outM2{1},2) );
    end

    
    
    % ricompongo immagine da soglie locali    
    ev(1:d1,divisX(j):divisX(j+1))   = ou1{1};
    ev(u2:end,divisX(j):divisX(j+1)) = ou2{1};
    
    em(1:d1,divisX(j):divisX(j+1))   = outM1{1};
    em(u2:end,divisX(j):divisX(j+1)) = outM2{1};
    
     
    pt(1:d1,  round(0.5*(divisX(j)+divisX(j+1)))) = mean(outM1{1},2);
    pt(u2:end,round(0.5*(divisX(j)+divisX(j+1)))) = mean(outM2{1},2);
    
    N(j) = N1(j) + N2(j);
    index{j} = [outM1{3},(outM2{3}+u2-1)];
  
end

% figure; imshow(ev)
% 
% figure; imshow(em)
O = uint8(zeros(size(em,1),size(em,2),3));
O(:,:,1) = b_f;
O(:,:,2) = pt;
O(:,:,3) = img_out;
