%% lavoro su immagine e salvo dati
clc
clear
close all

load_str = 'pulitav3.bmp';
TL = [240 230];
BR = [320 585];

th_lim = 43;

ff_x = [0.28,0.66,1;
        14 13 13 ;
        12 11 12];

[O,i,N,X,MATairfoil] = Calcolo_strisce(load_str,1.45,50,TL,BR,th_lim,ff_x);

% cambio canali
O_temp = O; O_temp(:,:,1) = O(:,:,2); O_temp(:,:,2) = 2*O(:,:,1);
O = O_temp; clear O_temp

figure(1)
imshow(O);
title('Immagine in seguito a elabolazione')


save('redo.mat')

%%
load('redo.mat');

% corda espressa in pixel
c_px = norm([(MATairfoil(end-1,1) - MATairfoil(2,1)),...
             (MATairfoil(end-1,2) - MATairfoil(2,2))]); %[px]

Origine = [MATairfoil(2,1);MATairfoil(2,2)] + ...
           0.25*[( MATairfoil(end-1,1)-MATairfoil(2,1) );
                 ( MATairfoil(end-1,2)-MATairfoil(2,2) )];
% corda espressa in metri
c_m  = 0.15; %[m]

K_mpx = c_m/c_px;


% punto mal rilevato durante elaborazione, correggo a mano
% deroga
i{15}(12) = 250;


dx = mean(diff(X));
true = max(size(i{1}));  % numero affidabile di fili di fumo
nX = size(i,2);          % numero di stazioni in X

% inizializzo con prima colonna (dove la misura � affidabile)
i_cell = cell(1,max(size(i{1})));
x_cell = i_cell;

for k = 1:max(size(i{1}))
    i_cell{1,k} = i{1}(k);
    x_cell{1,k} = X(1);
end

% catturo punti
for j = 2:nX %x
    
    for k = 1:size(i_cell,2) % fili
        
        j_temp = size(i_cell{1,k},2);
        % punto k-esimo � il punto con indice riga pi� simile
        % a punto (k-1)esimo, nella foto non ho rilevato tutti i 
        % fili di fumo,quindi i vettori i{j} hanno dimensioni diverse
        
        [err,imin] = min(abs(i{j} - i_cell{1,k}(j_temp)));
        
        % deroga 
        if k == 13 % filo di fumo sul dorso del profilo
            if j == 16
                i_cell{1,k} = [i_cell{k},259];
                x_cell{1,k} = [x_cell{k},270];
            end
            
            if j == 17
                i_cell{1,k} = [i_cell{k},257];
                x_cell{1,k} = [x_cell{k},287];
            end
            
            if j == 18
                i_cell{1,k} = [i_cell{k},256];
                x_cell{1,k} = [x_cell{k},304];
            end
            
            if j == 34
                i_cell{1,k} = [i_cell{k},262];
                x_cell{1,k} = [x_cell{k},582];
            end
            
            if j == 41
                i_cell{1,k} = [i_cell{k},264];
                x_cell{1,k} = [x_cell{k},703];
            end
            
            if j == 42
                i_cell{1,k} = [i_cell{k},264];
                x_cell{1,k} = [x_cell{k},720];
            end
        end
        
        if err > 10
            % non aggiorno perch� fuori tolleranza
        else
            i_cell{1,k} = [i_cell{k},i{j}(imin)];
            x_cell{1,k} = [x_cell{k},X(j)];
        end
    end
end

% interpolo per rilvare linee di corrente
so = size(O);
O2 = zeros(so(1:2)); clear so

x_vect = ceil(X(1):X(nX)); % vettore X sul quale interpolo

% valutando grandezze interpolate negli indici
% ricamp, ottengo valori su cordinate X
x_cam  = round(X);
ricamp = round(linspace(1,size(x_vect,2),nX));

i_spline_mat = zeros(size(i_cell,2),nX);

FILT_media = 1/50*ones(1,50);

for w = 1:size(i_cell,2)
    
    % indice riga
    y_vect = round(spline(x_cell{w},i_cell{w},x_vect));
    
    % filtro media ogni linea di corrente per
    y_vect = round(filter2(FILT_media,y_vect,'same'));
    
    i_spline_mat(w,:) = y_vect(ricamp);
    
    % controllo indici linea di corrente
    % dovendo formare un immagine devono essere
    % interi e nei limiti della dimensione dell'immagine
    % (eventuali correzione da implementare...)
    if sum(y_vect <= 0) > 0;
        error('campo di moto: y minore di zero passo %d',w);
    end    
    if sum(ceil(x_vect(k))<= 0)
        error('campo di moto: x minore di zero passo %d',w);
    end
    if sum(isreal(y_vect)-1) ~= 0
        error('campo di moto: y � complesso passo %d',w);
    end
    if sum(isreal(x_vect(k)) -1) ~= 0
        error('campo di moto: x � complesso passo %d',w);
    end
    
    % compongo immagine
    for k = 1:max(size(x_vect))
        O2(y_vect(k)-1:y_vect(k)+1,x_vect(k)) = 255;
    end
    for k = 1:20:max(size(x_vect))
        O2(y_vect(k)-2:y_vect(k)+2,x_vect(k)) = 255;
    end
    
end

O2 = uint8(O2);
O2 = uint8(repmat(O2,1,1,3));
figure; imshow(O2);

return
%indici da scartare (per via del filtraggio)
i_discard = 3;


% continuit� -> rho*A*V = C ( ip: rho costante)
C = diff(mean(i_spline_mat(:,(i_discard)),2))*1.5;

V_mat = zeros(size(i_cell,2)-1,nX);
y_cam = zeros(size(i_cell,2)-1,nX);

for rr = 1:nX
    % matrice velocit�
    V_mat(:,rr) = C./(diff(i_spline_mat(:,rr)));
    % matrice indici riga associati alla velocit� calcolata
    y_cam(:,rr) = i_spline_mat(1:(end-1),rr) + 0.5*diff(i_spline_mat(:,rr));
end

% elimino valori nel profilo 
V_mat(13,14:33) = NaN;

% riassemblo foto campo di moto in una foto rgb
x_px = round([x_cam - 0.5*(x_cam(2)-x_cam(1)),x_cam(end)+(x_cam(2)-x_cam(1))]);

foto_r = zeros(size(O2,1),size(O2,2));
foto_g = foto_r;
foto_b = foto_r;

red_map = 0.7*max(max(V_mat(2:end-1,i_discard:end-i_discard)));

blue_map = min(min(V_mat(2:end-1,i_discard:end-i_discard)));

for k = i_discard:size(x_cam,2)-i_discard
    
    for j = 2:size(i_spline_mat,1)
        
        [out] = r2rgb(V_mat(j-1,k),red_map,blue_map,128);
        
        foto_r(round(i_spline_mat(j-1,k):i_spline_mat(j,k)),x_px(k):x_px(k+1)) = out(1);        
        
        foto_g(round(i_spline_mat(j-1,k):i_spline_mat(j,k)),x_px(k):x_px(k+1)) = out(2);        
                
        foto_b(round(i_spline_mat(j-1,k):i_spline_mat(j,k)),x_px(k):x_px(k+1)) = out(3);        
        
    end
       
end

foto(:,:,1) = foto_r;
foto(:,:,2) = foto_g;
foto(:,:,3) = foto_b;
foto = uint8(foto);
foto_profilo = uint8(repmat(O(:,:,3),1,1,3));

%
figure(2)
imshow(foto+O2+foto_profilo)
title(sprintf('Campo di moto rilevato: Vmax = %2.2f, Vmin = %2.2f [m/s]',...
    max(max(V_mat(:,i_discard:size(x_cam,2)-i_discard))),...
    min(min(V_mat(:,i_discard:size(x_cam,2)-i_discard)))));
return
figure(3)
imshow(O+O2)
title('Confronto')


% misura scia
K_scia = 0.975;
foto_scia = [];

% assemblo fot con valori numerici (no colorband)

for k = 34:size(x_cam,2)-i_discard
    
    for j = 2:size(i_spline_mat,1)
 
        foto_scia(round(i_spline_mat(j-1,k):i_spline_mat(j,k)),x_px(k):x_px(k+1)) = V_mat(j-1,k);        
        
    end
       
end

foto_scia = foto_scia(:,x_px(34):end,:);
h_line_cor = diff(i_spline_mat(:,34:end));
% unica linea consistentemente appartenente alla scia � la #34

% riporto griglia concorde con CFD
x_cam  = (x_cam - Origine(1))*K_mpx; %[mm]
y_cam  = -((y_cam - Origine(2))*K_mpx); %[mm]

airfoil = K_mpx*[(MATairfoil(2:end-1,1) - Origine(1)),...
    -([(MATairfoil(2:end-1,2)),(MATairfoil(2:end-1,3))] - Origine(2)) ];

% CFD
% spessore scia
% carico dati sperimentali
a = load('dati_scia.dat');

[~,c] = min(a(:,1));

y_pc = a(1:c,1);
y = linspace(y_pc(1),y_pc(end),1000);
num_line = 6;

i_med = 350;

for j = 1:num_line
    is = (j-1)*c + 1;
    ie = j*c;
    u(:,j) = spline(y_pc,a(is:ie,2),y);
    
    [min_u(j),mw(j)] = min(u(:,j));
    
    % u ref
    u_ref_d(j) = K_scia*mean(u(1:i_med,j));
    u_ref_v(j) = K_scia*mean(u(end-i_med:end,j));
    
% % V REF diversa ogni stazione
%     u_up  = u(1:mw(j),j);
%     u_dwn = u(mw(j):end,j);
%     [dummy,up]  = min( abs(u_up-(u_ref_d(j))));
%     [dummy,dwn] = min( abs(u_dwn-(u_ref_v(j))));
%     
%     wp(j,:) = [up,(mw(j)+dwn-1)];
      
end

for j = 1:num_line   
    u_up = u(1:mw(j),j);
    u_dwn = u(mw(j):end,j);
    
% V REF uguale per ogni stazione e versante   
%     [dummy,up]  = min( abs(u_up-mean([u_ref_d,u_ref_v])));
%     [dummy,dwn] = min( abs(u_dwn-mean([u_ref_d,u_ref_v])));
    
% V REF diversa ogni stazione ma due versanti
    [~,up]  = min( abs(u_up-mean([u_ref_d])));
    [~,dwn] = min( abs(u_dwn-mean([u_ref_v])));


    wp(j,:) = [up,(mw(j)+dwn-1)];
        
end

x_m_abs = x_cam(34:end);
x_m_abs = x_m_abs(1:3:end);

%x_m = x_m_abs - x_m_abs(1);

clear a b c ie is j dummy

lato = ceil(sqrt(size(u,2)));

fprintf('Spessore medio CFD   = %f [mm]\n',1000*mean(y(wp(:,1))-y(wp(:,2))));
fprintf('Spessore medio video = %f [mm]\n',1000*mean(h_line_cor(13,:)*K_mpx));

t_cfd = mean(y(wp(:,1))-y(wp(:,2)));

% CONFRONTO CFD

cfdata = load('dati_scia.csv');

%indici regione da ricostruire
i_x = i_discard:nX-i_discard;
i_y = [1:size(V_mat,1)];

V_cfd = zeros(max(size(i_y)),max(size(i_x)));
for j = 1:max(size(i_x)) %x
    
    for k = 1:max(size(i_y)) %y
        
        pt_tgt = [x_cam(i_x(j)),y_cam(i_y(k),i_x(j))];
        
        if isnan(V_mat(k,j)) == 1
            V_cfd(k,j) = NaN;
        else
            [cfdata_err,imin] = min(sqrt((cfdata(:,1)-pt_tgt(1)).^2 + (cfdata(:,2)-pt_tgt(2)).^2));
            
            V_cfd(k,j) =  cfdata(imin,4);
        end
    end
end

figure(1000)
surf(repmat(x_cam(i_x),size(y_cam,1),1),y_cam(:,i_x),V_mat(:,i_x))
figure(1001)
surf(repmat(x_cam(i_x),size(y_cam,1),1),y_cam(:,i_x),V_cfd)
figure(1002)
surf(repmat(x_cam(i_x),size(y_cam,1),1),y_cam(:,i_x),(V_cfd-V_mat(:,i_x)))
title('Errore = V_{CFD} - V_{Visualizzazione}')
colorbar
err_moto = (V_cfd-V_mat(:,i_x));
EEE = [];
for w = 1:size(err_moto,1) 
    for j = 1:size(err_moto,2)
        if isnan(err_moto(w,j)) == 0
            EEE = [EEE,err_moto(w,j)]; 
        end;
    end;
end; 
%
figure(10)
for j = 1:num_line
    

    subplot(lato-1,lato,j)
    plot(y,u(:,j),'bo-');
    axis([-0.03 0.03 0 1.7])
    hold on
    grid on
    
    plot(y(wp(j,1)),u(wp(j,1),j),'ro','LineWidth',3)
    plot(y(wp(j,2)),u(wp(j,2),j),'ro','LineWidth',3)
    
% % V REF diversa ogni stazione
%     plot(y(1:mw(j)),u_ref_d(j).*ones(1,mw(j)),'r');
%     plot(y(mw(j):end),u_ref_v(j).*ones(size(y(mw(j):end))),'r');

% % V REF diversa ogni stazione e versante
%     plot(y(1:mw(j)),mean([u_ref_d,u_ref_v]).*ones(1,mw(j)),'r');
%     plot(y(mw(j):end),mean([u_ref_d,u_ref_v]).*ones(size(y(mw(j):end))),'r');

% V REF diversa ogni stazione ma due versanti
    plot(y(1:mw(j)),mean([u_ref_d]).*ones(1,mw(j)),'r');
    plot(y(mw(j):end),mean([u_ref_v]).*ones(size(y(mw(j):end))),'r');
    title(sprintf('x = %2.3f',x_m_abs(j)))

    
end

figure(11)
subplot(2,1,1)
scia_th = foto_scia.*(foto_scia < K_scia*1.5);
imshow(scia_th(240:320,:))
%title('Immagine soglia veloct�')

subplot(2,1,2)
plot(1000*x_cam(34:end),1000*h_line_cor(13,:)*K_mpx,'bo-')
hold on
plot(1000*x_m_abs,1000*(y(wp(:,1))-y(wp(:,2))),'ro-');
grid on
legend('Acquisizione','CFD','Location','best')
%title('Andamento scia')
grid on
xlabel('Coordinata X [mm] (rispetto a 25% corda)')
ylabel('Spessore scia [mm]')


save('campo_di_moto.mat')




