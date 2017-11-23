function [mat,MAT,img_out] = profile_id(a,TL,BR,fill,th_lim)
PLT = 0;
if nargin == 3
    fill = 0;
end

%% routines per ricerca soglia di th
if nargin < 5 
happy = 0;
    
th_vect = floor([mean(mean(double(a))):10:max(max(double(a)))]);
th_in = a;

while happy == 0
close all
for w = 1:size(th_vect,2)
   th_img = 255*(th_in > th_vect(w));
   th_img = uint8(th_img);

   figure; imshow(th_img); title(sprintf('Soglia a %d',th_vect(w)));

end

happy = input('soglia trovata? [1 = sì, 0 = no]\n');
if happy == 0 
    again = input('continuare ricerca? [1 = sì, 0 = no]\n');
    
    if again == 0
        return
    else
        alpha = input('inserire nuova soglia inferiore [numero intero]');
        omega = input('inserire nuova soglia superiore [numero intero]');
        delta = input('inserire nuovo delta [numero intero tc th_vect = alpha:delta:omega]');
        th_vect = floor([alpha:delta:omega]);  
    end
elseif happy == 1
    th_lim = input('soglia identificata a...');
end

end

end

th = uint8(255*( a > th_lim));
%figure; imshow(th); title('th per identificare profilo')
%%

%%% PARTE MODIFICA
thtemp = th(TL(1):BR(1),TL(2):BR(2),:);
th_line = sum(thtemp,1);
% figure; imshow(thtemp); title('th  secondo TL e BR')
i = 0; j = 1;
while i == 0
    i = th_line(j);
    j = j+1;
end;
jl = j;

i = 0; j = size(th_line,2);
while i == 0
    i = th_line(j);
    j = j-1;
end;
jr = j;

% fprintf('Profilo identificato: %d < Airfoil < %d \n',jl-1,jr+1);

th2 = thtemp(:,jl:jr,:);
%figure;imshow(th2); title('th  dove ho binaco')
camp = floor(linspace(1,size(th2,2),20));
%
jup = []; jdwn = [];

% 
for w = 1: size(camp,2)
    
    i = 0; j = 1;
    while i == 0
%         disp('-------------------')
%         camp(w)
        i = th2(j,camp(w)); % riga j che scorre, su colonna campione
        j = j+1;
    end;
%     jup(w) = j-1;
    jup(w) = j-2;
    i = 0; j = size(th2,1);
    while i == 0
        i = th2(j,camp(w));
        j = j-1;
    end;
%     jdwn(w) = j+1;
    jdwn(w) = j+2;
end

jup = jup + TL(1) ; jdwn = jdwn + TL(1);

x_cam = camp + TL(2) + jl -2;

MAT(:,1) = x_cam';
MAT(:,2) = jup';
MAT(:,3) = jdwn';

mat(:,1) = [MAT(1,1):MAT(end,1)]';
mat(:,2) = spline(MAT(:,1),MAT(:,2),mat(:,1));
mat(:,3) = spline(MAT(:,1),MAT(:,3),mat(:,1));

img_out = zeros(size(a));


% chiusura TE e LE
if MAT(1,2) == MAT(1,3) % bordo chiuso
else 
    MAT = [MAT(1,1)-1,0.5*(MAT(1,2) + MAT(1,3)),0.5*(MAT(1,2) + MAT(1,3));MAT];
end

if MAT(end,2) == MAT(end,3) % bordo chiuso
else 
    MAT = [MAT;MAT(end,1)+1,0.5*(MAT(end,2) + MAT(end,3)),0.5*(MAT(end,2) + MAT(end,3))];
end    

for j = 1:size(mat,1)
     q = mat(j,1);
    u = ceil(mat(j,2));
    d = floor(mat(j,3)-2);
%     u = floor(mat(j,2));
%     d = ceil(mat(j,3));   
    if fill == 0
        img_out(u,q) = 255;
        img_out(d,q) = 255;
    else
        img_out(u:d,q) = 255;
    end
end


  
img_out = uint8(img_out);

% divione
MAT1 = [1,MAT(1,2),MAT(1,2)];
MATEND = [size(a,2),MAT(end,2),MAT(end,2)];

MAT = [MAT1; MAT; MATEND];


mat = [MAT(1:2,:); mat; MAT(end-1:end,:)];


% disp('img_out e opzione fill non sono ancora stati implementati')


if PLT == 1
figure
plot(MAT(:,1),MAT(:,2),'o')
hold on
plot(MAT(:,1),MAT(:,3),'o')
axis equal
grid on
title(strcat('Contorno profilo rilevato',num2str(th_lim)))

hold on
plot(mat(:,1),mat(:,2),'-')
hold on
plot(mat(:,1),mat(:,3),'-')
axis equal
grid on
end



end

