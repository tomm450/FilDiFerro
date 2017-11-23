function [dritto] = correction_botte(f,cup_export,cdwn_export,point,ref)
close all
cup = cup_export.*norm(ref);
cdwn = cdwn_export.*norm(ref);

cup  = round( [(cup(:,1) + point(1)),(cup(:,2) + point(2))] );
cdwn = round( [(cdwn(:,1) + point(1)),(cdwn(:,2) + point(2))] );

X = size(f,2);
Y = size(f,1);

if sum(cup(:,1)<1) > 0
    % devo partire dopo
    [err,imin] = min(abs(cup(:,1)-1));
    cup  = cup(imin:end,:);
    cdwn = cdwn(imin:end,:);
end
if sum(cup(:,1)>X) > 0
    % devo troncare
    [err,imin] = min(abs(cup(:,1)-X));
    cup  = cup(1:imin-1,:);
    cdwn = cdwn(1:imin-1,:);
end

if sum(cup(:,2)<1) > 0
    % devo partire dopo
    delta = 1+abs(min(cup(:,2)));
    f = [zeros(delta,size(f,2));f];
    cup(:,2) = cup(:,2) + delta;
    cdwn(:,2) = cdwn(:,2) + delta;
end
if sum(cup(:,2)>Y) > 0
    % devo troncare
    error
end

if sum(cdwn(:,2)<1) > 0
   % devo partire dopo
   error
end
if sum(cdwn(:,2)>Y) > 0
    % devo troncare
    
    delta = 1+abs(-Y+max(cdwn(:,2)));
    f = [f;zeros(delta,size(f,2))];

end


q = zeros(size(f));

for k = 1:size(cup,1)
    q(cup(k,2),cup(k,1)) = 255;
end

for k = 1:size(cdwn,1)
    q(cdwn(k,2),cdwn(k,1)) = 255;
end


rgb(:,:,2) = f;
rgb(:,:,1) = uint8(q);
rgb(:,:,3) = uint8(q);

% figure; imshow(rgb)
% title('botti')

h_vect = (cdwn(:,2) - cup(:,2));
h_new = min(h_vect);
h_old = max(h_vect);
l = cup(end,1)-cup(1,1);

K = h_new/h_old;

dritto = [];
%sporco = zeros(size(f));
y = [1:(h_new)];

w = 1;
x_old = 0;
for k = 1:size(cup,1)
    if cup(k,1) == x_old
        % do nothing
    else
    temp = double(f(cup(k,2):cdwn(k,2),cup(k,1)));
%    sporco(cup(k,2):cdwn(k,2),cup(k,1)) = temp;
    y_temp = 1:size(temp,1);
    y_yemp = (y_temp./max(y_temp)).*h_new;
    new = spline(y_temp,temp,y');
    
    dritto(:,w) = (new);
    w = w+1;
    x_old = cup(k,1);
    end
    
end
figure
imshow(histeq(uint8(dritto)));
title('dritto')


dritto_prop = imresize(uint8(dritto),[size(dritto,1),round(h_new/h_old*size(dritto,2))]);
figure
imshow(dritto_prop)
title('proporzioni originali')

end