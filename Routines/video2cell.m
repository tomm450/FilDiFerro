addpath ./Routines

folder_in = './IN';
folder_out = './OUT';

addpath(folder_in);

fta = 'MVI_9900.MOV';

[A] = importfile(strcat(folder_in,'/',fta));

A = A(140:end,450:end,:,:);

a = sum(A,4);
k = 124/mean(mean(mean(a)));
ak = k*a;

disp('max3(a)) =')
max(max(max(a)))

disp('min3(a)) =')
max(max(max(a)))

imshow(uint8(ak));