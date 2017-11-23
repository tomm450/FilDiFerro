function [ v_he ] = visteq( v_in )
% HIST EQ PER VIDEO

x = size(v_in,2);
y = size(v_in,1);
t = size(v_in,3);

v_he = zeros(size(v_in));
temp = zeros(y,t*x);
for k = 1:t

temp(:,(k-1)*x+1:k*x) = v_in(:,:,k);

end

temp = histeq(uint8(temp));

for k = 1:t

v_he(:,:,k) = temp(:,(k-1)*x+1:k*x); 


end

v_he = uint8(v_he);