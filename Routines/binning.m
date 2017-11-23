function [a2] = binning( a,M,r )

% Binning
%M = 0 binning puro
%M = 1 mediato (per riduzione risoluzione)

% r = 0 faccio una sola volta
% r = n faccio n+1 volte in maniera ricorsiva
if nargin == 1
    M = 0;
    r = 0;
elseif nargin == 2
    r = 0;
end

z = size(a,3);

if M == 0
    C = 1;
elseif M == 1
    C = 0.25;
else
    error('Controllare input M')
end
a2 = a;
% Binning
for w = 1:r+1
    
    a_old = a2;
    clear a2
    x = size(a_old,2);
    y = size(a_old,1);
    
    for k = 1:z
        a2(:,:,k) = C*a_old(1:2:y-1,1:2:x-1,k)...
                   +C*a_old(2:2:y,1:2:x-1,k)...
                   +C*a_old(1:2:y-1,2:2:x,k)...
                   +C*a_old(2:2:y,2:2:x,k);
    end
end
