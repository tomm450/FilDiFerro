function [ n,s_v ] = contatore( vect )

% il vettore è binario ( ==0 o ~= 0) voglio contare il numero di bande e il
% loro spessore

s_v = [];
v = 1;
s = 0;
old = 0;

for j = 1: size(vect,1)
    
    if vect(j) == old % mantengo 
        if old == 0 
            % ero nel vuoto e rimango nel vuoto
        else 
            % ero in banda e rimando in banda
            s = s+1;
        end
        
    else  % scambio
        if old == 0
            % ero nel vuoto e passo a banda
            s = 1;
            v = v+1;
        else
            % ero in bada e passo a vuoto 
            s_v(v) = s;
        end
    end
    
    old = vect(j);

end

n = size(s_v,2);
