function [dummy] = process_bar(i,iend,info)

if nargin == 2
    info = '';
end

if i == 1
    clc
    fprintf('>                    |  %s \n',info);
elseif i == round(iend/20)
    clc
    fprintf('->                   |  %s \n',info);
elseif i == round(2*(iend/20))
    clc
    fprintf('-->                  |  %s \n',info);
elseif i == round(3*(iend/20))
    clc
    fprintf('--->                 |  %s \n',info);
elseif i == round(4*(iend/20))
    clc
    fprintf('---->                |  %s \n',info);
elseif i == round(5*(iend/20))
    clc
    fprintf('----->               |  %s \n',info);
elseif i == round(6*(iend/20))
    clc
    fprintf('------>              |  %s \n',info);
elseif i == round(7*(iend/20))
    clc
    fprintf('------->             |  %s \n',info);
elseif i == round(8*(iend/20))
    clc
    fprintf('-------->            |  %s \n',info);
elseif i == round(9*(iend/20))
    clc
    fprintf('--------->           |  %s \n',info);
elseif i == round(10*(iend/20))
    clc
    fprintf('---------->          |  %s \n',info);
elseif i == round(11*(iend/20))
    clc
    fprintf('----------->         |  %s \n',info);
elseif i == round(12*(iend/20))
    clc
    fprintf('------------>        |  %s \n',info);
elseif i == round(13*(iend/20))
    clc
    fprintf('------------->       |  %s \n',info);
elseif i == round(14*(iend/20))
    clc
    fprintf('-------------->      |  %s \n',info);
elseif i == round(15*(iend/20))
    clc
    fprintf('--------------->     |  %s \n',info);
elseif i == round(16*(iend/20))
    clc
    fprintf('---------------->    |  %s \n',info);
elseif i == round(17*(iend/20))
    clc
    fprintf('----------------->   |  %s \n',info);
elseif i == round(18*(iend/20))
    clc
    fprintf('------------------>  |  %s \n',info);
elseif i == round(19*(iend/20))
    clc
    fprintf('-------------------> |  %s \n',info);
elseif i == round(20*(iend/20))
    clc
    fprintf('-------------------->|  %s \n',info);
    
end

dummy = 1;

end

