
function range = decision_range(name,dim)

range = ones(dim,2);      %行代表变量数，列代表上下限

switch name
    case {'UF1','UF2','UF5','UF6','UF7','CF2'}
        range(1,1)      =  0;
        range(2:dim,1)  = -1;
    case {'UF3','LIRCMOP1','LIRCMOP2','LIRCMOP3','LIRCMOP4','LIRCMOP5','LIRCMOP6','LIRCMOP7','LIRCMOP8','LIRCMOP9','LIRCMOP10','LIRCMOP11','LIRCMOP12','LIRCMOP13','LIRCMOP14'}
        range(:,1)      =  0;
    case {'UF4','CF3','CF4','CF5','CF6','CF7'}
        range(1,1)      =  0;
        range(2:dim,1)  = -2;
        range(2:dim,2)  =  2;
    case {'UF8','UF9','UF10','CF9','CF10'}
        range(1:2,1)    =  0;
        range(3:dim,1)  = -2;
        range(3:dim,2)  =  2;
    case 'CF1'
        range(:,1)      =  0;
    case {'CF8'}
        range(1:2,1)    =  0;
        range(3:dim,1)  = -4;
        range(3:dim,2)  =  4;
    case{'CTP1','CTP2','CTP3','CTP4','CTP5','CTP6','CTP7','CTP8'}
        range(1,1)      =  0;
        range(2:dim,1)  = -5;
        range(2:dim,2)  = 5;
    case 'BNH'
        range([1,2],1) = 0;
        range(1,2) = 5;
        range(2,2) = 3;
    case 'TNK'
        range(:,1)      = 0;
        range(:,2)      = pi;
    case 'SRN'
        range(:,1)      = -20;
        range(:,2)      = 20;
    case 'OSY'
        range([1,2,4,6],1) = 0;
        range([1,2,6],2)   = 10;
        range([3,5],2)     = 5;
        range(4,2)         = 6;
    case 'CONSTR'
        range(1,1) = 0.1;
        range(2,1) = 0;
        range(2,2) = 5;
        
end
end