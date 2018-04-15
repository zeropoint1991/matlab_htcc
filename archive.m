function arch = archive(arch,population,popsize,M)

temp_fea=population(:,end)==0;
temp_population=population(temp_fea,:);
if isempty(arch)
    arch=temp_population;
elseif ~isempty(arch)&&~isempty(temp_population)
    arch=union(arch,temp_population,'rows');
end
% parato dominance comparison
if size(arch,1)>1
    arch=pare(arch,M);
end
%if archive size larger than population size, than use crowding distance to
%compare
if size(arch,1)>popsize
    arch=cd_count(arch,popsize,M);
end
end


%parato dominance
function arch=pare(arch,M)
popsize=size(arch,1);
Fit=arch(:,end-M:end-1);
nomdom_label=true(popsize,1);
for i = 1:popsize
    for j = 1:popsize
        
        % If both of the individual is feasible    CDP规则判断
        less = 0;
        equal = 0;
        more = 0;                         %这三个记录关系
        
        for k = 1:M
            if Fit(i,k) < Fit(j,k)
                less = less+1;
            elseif Fit(i,k) == Fit(j,k)
                equal = equal+1;
            else
                more = more+1;
            end
        end
        
        % Individual i is dominated
        if less == 0 && equal ~= M       %个体i被j支配的情况
            nomdom_label(i)=false;
            % Individual j is dominated by i
        elseif more == 0 && equal ~= M  %个体i支配j支配的情况
            nomdom_label(j)=false;
        end
    end
end
arch=arch(nomdom_label,:);


end

function arch=cd_count(arch,popsize,M)
if M==2
    [~,rank]=sort(arch(:,end-M));
    arch=arch(rank,:);
    Fit=arch(:,end-M:end-1);
    archsize=size(arch,1);
    distance=zeros(archsize,1);
    distance(1)=inf;
    distance(end)=inf;
    
    for i=2:archsize-1
        for c=1:M
            distance(i)=distance(i)+abs(Fit(i-1,c)-Fit(i+1,c));
        end
    end
    %%%%%%%%%%%%%%%%%% 归一化 %%%%%%%%%%%%%%%%%%%%%%%%%%
    for c=1:M
        fix_max=max(Fit(:,c));
        fix_min=min(Fit(:,c));
        Fit(:,c)=(Fit(:,c)-fix_min)/(fix_max-fix_min);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while archsize>popsize
        [~,d_rank]=sort(distance);
        label=d_rank(1);
        for c=1:M
            distance(label-1)=distance(label-1)+abs(Fit(label-1,c)-Fit(label,c));
            distance(label+1)=distance(label+1)+abs(Fit(label,c)-Fit(label+1,c));
        end
        arch(d_rank(1),:)=[];
        Fit(d_rank(1),:)=[];
        distance(label)=[];
        archsize=archsize-1;
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%% M>=3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    Fit=arch(:,end-M:end-1);
    archsize=size(arch,1);
    distance=zeros(archsize,1);
    new_arch=[arch,distance];
    for c=1:M
        [~,rank]=sort(Fit(:,c));
        fix_max=Fit(rank(end),c);
        fix_min=Fit(rank(1),c);
        Fit(:,c)=(Fit(:,c)-fix_min)/(fix_max-fix_min);
        new_arch(rank(1),end)=inf;
        new_arch(rank(end),end)=inf;
        Fit=Fit(rank,:);
        new_arch=new_arch(rank,:);
        temp_fit=Fit(:,c);
        for t=2:archsize-1
        new_arch(t,end)=new_arch(t,end)+abs(temp_fit(t-1)-temp_fit(t+1));
        end
    end
    [~,d_rank]=sort(new_arch(:,end),'descend');
    arch=new_arch(d_rank(1:popsize),1:end-1);
        
    
end


end