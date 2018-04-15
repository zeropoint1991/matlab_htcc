

function chromosome = update_neighbors_NPOP(chromosome,offspring,p,nr,lambda,z,v,m)

%p is the neighborhood pop or the whole pop
%this function is to use offspring to update the chrommosome with nr times

if (offspring(end)==0) %if the offspring is feasible, then replace the feasible population
    c = 0 ;
    while (c < nr)&&~isempty(p) 
        p_len = length(p);
        rand_index = randi(p_len,1,1);  %randomly select a subproblem to replace
        offspring_fit=teFit(offspring(v+1:v+m),p(rand_index)); % for offspring
        if(~isempty(chromosome{p(rand_index)}.feasiblepop))
            parent_fit=teFit(chromosome{p(rand_index)}.feasiblepop(v+1:v+m),p(rand_index));
            if(offspring_fit<parent_fit)
                c=c+1;
                chromosome{p(rand_index)}.feasiblepop=offspring;
            end
        else
            c=c+1;
            chromosome{p(rand_index)}.feasiblepop=offspring;
        end
        p(rand_index) = [];
    end
else %if the offspring is infeasible, then replace the infeasible population
    c = 0 ;
    tempp=p;
    while (c < nr)&&~isempty(p) 
        p_len = length(p);
        rand_index = randi(p_len,1,1);  %randomly select a subproblem to replace
        offspring_fit=teFit(offspring(v+1:v+m),p(rand_index)); % for offspring
        if(~isempty(chromosome{p(rand_index)}.infeasiblepop1))
            parent_fit=teFit(chromosome{p(rand_index)}.infeasiblepop1(v+1:v+m),p(rand_index));
            if(offspring_fit<parent_fit)
                c=c+1;%0.001*offspring(end)+0.001*chromosome{p(rand_index)}.infeasiblepop1(end)
                chromosome{p(rand_index)}.infeasiblepop1=offspring;
            end
        else
            c=c+1;
            chromosome{p(rand_index)}.infeasiblepop1=offspring;
        end
        p(rand_index) = [];
    end
    
    c=0;
    p=tempp;
    while (c < nr)&&~isempty(p) 
        p_len = length(p);
        rand_index = randi(p_len,1,1);  %randomly select a subproblem to replace
        offspring_fit=teFit(offspring(v+1:v+m),p(rand_index)); % for offspring
        if(~isempty(chromosome{p(rand_index)}.infeasiblepop2))
            parent_fit=teFit(chromosome{p(rand_index)}.infeasiblepop2(v+1:v+m),p(rand_index));
            if(offspring(end)+0.01*offspring_fit<chromosome{p(rand_index)}.infeasiblepop2(end)+0.01*parent_fit)
                c=c+1;
                chromosome{p(rand_index)}.infeasiblepop2=offspring;
            end
        else
            c=c+1;
            chromosome{p(rand_index)}.infeasiblepop2=offspring;
        end
        p(rand_index) = [];
    end
%     while (c < nr)&&~isempty(p) 
%         p_len = length(p);
%         rand_index = randi(p_len,1,1);  %randomly select a subproblem to replace
% 
%         if(isempty(chromosome{p(rand_index)}.infeasiblepop))
%             chromosome{p(rand_index)}.infeasiblepop=offspring(:,1:v+m+1);
%             chromosome{p(rand_index)}.size=1;
%             c=c+1;p(rand_index) = [];
%             continue;
%         end
%         subpop=chromosome{p(rand_index)}.infeasiblepop;
%         subpopsize=chromosome{p(rand_index)}.subpopsize;
%         %calculate two objectives for 1. the violence: offspring(end)
%         %2. the tchebycheff function: teFit(offspring(v+1:v+m))
%         offspring(1,v+m+2)=teFit(offspring(v+1:v+m),p(rand_index)); % for offspring
%         for i=1:size(subpop,1)
%             subpop(i,v+m+2)=teFit(subpop(v+1:v+m),p(rand_index));%for chromosome{rand_index}.pop
%         end
%         % if the offspring is dominated by the original solution, than continue
%         % if the offspring is dominate the original solution, than replace it
%         % and add nr by 1
%         % if the offspring is non-dominated with the original solutions, than check crowding distance 
%         % parato dominance comparison
%         idx=1;mark=[];
%         for i=1:size(subpop,1)
%             flagd=dominated_relationship(offspring,subpop,2,v+2);
%             if flagd==1 % offspring dominates subpop
%                 mark(idx)=i;
%                 idx=idx+1;%mark subpop(i) as a dominated solution
%             elseif (flagd==2 || flagd==3) % subpop dominates (or equal to) offspring
%                 break;% do not replace in this subproblem
%             elseif flagd==4 % subpop non-dominated with offspring
%                 continue;
%             end
%         end
%         if (flagd==2 || flagd==3 )
%             p(rand_index) = [];continue; % do not replace this subproblem
%         end
%         %delete the marked dominsted solution in subproblem
%         for i=1:size(mark)
%             subpop(mark(i),:)=[];
%         end
%         %add the offspring to subproblem
%         subpop(size(subpop,1)+1,:)=offspring;
%         c=c+1;
%         %if archive size larger than population size, than use crowding distance to
%         %compare
%         if size(subpop,1)>subpopsize
%             %%
%             M=2;
%             [~,rank]=sort(subpop(:,end-M)); %sort by violence
%             subpop=subpop(rank,:);
%             Fit=subpop(:,end-M+1:end);
%             archsize=size(subpop,1);
%             distance=zeros(archsize,1);
%             distance(1)=inf;
%             distance(end)=inf;
%             %normalize the objective function
%             for c=1:M
%                 fix_max=max(Fit(:,c));
%                 fix_min=min(Fit(:,c));
%                 Fit(:,c)=(Fit(:,c)-fix_min)/(fix_max-fix_min);
%             end
%             %calculate distance
%             for i=2:archsize-1
%                 for c=1:M
%                     distance(i)=distance(i)+abs(Fit(i-1,c)-Fit(i+1,c));
%                 end
%             end
%             %delete the small crowding distance
%             [~,d_rank]=sort(distance);
%             subpop(d_rank(1),:)=[];
%             
%         end
%         chromosome{p(rand_index)}.infeasiblepop=subpop(:,1:end-1);
%         chromosome{p(rand_index)}.size=size(subpop,1);
%         p(rand_index) = [];
%     end
end
% tch method
function y = teFit(fit,i)
    lambda(i,:) = max(lambda(i,:),repmat(1e-6,size(lambda(i,:))));
    y =max(1./lambda(i,:).*abs(fit - z')); 
end
% function y = teFit_diversity(fit,i)
%     lambda(i,:) = max(lambda(i,:),repmat(1e-6,size(lambda(i,:))));
%     y =max(1./lambda(i,:).*abs(fit - z')); 
% end
function y = teFit_diversity(fit,i)
    d1=sum(1./lambda(i,:).*abs(fit - z'));
    y = abs(d1) / sqrt(sum(lambda(i,:).^2));
end
end