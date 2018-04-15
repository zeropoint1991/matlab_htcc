function [y,P] = select_and_reproduct(pop,B,i,delta,n_dim,lu,pcrossover,de_factor,miu,pmutation)
%% with the probability of delta to select the whole population as the neighbor of individual i
if rand < delta
    P = B(i,:);
else
    P = 1:size(pop,1);
end

p_len = length(P);

%% select two individuals based on tournament selection
%  select two individuals randomly
parent = zeros(2,pop{1}.dim);
rand_index = randi(p_len,1,2);
while rand_index(1) == rand_index(2)
    rand_index = randi(p_len,1,2);
end

if((rand<0.5 || isempty([pop{P(rand_index(1))}.infeasiblepop1;pop{P(rand_index(1))}.infeasiblepop2])) && ~isempty(pop{P(rand_index(1))}.feasiblepop))
    parent(1,:) = pop{P(rand_index(1))}.feasiblepop;
else
    if((rand<0.5 || isempty(pop{P(rand_index(1))}.infeasiblepop2)) && ~isempty(pop{P(rand_index(1))}.infeasiblepop1))
        parent(1,:) = pop{P(rand_index(1))}.infeasiblepop1;
    else
        parent(1,:) = pop{P(rand_index(1))}.infeasiblepop2;
    end
end

if((rand<0.5 || isempty([pop{P(rand_index(2))}.infeasiblepop1;pop{P(rand_index(2))}.infeasiblepop2])) && ~isempty(pop{P(rand_index(2))}.feasiblepop))
    parent(2,:) = pop{P(rand_index(2))}.feasiblepop;
else
    if((rand<0.5 || isempty(pop{P(rand_index(2))}.infeasiblepop2)) && ~isempty(pop{P(rand_index(2))}.infeasiblepop1))
        parent(2,:) = pop{P(rand_index(2))}.infeasiblepop1;
    else
        parent(2,:) = pop{P(rand_index(2))}.infeasiblepop2;
    end
end
%% differential evolution
if ((rand < 0.6 || isempty([pop{i}.infeasiblepop1;pop{i}.infeasiblepop2])) && ~isempty(pop{i}.feasiblepop)) % generate child solution by feasiblepop
%     if(~isempty(pop{P(rand_index(1))}.feasiblepop))
%         parent(1,:) = pop{P(rand_index(1))}.feasiblepop;
%     end
% 
%     if(~isempty(pop{P(rand_index(2))}.feasiblepop))
%         parent(2,:) = pop{P(rand_index(2))}.feasiblepop;
%     end
    y = de_crossover(pop{i}.feasiblepop,parent(1,:), parent(2,:), n_dim, lu, pcrossover, de_factor);%pop{1}.flag=0;
else
    if((rand<0.5 || isempty(pop{i}.infeasiblepop2)) && ~isempty(pop{i}.infeasiblepop1))
%         if(~isempty(pop{P(rand_index(1))}.infeasiblepop1))
%             parent(1,:) = pop{P(rand_index(1))}.infeasiblepop1;
%         end
% 
%         if(~isempty(pop{P(rand_index(2))}.infeasiblepop1))
%             parent(2,:) = pop{P(rand_index(2))}.infeasiblepop1;
%         end
        y = de_crossover(pop{i}.infeasiblepop1,parent(1,:), parent(2,:), n_dim, lu, pcrossover, de_factor);%pop{1}.flag=1;
    else
%         if(~isempty(pop{P(rand_index(1))}.infeasiblepop2))
%             parent(1,:) = pop{P(rand_index(1))}.infeasiblepop2;
%         end
% 
%         if(~isempty(pop{P(rand_index(2))}.infeasiblepop2))
%             parent(2,:) = pop{P(rand_index(2))}.infeasiblepop2;
%         end
        y = de_crossover(pop{i}.infeasiblepop2,parent(1,:), parent(2,:), n_dim, lu, pcrossover, de_factor);%pop{1}.flag=2;
    end
end

%% mutation operators
y1 = y;
rand_value = rand;
if rand_value < 0.5
    theta = (2*rand_value)^(1/(1+miu)) - 1;
else
    theta = 1 - (2 - 2*rand_value)^(1/(1+miu));
end
rand_mut = rand(n_dim,1);
muted = rand_mut < pmutation;
y = y1 + theta*(lu(2,:) - lu(1,:));
y(~muted) = y1(~muted);


%% repair
%% make sure that generated offspring in their ranges.
y = min(max(y,lu(1,:)),lu(2,:));


end

