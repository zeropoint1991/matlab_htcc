function child = de_crossover(parent1, parent2,parent3, n, lu, pcrossover, de_factor)

parent_1 = parent1(1 : n);
parent_2 = parent2(1 : n);
parent_3 = parent3(1 : n);

% pcrossover=0.9 is the probability of crossover
cross_idx1 = rand(1,n) < pcrossover;
cross_idx1(randi(n)) = 1; % at least one decision component to crossover

child = parent_1 + de_factor * (parent_3 - parent_2);
child(~cross_idx1) = parent_1(~cross_idx1);
child =  max(min(lu(2,:),child),lu(1,:));
end
