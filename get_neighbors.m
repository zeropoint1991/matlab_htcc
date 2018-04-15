function B = get_neighbors(n,lamda,t)
%邻居是包括自己的，n是向量个数，t是邻居个数
distance_matrix = zeros(n,n); %放距离

B = zeros(n,t);               %放邻居标识

if n <= t
    error('N is smaller than T');
end

for i = 1:n
    for j = (i+1):n
        distance_matrix(i,j) = norm(lamda(i,:)-lamda(j,:));
        distance_matrix(j,i) = distance_matrix(i,j);
    end
    [~,near_index] = sort(distance_matrix(i,:));
    B(i,:) = near_index(1:t);
end
end
    

