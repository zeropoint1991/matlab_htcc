function B = get_neighbors(n,lamda,t)
%�ھ��ǰ����Լ��ģ�n������������t���ھӸ���
distance_matrix = zeros(n,n); %�ž���

B = zeros(n,t);               %���ھӱ�ʶ

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
    

