function result = get_lambda(m,n)
    file_name = strcat('W',num2str(m),'D_',num2str(n),'.dat');
    result = importdata(file_name);
end