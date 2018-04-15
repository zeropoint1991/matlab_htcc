function result = overall_cv(cv)
cv(cv > 0) = 0;cv = abs(cv);
result = sum(cv,1); %ÐÐÏà¼Ó
end


