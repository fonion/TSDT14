function [ACF_matrix] = ACF_bart(y,k)
    ACF_matrix = zeros(1,2*k);
    counter=1;
    for i = -k:k
        ACF_matrix(counter) = ACF_help(y,i);
        counter = counter + 1;
    end
end