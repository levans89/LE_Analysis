function z = getRows(unique_1, unique_2, unique_3, unique_4)

    switch nargin
        case 1
            z = unique_1;
        case 2
            z = ismember(unique_1, unique_2);
        case 3
            a = ismember(unique_1, unique_2);
            b = ismember(unique_1, unique_3);
            z = a & b;
        case 4
            a = ismember(unique_1, unique_2);
            b = ismember(unique_1, unique_3);
            c = ismember(unique_1, unique_4);
            z = a & b & c;
    end
    
end
