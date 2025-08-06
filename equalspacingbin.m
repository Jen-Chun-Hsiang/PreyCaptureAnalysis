function boundaries = equalspacingbin(data, num_bin)
    data = data(:);
    bin_size = floor(length(data)/num_bin);
    sdata = sort(data, 'ascend');
    boundaries = nan(num_bin, 3);
    for i = 1:num_bin
        if i == 1
            boundaries(i, :) = [sdata(1), sdata(bin_size*i), bin_size];
        elseif i == num_bin
            boundaries(i, :) = [sdata(bin_size*(i-1)+1), sdata(end), length(data)-bin_size*(i-1)];
        else
            boundaries(i, :) = [sdata(bin_size*(i-1)+1), sdata(bin_size*i), bin_size];
        end
    end
end