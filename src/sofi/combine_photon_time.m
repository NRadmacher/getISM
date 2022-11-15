function y = combine_photon_time(x, tol)
%COMBINE_PHOTON_TIME removes scan time from photen stream of ISM images
%with multiple scans via recursion
%   find first index where differenc is bigger than tol

time_diff = diff(x);
ind = find(time_diff > tol);

if isempty(ind)
    y = x;
    return
else
    x(ind(1)+1:end) =  x(ind(1)+1:end) - time_diff(ind(1));
    y = combine_photon_time(x, tol);
    return
end

end