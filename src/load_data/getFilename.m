function name = getFilename(fname)

%Image title and name from file name
tmp_name        = strsplit(fname, '\');

%find date in filepath
i = 1;
date = [];
datePat = digitsPattern(6);
while isempty(date)
    if i == size(tmp_name,2)
        ME = MException('Filemane:nodate','file %s \ncontains no date', fname);
        date = {'000000'};
        fprintf('Warning no date found in file name!\n')
        break;
    end
    date = extract(tmp_name{end-i},datePat);
    i = i+1;
end
date = date{1};

%check if file was recorded by Symphotime and add ws name to image name
if contains(tmp_name{end-i+2},'.')
    wsName  = strsplit(tmp_name{end-i+2},'.');
    date    = append(date,' ',wsName{1});
    %check if file is part of Group Measurment and add name of GM to image
    %name
    if (i-3) > 0
        gmName  = tmp_name{end-i+3};
        date    = append(date,' ',gmName);
    end
end

name        = tmp_name{end};
name        = strsplit(name, '.');
name        = name{end-1};
name        = append(date,' ',name);

end