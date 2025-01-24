function SOFIimag = getSOFI(img, cumulant,batchSize)


[sof, ~, ~] = SOFIAnalysis(img,cumulant,batchSize, 1);

SOFIimag = abs(sof);
% SOFIimag = sof-min(sof,[],"all");
end