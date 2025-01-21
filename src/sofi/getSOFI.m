function SOFIimag = getSOFI(img, cumulant,batchSize)


[sof, im0, soffull] = SOFIAnalysis(img,cumulant,batchSize, 1);

sof = abs(sof);
SOFIimag = sof-min(sof,[],"all");
end