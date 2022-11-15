function im = LT_image(taus,tags, spectrum, threshold, gamma, lims)
% Scale lifetime image by intensity
%taus                   Lifetime
%tags                   Intensity
%spectrum               colormap for lifetimes
%threshold              intensity threshold lower bound
%gamma                  gamma for intensity scaling
%lims                   limits for lifelime colormap

ind = (tags < threshold);    % tags: pixelwise intensity 
taus(ind) = 0;                % set lt to where intensity is too low

ind = (tags >= threshold);
intens = squeeze((ind.*tags).^gamma); 
intens = intens./max(intens(:));

%intensity in every rgb cannel
tmp = repmat((intens), [1 1 3]);

%set nans to zero and applay lifetime limits
val  = taus(:);
val(isnan(val))  = 0;
val(val<lims(1)) = 0;
val(val>lims(2)) = lims(2);

% assign colormap value according to lifetime
k      = 1 + round((val-lims(1))./(lims(2)-lims(1)).*(size(spectrum,1)-1));
% set values resultion vom lifetimes smalder than lim to 1
k(k<1) = 1;

%assigne rgb value to scaled value k;
im = spectrum(k,:);                 

im = reshape(im,size(tmp,1),size(tmp,2),3); % reshape color list back into 2D image
im = im.*tmp;                               % scale colors with intensity

end



