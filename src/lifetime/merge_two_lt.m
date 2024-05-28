function merge_two_lt(lt_amp_img, names, LT_name, pattern_tau, IM_R)
    
    lt_rgb      = zeros(size(lt_amp_img));
    % adjust contrast
    lt_rgb(:,:,1) = adapthisteq(lt_amp_img(:,:,1)./max(lt_amp_img(:,:,1),[],'all'),'ClipLimit',0.005);%red
    lt_rgb(:,:,2) = adapthisteq(lt_amp_img(:,:,2)./max(lt_amp_img(:,:,2),[],'all'),'ClipLimit',0.005);%green

    save_name = compose('%s multicolor %0.1f %0.1f', LT_name, pattern_tau(1), pattern_tau(2));
    title = append(names{1,1}, ': ', names{2,1}, ', ', names{1,2}, ': ', names{2,2});

    rgb_img_plot(lt_rgb, save_name{1}, title, 4, IM_R, 1)
end