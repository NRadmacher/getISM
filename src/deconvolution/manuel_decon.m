function manuel_decon(img, psf, options)
%MAUEL_DECON Decolvolves image img with Lucy-Richardson. Iterations are
%maulaly added via text filed  c_map, name, title_n, sb_lenght, IM_R, reso, save

arguments
    img (:,:) double = checkerboard(256);
    psf (:,:) double = ones(11,11)
    %Imgae Colormap
    options.c_map (256,3) = hot
    %name of saved image
    options.s_name = 'decon'
    % Image titel name
    options.t_name = 'ISM + deconvolution'
    %sclae bar lenght
    options.sb_lenght = 2
    %Pixel size in µm
    options.IM_R = 0.05
    %pixel binning
    options.pix_bin = 1;
    %do resolution analysis
    options.reso = 0
    %save image to 600 dpi png
    options.save = 1
end

%Deconvolution Parameters
params.dampener = 0;
params.iter = 1;

%Store starting image
data.Image       = img;
%Current deconvolved image
data.decon_img   = img;
%PSF
data.psf = psf;

% Opten UiFigure to show deconvolution
handles.fig = uifigure('HandleVisibility', 'on');
% Grid to place UI elemnts
handles.gl = uigridlayout(handles.fig,[3 3]);
handles.gl.RowHeight = {22, 22, 22, '1x'};
handles.gl.ColumnWidth = {'1x', 'fit', 50};

%Axes for plotting image
handles.ax = uiaxes(handles.gl);
% Positon of image on the grid
handles.ax.Layout.Row =  [1 4];
handles.ax.Layout.Column = 1;

%Text field to enter numner of LR iterations
handles.iter  = uieditfield(handles.gl,'numeric','ValueChangedFcn',@(txt,event)set_iter(options, event));
handles.iter.Layout.Row = 1;
handles.iter.Layout.Column = 3;
handles.iter_label = uilabel(handles.gl, 'Text','Iterations');
handles.iter_label.Layout.Row = 1;
handles.iter_label.Layout.Column = 2;

%Text field to enter dampener value of LR
handles.dampen  = uieditfield(handles.gl,'numeric','ValueChangedFcn',@(txt,event)set_dampener(options, event));
handles.dampen.Layout.Row = 2;
handles.dampen.Layout.Column = 3;
handles.dampen_label = uilabel(handles.gl, 'Text','Dampener');
handles.dampen_label.Layout.Row = 2;
handles.dampen_label.Layout.Column = 2;

%Button to save and plot image
handles.btn = uibutton(handles.gl,'push', 'ButtonPushedFcn', @(btn,event) print_image(options));
handles.btn.Layout.Row = 3;
handles.btn.Layout.Column = 3;

%At the start of the funtion show ISM image
colormap(handles.fig, options.c_map)
imagesc(handles.ax, data.Image);
set(handles.ax,'DataAspectRatio', [1,1,1], ...
'PlotBoxAspectRatio',[1 1 1], ...
'XDir','normal', ...
'YDir','reverse');
set(findall(handles.fig,'-property','FontSize'),'FontSize',13)
c = colorbar(handles.ax);
c.Label.String = 'Count';
handles.ax.Visible = 'off';


function deconvolve(options)
    % Deconvolution parameters
    iter = params.iter;
    dampar = params.dampener;
    
    % Actual devonvolution with PSFS
    new_img             = deconvlucy(data.Image, data.psf, iter, dampar);
%     new_img             = deconvwnr(handles.Image, psf, iter);
%     [new_img,~]         = deconvblind(handles.Image, psf);
%     new_img             = deconvreg(handles.Image, psf);
    handles.decon_img   = new_img;
    
    %Plot deconvolved image in figure
    colormap(handles.fig, options.c_map)
    imagesc(handles.ax, handles.decon_img);
    set(handles.ax,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1], ...
    'XDir','normal', ...
    'YDir','reverse');
    c_bar = colorbar(handles.ax);
    c_bar.Label.String = 'Count';
    set(findall(handles.fig,'-property','FontSize'),'FontSize',13)
    handles.ax.Visible = 'off';
end

function print_image(options)
    reso_line = [[42 92]; [98 98]];
    % Plot and save using img_plot
    img_plot(handles.decon_img, options.c_map,...
            options.s_name, options.t_name, ...
            options.sb_lenght, options.IM_R, ...
            options.pix_bin, reso_line, options.reso, options.save)
end

function set_dampener(options, event)
   params.dampener = event.Value;
   deconvolve(options);
end

function set_iter(options, event)
   params.iter = event.Value;
   deconvolve(options);
end

end

