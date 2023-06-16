function [eid, psf] = get_psf(rhofield, pix)
%% PSF calculation
NA = 1.49;
lamex = 0.64; % excitation wavelength in mum
lamem = 0.70; % emission wavelength in mum
fd = 1.8e3; % objective focal length in mum
tubelens = 180e3; % tube lens focal length in mum
mag = 2*tubelens/fd; % magnification in pinhole plane
av = 1000/2; % pinhole radius in mum
zpin = 0; % pinhole position raltive to focal plane
% pix = 0.05; % step size of scanning in mum
% rhofield = [0 0.8]; % radial range in mum to be considered
zfield = [-1 1]; % axial range in mum to be considered
n0 = 1.5; % ref index of immersion medium = oil/glass
n = n0;
n1 = n0;
d0 = [];
d = 1;
d1 = [];
over =  inf; % diffraction-limited focusing
focpos = 0; % focusing on surface

exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos);
exc1 = RotateEMField(exc, pi/2);
exc.fxc = exc.fxc + 1i*exc1.fxc;
exc.fxs = exc.fxs + 1i*exc1.fxs;
exc.fyc = exc.fyc + 1i*exc1.fyc;
exc.fys = exc.fys + 1i*exc1.fys;
exc.fzc = exc.fzc + 1i*exc1.fzc;
exc.fzs = exc.fzs + 1i*exc1.fzs;

mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, [], 0);

%% cross section of excitation intensity distribution in focal plane
[feld, phi, rr, ~] = FocusImage3D(exc.rho,exc.z,cat(4,cat(3,exc.fxc,exc.fxs),cat(3,exc.fyc,exc.fys),cat(3,exc.fzc,exc.fzs)));
[~,ind] = min(abs(exc.z(1,:)));

%% pixelated cross section of excitation intensity distribution in focal plane
xmax = floor(max(rr(:))/sqrt(2)/pix);
[xx,yy]= meshgrid(pix*(-xmax:xmax),pix*(-xmax:xmax));
psi = angle(xx+1i*yy)+pi;
rad = sqrt(xx.^2+yy.^2);
eid = interp2(squeeze(phi(:,ind,:)),squeeze(rr(:,ind,:)),...
    squeeze(sum(abs(feld(:,ind,:,:)).^2,4)),psi,rad,'cubic');

%% cross section of PSF in focal plane
[feld, phi, rr, ~] = FocusImage3D(exc.rho,exc.z,mdf.volx+mdf.voly);
[~,ind] = min(abs(exc.z(1,:)));

%% pixelated image of PSF cross section in focal plane
xmax = floor(max(rr(:))/sqrt(2)/pix);
[xx,yy]= meshgrid(pix*(-xmax:xmax),pix*(-xmax:xmax));
psi = angle(xx+1i*yy)+pi;
rad = sqrt(xx.^2+yy.^2);
psf = interp2(squeeze(phi(:,ind,:)),squeeze(rr(:,ind,:)),...
    squeeze(sum(abs(feld(:,ind,:,:)).^2,4)),psi,rad,'cubic');
end
