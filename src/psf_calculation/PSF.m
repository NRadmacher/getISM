function int = PSF(NA,fd,lamex,pix,over)

    rhofield = [0 2];
    zfield = 0;
    %refractiv index
    n0 = 1.52;
    n = n0;
    n1 = n0;
    d0 = [];
    d = max(zfield);
    d1 = [];
    focpos = 0;
    atf = [];
    resolution = [];
    ring = [];
    maxm = [];

    exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, ...
        over, focpos, atf, resolution, ring, maxm);
    exc1 = RotateEMField(exc, pi/2);
    exc.fxc = exc.fxc + 1i*exc1.fxc;
    exc.fxs = exc.fxs + 1i*exc1.fxs;
    exc.fyc = exc.fyc + 1i*exc1.fyc;
    exc.fys = exc.fys + 1i*exc1.fys;
    exc.fzc = exc.fzc + 1i*exc1.fzc;
    exc.fzs = exc.fzs + 1i*exc1.fzs;
    
    nn = floor(max([0 5])/pix/sqrt(2));
    [xx, yy] = meshgrid(pix*(-nn:nn),pix*(-nn:nn));
    pp = atan2(yy,xx);
    rr = sqrt(xx.^2+yy.^2);
    int = interp1(exc.rho,abs(exc.fxc(:,1,1)).^2+abs(exc.fyc(:,1,1)).^2+abs(exc.fzc(:,1,1)).^2,rr,'pchip',0);
    for j=1:size(exc.fxs,3)
        fx = interp1(exc.rho,exc.fxc(:,1,j+1),rr,'pchip',0).*cos(j*pp) + interp1(exc.rho,exc.fxs(:,1,j),rr,'pchip',0).*sin(j*pp);
        fy = interp1(exc.rho,exc.fyc(:,1,j+1),rr,'pchip',0).*cos(j*pp) + interp1(exc.rho,exc.fys(:,1,j),rr,'pchip',0).*sin(j*pp);
        fz = interp1(exc.rho,exc.fzc(:,1,j+1),rr,'pchip',0).*cos(j*pp) + interp1(exc.rho,exc.fzs(:,1,j),rr,'pchip',0).*sin(j*pp);
        int = int + abs(fx).^2 + abs(fy).^2 + abs(fz).^2;
    end
end