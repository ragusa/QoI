function do_plot(phi,varargin)

global npar dat

figure(1)
plot(npar.xf,reshape(phi,npar.ndofs,1),'-','LineWidth',1);
title(sprintf('scalar flux, problem %d',dat.pb_ID));
xlabel('x'); ylabel('scalar flux');
% filename=sprintf('scal_%s_sn%i.png',tit,sn);
% print('-dpng',filename);

if length(varargin)==1
    for idir=1:snq.n_dir
        figure(10+idir)
        plot(npar.xf,reshape(psi(:,:,idir),npar.ndofs,1),'-','LineWidth',1);
        title(sprintf('angular flux %d, problem %d',idir,dat.pb_ID));
        xlabel('x'); ylabel('angular flux');
    end
end

return
end
