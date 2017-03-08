function do_plot(phi,figID,varargin)

global npar dat snq

figure(1+figID)
plot(npar.xf,reshape(phi,npar.ndofs,1),'-','LineWidth',1);
title(sprintf('scalar flux, problem %d',dat.pb_ID));
xlabel('x'); ylabel('scalar flux');
% filename=sprintf('scal_%s_sn%i.png',tit,sn);
% print('-dpng',filename);

if length(varargin)>=1
    
    figure(2+figID)
    E = varargin{1};
    plot(npar.xf,reshape(E,npar.ndofs,1),'-','LineWidth',1);
    title(sprintf('Eddington tensor, problem %d',dat.pb_ID));
    xlabel('x'); ylabel('E');
    
    if length(varargin)==2
        psi = varargin{2};
        for idir=1:snq.n_dir
            figure(figID+2+idir)
            plot(npar.xf,reshape(psi(:,:,idir),npar.ndofs,1),'-','LineWidth',1);
            title(sprintf('angular flux %d, problem %d',idir,dat.pb_ID));
            xlabel('x'); ylabel('angular flux');
        end
    end
end

return
end
