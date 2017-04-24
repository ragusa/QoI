function do_plot(phi,my_legend,figID,forward,varargin)

global npar dat snq

if forward
    forward_txt='Forward';
else
    forward_txt='Adjoint';
end
% dealing with legend
init=false;
if ~isfield(dat,'leg')
    dat.leg{(figID+1)}=char(my_legend);
    init=true;
else
    try
        a = length(dat.leg{(figID+1)});
    catch ME
        switch ME.identifier
            case 'MATLAB:badsubscript'
                dat.leg{(figID+1)}=char(my_legend);
                init=true;
            otherwise
                ME
                ME.identifier
                error('ssss')
        end
    end
end
if ~init
    dat.leg{(figID+1)}=char(dat.leg{(figID+1)},my_legend);
end
% markers
my_markers = ['+', 'o', '*', 's', 'd', 'v', '.', '^', '<', '>', 'p', 'h'];
[n1,n2]=size(dat.leg{(figID+1)});
sss=sprintf('%s-',my_markers(n1));
figure(1+figID); hold all;
plot(npar.xf,reshape(phi,npar.ndofs,1),sss,'LineWidth',1);
title(sprintf('%s scalar flux, problem %d',forward_txt,dat.pb_ID));
xlabel('x'); ylabel('scalar flux');
% filename=sprintf('scal_%s_sn%i.png',tit,sn);
% print('-dpng',filename);
legend(dat.leg{figID+1},'Location','Best');
hold off;

if length(varargin)>=1
    
    figure(2+figID)
    E = varargin{1};
    plot(npar.xf,reshape(E,npar.ndofs,1),'+-','LineWidth',1);
    title(sprintf('%s Eddington tensor, problem %d',forward_txt,dat.pb_ID));
    xlabel('x'); ylabel('E');
    
    if length(varargin)==2
        psi = varargin{2};
        for idir=1:snq.n_dir
            figure(figID+2+idir)
            plot(npar.xf,reshape(psi(:,:,idir),npar.ndofs,1),'-','LineWidth',1);
            title(sprintf('%s angular flux %d, problem %d',forward_txt,idir,dat.pb_ID));
            xlabel('x'); ylabel('angular flux');
        end
    end
end

return
end
