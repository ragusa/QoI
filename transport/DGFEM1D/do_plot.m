function do_plot(phi,psi,x,spa_disc)

if(spa_disc<2)
    phi(:,2)=phi(:,1);
    psi(:,:,2)=psi(:,:,1);
end

switch spa_disc
    case 0
        tit='SD';
    case 1
        tit='DD';
    case 2
        tit='LD';
end
sn=length(psi(1,:,1));
    
figure(1)
for i=1:length(x)-1,
    plot(x(i:i+1),phi(i,1:2),'-','LineWidth',2)
    hold on
end
filename=sprintf('scal_%s_sn%i.png',tit,sn);
print('-dpng',filename);

figure(2)
xm=diff(x)/2+x(1:end-1);
plot(xm,sum(phi,2),'LineWidth',2)
filename=sprintf('scal2_%s_sn%i.png',tit,sn);
print('-dpng',filename);

for j=1:sn
    figure(10+j)
    for i=1:length(x)-1,
        plot(x(i:i+1),shiftdim(psi(i,j,1:2)),'-','LineWidth',2)
        hold on
    end
    filename=sprintf('angu%s_%i_sn%i.png',tit,j,sn);
    print('-dpng',filename);
end

for j=1:sn
    figure(40+j)
    plot(xm,sum(psi(:,j,:),3),'-','LineWidth',2);
    filename=sprintf('ang%s_%i_sn%i.png',tit,j,sn);
    print('-dpng',filename);
end

return
end
