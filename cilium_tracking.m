


%%%%% script to track center of drag of a defined cilium 

%% 
clear all
data_dir='/media/np451/Seagate Backup Plus Drive/DATA/' ;
cd(data_dir)
d=dir('*11.58*movie');
mo=moviereader(d(1).name);
fs=mo.read();
%%
%%%select roi
s= std(double(fs(:,:,1:300)),[],3);
fs_pic= double(fs(:,:,2:300))- movmean(fs(:,:,2:300),300,3);
scroll_stack(fs_pic);
rect_out= getrect();
close all
[rows,cols]= rect2sub(rect_out);
fs_roi=fs(rows(1):rows(end),cols(1):cols(end),:);

fsbw= double(fs_roi(:,:,2:300))- movmean(fs_roi(:,:,2:300),300,3);
    for kk=1:size(fsbw,3)
    fsbw2(:,:,kk)= wiener2(fsbw(:,:,kk),[4,4]);
    end
   % fsbw2=imadjustn(fsbw2);
    k=0;
    scroll_stack(fsbw2)
    cc=1;
    while k==0
    k = waitforbuttonpress;
    [x,y] = getpts;
    xx(cc)=x;
    yy(cc)=y;
    cc=cc+1;
    end

%% trajectory fitting    
fps=mo.FrameRate;
px2mu= 0.14/40*60;

xdot=smooth(diff(xx))*fps*px2mu;  
ydot=smooth(diff(xx))*fps*px2mu; 
xmid= 0.5*(xx(1:end-1)+xx(2:end))*px2mu;
ymid= 0.5*(yy(1:end-1)+yy(2:end))*px2mu;


plot(xmid,xdot,'o');
xlabel('$x[\mu m]$','interpreter','latex');
ylabel('$\dot{x} [\mu m/s]$','interpreter','latex');

ydotm=ydot(ydot<0);
ymidm= ymid(ydot<0);
ydotu=ydot(ydot<0);
ymidu= ymid(ydot<0);

xdotm=xdot(xdot<0);
xmidm= xmid(xdot<0);
xdotu=xdot(xdot>0);
xmidu= xmid(xdot>0);

pxm=polyfit(xmidm,xdotm',3);
pxu=polyfit(xmidu,xdotu',3);
hold on;
xx_array=linspace(min(xmid(:)),max(xmid(:)));
plot(xx_array,polyval(pxm,xx_array),'k-','MarkerSize',7,'LineWidth',1);
plot(xx_array,polyval(pxu,xx_array),'k-','MarkerSize',7,'LineWidth',1);

%%
clear fs;
save('cilium_tracking.mat')
save('cilia_tracking_data_P11.mat','pxu','pxm','xmidu','xmidm','xdotm','ydotm');