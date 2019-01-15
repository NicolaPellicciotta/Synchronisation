%data_dir='/u/homes/np451/Desktop/ependymal/sync/28.6.18/P3/';
close all;
clear all
data_dir='/media/np451/Seagate Expansion Drive/DATA/Cyto_D.23.10.18/P2_24hr/';
%data_dir='/media/np451/Seagate Expansion Drive/DATA/22.10.18/P2';
cd(data_dir);


%%% select a mask from the movie without external flow %%%%
d=dir('*P0*movie');
mo=moviereader(d(end).name);
fs=mo.read();
N_frames= size(fs,3);
s=std(double(fs(:,:,2:300)),[],3);
s=mat2gray(s);


fsbw= double(fs(:,:,2:300))- mean(fs(:,:,2:300),3);
k=0;
scroll_stack(fsbw)
%nc=1;
nc=11;
while k==0;
k = waitforbuttonpress;    
BW_temp{nc}=roipoly();
nc=nc+1;
end

cc=1;
for nc=1:numel(BW_temp)
if ~isempty(BW_temp{nc}); BW{cc}=BW_temp{nc}; cc=cc+1;end
end
%close all;
mkdir('plots')
%figure();
Reds= cat(3,s(:,:,1),0*s(:,:,1),0*s(:,:,1));  
    imshow(imadjust(s));
    BWT=zeros(size(s));
    for nc=1:numel(BW);
        BWT= BWT+BW{nc};
        [maskx,masky]=find(BW{nc});
        text(masky(1)-30,maskx(1)-30,num2str(nc),'Color','red','FontSize',24)
    end
    hold on
    h = imshow(Reds); % Save the handle; we'll need it later
    hold off        
    set(h, 'AlphaData', BWT); 
    saveas(gca,'plots/cell_analysed.jpg')


%%
Ncell=numel(BW);
% %
% prompt= 'how many cells do you want to check? should be a number  ';
% Ncell = (input(prompt));
% for nc=1:Ncell
% BW{nc}=roipoly(s);
% end
%fs_roi= fs(repmat(BW,[1,1,N_frames]));

d=dir('*_V*movie');
for jj=1:numel(d);
    disp(jj)
    filename=d(jj).name;
%%%%% load data and the roi
    mo=moviereader(filename);
    fs=mo.read();
    N_frames= size(fs,3);
    FR =mo.FrameRate;
    filename=d(jj).name;
    
    for nc=1:Ncell
    
    R{nc,jj}.filename=filename;
    Vind=strfind(R{jj}.filename,'V');
    Vend=strfind(R{jj}.filename,'_O');
    Rc{nc,jj}.Volt= str2num(R{jj}.filename(Vind+1:Vend-1));
    Offind=strfind(filename,'O');
    R{nc,jj}.Offset= str2num(filename(Offind+1));
    Pind=strfind(filename,'P');
    Endind=strfind(R{nc,jj}.filename,'.29Jul');
    R{nc,jj}.Period= str2num(filename(Pind+1:Endind-1));
    
    R{nc,jj}.FR = FR;
    temp_BW=BW{nc};
    fs_roi= fs(repmat(temp_BW,[1,1,N_frames]));
    fs_roi=reshape(fs_roi,[sum(temp_BW(:)),N_frames]);
    %roi=double(fs_roi)- mean(fs_roi,2);
    roi=double(fs_roi)- movmean(fs_roi,30,2);

%%%% periodogram needs time on the row and pixel on the coloun; 
    window = hann(floor(N_frames));
    window= repmat(window,[1,size(roi,1)])';
    n= floor(N_frames/2);
    if mod(n,2)==0; n= n-1;end
    pxx= abs(fft(double(roi).*window,n,2)).^2;
    R{nc,jj}.m_pxx= mean(pxx(:,1:floor(n/2)),1);
    fq= (0:(FR./n):(FR./2-FR./n));    
    R{nc,jj}.fq= fq;
%     f_range=fq>5 & fq<40;
%     [pks,locs,w,p] = findpeaks(R{nc,jj}.m_pxx(f_range),fq(f_range));
%     [~,ind_sort]= sort(p);
%     pks=pks(ind_sort); locs=locs(ind_sort);w=w(ind_sort);p=p(ind_sort);
%     R{nc,jj}.pks=pks; R{nc,jj}.locs=locs; R{nc,jj}.w=w; R{nc,jj}.p=p;
%     [~,where] = max(pks);
%     R{nc,jj}.fp1 = locs(end);
%     R{nc,jj}.fp2 = locs(end-1);
%     
    end
 
end
    
    
    save('fft_results_fft.mat','R','BW','s');
    