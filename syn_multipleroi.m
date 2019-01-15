%data_dir='/u/homes/np451/Desktop/ependymal/sync/28.6.18/P3/';
data_dir='/media/np451/Seagate Expansion Drive/ependymal/syn/10.7.18/P4';
cd(data_dir);

%%% select a mask from the movie without external flow %%%%
d=dir('*P0*movie');
mo=moviereader(d(1).name);
fs=mo.read();
N_frames= size(fs,3);
s=std(double(fs(:,:,1:100)),[],3);
s=mat2gray(s);

BW=roipoly(s); 
%fs_roi= fs(repmat(BW,[1,1,N_frames]));

d=dir('*movie');
for jj=1:numel(d);
    filename=d(jj).name;
    R{jj}.filename=filename;
    Vind=strfind(filename,'V');
    R{jj}.Volt= str2num(filename(Vind+1));
    Offind=strfind(filename,'O');
    R{jj}.Offset= str2num(filename(Offind+1));
    Pind=strfind(filename,'P');
    Endind=strfind(R{jj}.filename,'.28Jun');
    R{jj}.Period= str2num(filename(Pind+1:Endind-1));
    
    
%%%%% load data and the roi
    mo=moviereader(filename);
    fs=mo.read();
    N_frames= size(fs,3);
    FR =mo.FrameRate;
    R{jj}.FR = FR;
    fs_roi= fs(repmat(BW,[1,1,N_frames]));
    fs_roi=reshape(fs_roi,[sum(BW(:)),N_frames]);
    %roi=double(fs_roi)- mean(fs_roi,2);
    roi=double(fs_roi)- movmean(fs_roi,3,2);

%%%% periodogram needs time on the row and pixel on the coloun; 
    N_seg=3;
    window = hann(floor(N_frames)/N_seg);

    c_roi=roi';
    [pxx, fq]=pwelch(double(c_roi), window,50,size(c_roi,1),FR);
    R{jj}.m_pxx= mean(pxx,2);
    R{jj}.fq= fq;
    f_range=fq>7 & fq<40;
    [pks,locs,w,p] = findpeaks(R{jj}.m_pxx(f_range),fq(f_range));
    [~,ind_sort]= sort(p);
    pks=pks(ind_sort); locs=locs(ind_sort);w=w(ind_sort);p=p(ind_sort);
    R{jj}.pks=pks; R{jj}.locs=locs; R{jj}.w=w; R{jj}.p=p;
    [~,where] = max(pks);
    R{jj}.fp1 = locs(end);
    R{jj}.fp2 = locs(end-1);
    
    
 
end
    
    save('fft_results.mat','R','BW','s');
 %% repopulate results in the right range of frq
 
 fq_guess=23;

  for jj=1:(numel(R)-1)
    fq=R{jj}.fq;  
    f_range=R{jj}.fq>(fq_guess-8) & R{jj}.fq<(fq_guess+8);
    [pks,locs,w,p] = findpeaks(R{jj}.m_pxx(f_range),fq(f_range));
    [~,ind_sort]= sort(p);
    pks=pks(ind_sort); locs=locs(ind_sort);w=w(ind_sort);p=p(ind_sort);
    R{jj}.pks=pks; R{jj}.locs=locs; R{jj}.w=w; R{jj}.p=p;
    [~,where] = max(pks);
    R{jj}.fp1c = locs(end);
    if numel(locs)~=1
        R{jj}.fp2c = locs(end-1);
    else
        R{jj}.fp2c = nan;
    end
  end
 
 f_rest=[];
 cc=0;
 for jj=1:(numel(R)-1)
      if R{jj}.Period==0; cc=cc+1; f_rest(cc)=R{jj}.fp1c;end
 end
 
 F_rest=mean(f_rest);  
 %% plots Arnaud tongue 
 
 
 
 Periods=36:70;
 Volts=2:8;
 Amplitude= 2:8;  %%% in um
 Mat1= zeros([numel(Periods),numel(Volts)]);
 Mat2= zeros([numel(Periods),numel(Volts)]);
 figure(4);
 for jj=1:(numel(R)-1);

 Pind=strfind(R{jj}.filename,'P');     
 Endind=strfind(R{jj}.filename,'.10Jul');
 R{jj}.Period= str2num(R{jj}.filename(Pind+1:Endind-1));    
 v=R{jj}.Volt;
 p=R{jj}.Period;
 
 fp1=R{jj}.fp1c;   %%% frequency peak corrected in a smaller range
 fp2=R{jj}.fp2c;
 Mat1(p==Periods,v==Volts)= fp1-(1000/p);
 a=Amplitude(v==Volts);
 if abs(fp1-(1000/p))<= 0.15  & ( abs(fp2-F_rest )> 1) 
     plot(1000./p,a*1000./p,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[0,0,0]);hold on;
 else
     plot(1000./p,a*1000./p,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[1,1,1]);hold on;
 end
 
 end
 xlabel('frequency [Hz]');ylabel('flow velocity [um/s]');
% ylim([1.5,9]);
% figure(2);imagesc(Volts,1000./Periods,Mat1);
% figure();imagesc(Mat1);
 
 
  %% plot different spectra
close all 
Hz_toplot=[18,20,21,22,23]; 
Volt_toplot=[8];
Period_toplot= floor(1000./Hz_toplot);

ind=[];
cc=1;
for jj=1:(numel(R)-1) 
if any((R{jj}.Period)== Period_toplot) & any((R{jj}.Volt)== Volt_toplot); ind(cc)=jj;cc=cc+1;end
end

figure(1);
cleg=1;
 for jj=ind;
 figure(1);plot(R{jj}.fq(5:end),R{jj}.m_pxx(5:end),'LineWidth',2);hold on;
 %figure(2); plot(str2num(Hz{jj}),R{jj}.fp-str2num(Hz{jj}),'o'); hold on;
 leg{cleg}=num2str(1000/R{jj}.Period);cleg=cleg+1;
 end
 figure(1); legend(leg)
 xlabel('frequecy [Hz]'),ylabel('periodogram');
 xlim([18,27])
 
 
 
 
 %%
    %%%% periodogram for BG 
BG=roipoly(s);
fs_bg= fs(repmat(BG,[1,1,N_frames]));  
roi_bg=reshape(fs_bg,[sum(BG(:)),N_frames]);
roi_bg=double(roi_bg)- mean(roi_bg,2);
 
roi_bg=roi_bg';
[pxx, fq]=pwelch(double(roi_bg), window,50,size(roi_bg,1),FR);
R{1}.bg_pxx= mean(pxx,2);
 
R{1}.pxx= R{1}.m_pxx - R{1}.bg_pxx;   

f_range=fq>15 & fq<35
[pks,locs] = findpeaks(R{1}.pxx(f_range),fq(f_range));
[~,where] = max(pks);
R{1}.fp = locs(where); 

 
 
 
 for jj=2:numel(filename)
 
  mo=moviereader(filename{jj});
  FR=mo.FrameRate; R{jj}.FR = FR;
  fs=mo.read();
  N_frames= size(fs,3);
 
  fs_roi= fs(repmat(BW,[1,1,N_frames]));
  roi=reshape(fs_roi,[sum(BW(:)),N_frames]);
  roi=double(roi)- mean(roi,2);


%%%% periodogram needs time on the row and pixel on the coloun; 
c_roi=roi';
N_seg=1;
window = hann(floor(N_frames)/N_seg);

[pxx, fq]=pwelch(double(c_roi), window,50,size(c_roi,1),FR);
 R{jj}.fq= fq;
 R{jj}.m_pxx= mean(pxx,2);
          
%%%% background %%%%%
fs_bg= fs(repmat(BG,[1,1,N_frames]));  
roi_bg=reshape(fs_bg,[sum(BG(:)),N_frames]);
roi_bg=double(roi_bg)- mean(roi_bg,2);
 
roi_bg=roi_bg';
[pxx, fq]=pwelch(double(roi_bg), window,50,size(roi_bg,1),FR);
R{jj}.bg_pxx= mean(pxx,2);
 
R{jj}.pxx= R{jj}.m_pxx - R{jj}.bg_pxx; 
 
f_range=fq>15 & fq<35;
[pks,locs] = findpeaks(R{jj}.pxx(f_range),fq(f_range));
[~,where] = max(pks);
R{jj}.fp = locs(where); 
 
 
 end
 
 %%
 
Hz_toplot=[21];  

ind=[];
cc=1;
for jj=1:numel(R) 
if any(floor(str2num(Hz{jj}))== Hz_toplot); ind(cc)=jj;cc=cc+1;end
end

figure(1);figure(2);
%ind=2:numel(d);
 for jj=ind;
 figure(1);plot(R{jj}.fq(5:end),R{jj}.m_pxx(5:end),'LineWidth',2); hold on;
 figure(2); plot(str2num(Hz{jj}),R{jj}.fp-str2num(Hz{jj}),'o'); hold on;
 end
 figure(1); legend(Hz{ind})
 xlabel('frequecy [Hz]'),ylabel('periodogram');
 xlim([18,30])
 
 
 