data_dir='/u/homes/np451/Desktop/ependymal/sync/19.6.18/P4/'
cd(data_dir);

d=[dir('*0V*movie');dir('*7V*movie')];

for jj=1:numel(d);
    filename{jj}=d(jj).name;
    Vind=strfind('V_',d(jj).name);
    
    Hz{jj}=d(jj).name(8:12); 
end

% filename={'60X_P2_00.00Hz_0V_0O.05Jun2018_18.13.43.movie',...
%     '60X_P2_21.22Hz_7V_9O.05Jun2018_18.46.26.movie',...
%     '60X_P2_21.28.00Hz_4V_5O.05Jun2018_18.15.36.movie',...
%      '60X_P2_21.74.00Hz_4V_5O.05Jun2018_18.16.19.movie',...
%      '60X_P2_22.22.00Hz_4V_5O.05Jun2018_18.17.07.movie',...
%      '60X_P2_22.73.00Hz_4V_5O.05Jun2018_18.17.43.movie',...
%      '60X_P2_23.26.00Hz_4V_5O.05Jun2018_18.18.16.movie',...
%      '60X_P2_23.81.00Hz_4V_5O.05Jun2018_18.18.58.movie'};
%   
%   
  
  mo=moviereader(filename{1});
  FR =mo.FrameRate;
  R{1}.FR = FR;
  fs=mo.read();
  N_frames= size(fs,3);
  s=std(double(fs(:,:,1:100)),[],3);
  s=mat2gray(s);


  BW=roipoly(s); 
  fs_roi= fs(repmat(BW,[1,1,N_frames]));
  roi=reshape(fs_roi,[sum(BW(:)),N_frames]);
  roi=double(roi)- mean(roi,2);
  

%%%% periodogram needs time on the row and pixel on the coloun; 
N_seg=3;
window = hann(floor(N_frames)/N_seg);

c_roi=roi';
[pxx, fq]=pwelch(double(c_roi), window,50,size(c_roi,1),FR);
 R{1}.m_pxx= mean(pxx,2);
 R{1}.fq= fq;

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
for jj=1:numel(d) 
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
 
 
 