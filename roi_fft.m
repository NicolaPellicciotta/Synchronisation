  


 filename='60X_P2_00.00Hz_0V_0O.05Jun2018_16.58.57.movie';

  
  
  
  
  mo=moviereader(filename);
  FR=mo.FrameRate;
  fs=mo.read();
  N_frames= size(fs,3);
  s=std(double(fs(:,:,1:100)),[],3);
  %%
  %%%%% background noise %%%%
  BW=logical(ones(size(s)));
   fs_roi= fs(repmat(BW,[1,1,N_frames]));
  roi=reshape(fs_roi,[sum(BW(:)),N_frames]);
  roi=double(roi)- mean(roi,2);


 % pxx is a matrix with the spectrum of each pixel as a column
%c_roi = AutoCorr(double(roi'),floor(N_frames/2));

%%%% periodogram needs time on the row and pixel on the coloun; 

c_roi=roi';

clear bg_pxx;
bg_pxx=zeros([(floor(N_frames/2) +1),1]);
Npx=1000;
N_rep= floor(size(c_roi,2)/Npx);
for jj=1:N_rep
    ind_l=(jj-1)*Npx +1;ind_r=(jj)*Npx; 
    temp_roi=c_roi(:,ind_l:ind_r);
    N_seg=3;
    window = hann(floor(N_frames)/N_seg);

    [pxx, fq]=pwelch(double(temp_roi), window,50,size(temp_roi,1),FR);
%    [pxx, fq] = periodogram(double(temp_roi),window,size(temp_roi,1),FR);
%     [pwxx, frequencies] = pwelch(temp_fs,wwindow,[],numel(wwindow),FR);

  % average points of the spectra with same box index and same frequency
%  boxavg_pxx = accumarray({freq_ind(:),temp_ind_mat(:)},pxx(:)) ./ (bsz*bsz);
    bg_pxx= bg_pxx+ mean(pxx,2);
end
bg_pxx=bg_pxx./N_rep;
figure(); plot(fq(5:end),bg_pxx(5:end)); xlabel('[Hz]'),ylabel('fft');
  
  
%%  
  %%%% on single cell %%%%
  
  s=std(double(fs(:,:,1:100)),[],3);
  BW=roipoly(mat2gray(s));
  fs_roi= fs(repmat(BW,[1,1,N_frames]));
  roi=reshape(fs_roi,[sum(BW(:)),N_frames]);
  roi=double(roi)- mean(roi,2);


 % pxx is a matrix with the spectrum of each pixel as a column
%c_roi = AutoCorr(double(roi'),floor(N_frames/2));

%%%% periodogram needs time on the row and pixel on the coloun; 
c_roi=roi';
N_seg=3;
window = hann(floor(N_frames)/N_seg);

[pxx, fq]=pwelch(double(c_roi), window,50,size(c_roi,1),FR);
%[pxx, fq] = periodogram(double(c_roi),window,size(c_roi,1),FR);
%     [pwxx, frequencies] = pwelch(temp_fs,wwindow,[],numel(wwindow),FR);

  % average points of the spectra with same box index and same frequency
%  boxavg_pxx = accumarray({freq_ind(:),temp_ind_mat(:)},pxx(:)) ./ (bsz*bsz);
  m_pxx= mean(pxx,2);
  figure(); plot(fq(5:end),m_pxx(5:end)); xlabel('[Hz]'),ylabel('fft');
  %%%-bg_pxx(5:end)%%%
 