%%%%%% with this I am trying to rewrite the code for image correlation
%%%%%% using a rectangular window... lets see
 filename='60X_P3_00.00_.00Hz_7V_9O.05Jun2018_19.30.30.movie';


  filename='60X_P3_20.00Hz_7V_9O.05Jun2018_19.17.48.movie';

  
  mo=moviereader(filename);
  FR=mo.FrameRate;
  fs=mo.read();
  N_frames= size(fs,3);
  s=std(double(fs(:,:,1:100)),[],3);
  imagesc(s);
  rect=imrect;
  wait(rect);
  xy = rect.getPosition;
  [rows,cols]=rect2sub(xy);
  
  roi= fs(rows(1):rows(end),cols(1):cols(end),:);
  
  Iqtau= DDM_correlation(roi);
  
  temp=Iqtau(1,1:floor(end/2));
  [pks,locs,w,p] = findpeaks(smooth(smooth(temp)));
  p_range= p
  N_seg=3;
  window = hann(floor(numel(temp))/N_seg);

 [pxx, fq]=pwelch(temp',window,50,numel(temp),FR);

  figure();plot(pxx);
  xlabel('frequecy [Hz]'),ylabel('periodogram');