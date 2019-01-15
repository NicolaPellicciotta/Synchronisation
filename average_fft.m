    function [fq,m_pxx] = average_fft(fs,mask,FR)
    
    N_frames= size(fs,3);
    fs_roi= fs(repmat(mask,[1,1,N_frames]));
    fs_roi=reshape(fs_roi,[sum(mask(:)),N_frames]);
    %roi=double(fs_roi)- mean(fs_roi,2);
    roi=double(fs_roi)- movmean(fs_roi,60,2);

%%%% periodogram needs time on the row and pixel on the coloun; 
    window = hann(floor(N_frames));
    window= repmat(window,[1,size(roi,1)])';
    n= floor(N_frames/2);
    if mod(n,2)==0; n= n-1;end
    pxx= abs(fft(double(roi).*window,n,2)).^2;
    m_pxx= mean(pxx(:,1:floor(n/2)),1);
    fq= (0:(FR./n):(FR./2-FR./n));    
    
    end