function [Meta] = spatial_noise_fft_syn(f_guess,box_size,movie,Sp_rest)
% 

if nargin < 3 || isempty(movie)
    movie='P0';
end

% d=dir('*noise*movie');
% Meta.movie='noise';
% 
% if isempty(d); d=dir('*P0*movie');Noise.movie='P0';end
% mo=moviereader(d(1).name);
d=dir(strcat('*',movie,'*movie'));
mo=moviereader(d(1).name);
Meta.mo=mo;

fq_max=f_guess+7;fq_min=f_guess-7;


fs=mo.read();
N_frames=size(fs,3);
FR=mo.FrameRate;
s= std(double(fs(:,:,1:300)),[],3);
[level,EM] = graythresh(s);
s_bin = imbinarize(s,level);   %%%%% standard deviation of the movie binarised between 1 and 0

BW=imbinarize(s,level);
rect_out= Sp_rest.rect_out;
[rows,cols]= rect2sub(rect_out);
fs_roi= fs(rows(1):rows(end),cols(1):cols(end),:);
fsm_roi=double(fs_roi)-movmean(fs_roi,100,3);
row_dim= size(fsm_roi,1);
col_dim= size(fsm_roi,2);

s_roi=s(rows(1):rows(end),cols(1):cols(end));
mask_tot= zeros(size(s_roi));
[level,EM] = graythresh(s_roi);


s_bin = Sp_rest.s_bin;
%box_size=8;   %%%% setting the box size 
Nbox=floor((row_dim*col_dim)/box_size^2);  %%%% number of boxes based on the total area
nboxes_row=floor(row_dim/box_size);
nboxes_col=floor(col_dim/box_size);
[X,Y]=meshgrid(1: box_size: nboxes_row*box_size, 1:box_size:nboxes_col*box_size );

F= zeros(size(X)-1);
P= zeros(size(X)-1);


for xx=1:(size(X,1)-1);
    for yy=1:(size(X,2)-1);    

        std_value= s_bin(X(xx,yy):X(xx+1,yy+1),Y(xx,yy):Y(xx+1,yy+1));
        if mean(std_value(:))> 0.7; good=1;
            mask_tot(X(xx,yy):X(xx+1,yy+1),Y(xx,yy):Y(xx+1,yy+1))=1;
        else good=0;
        end
        Meta.Box(xx,yy).good=good;
        
        roi=fsm_roi(X(xx,yy):X(xx+1,yy+1),Y(xx,yy):Y(xx+1,yy+1),:);
        roi=reshape(roi,[size(roi,1)*size(roi,2),N_frames]);
        
        Meta.Box(xx,yy).vec_pos=[X(xx,yy),Y(xx,yy),X(xx+1,yy+1),Y(xx+1,yy+1)];
        


        %%%%%% calculate fft on the roi
        window = hann(floor(N_frames));
        window= repmat(window,[1,size(roi,1)])';
        %n= floor(N_frames/2);
        n=N_frames;
        if mod(n,2)==0; n= n-1;end
        
        f_spectrum= fft(double(roi).*window,n,2);
        phase=angle(f_spectrum);
        pxx= abs(f_spectrum).^2;
        m_pxx= mean(pxx(:,1:floor(n/2)),1);
        fq= (0:(FR./n):(FR./2-FR./n));    
        Meta.fq=fq;
%%%        Meta.Box(xx,yy).m_pxx=m_pxx;   %%% commented 
%%%        Meta.Box(xx,yy).f_spectrum=f_spectrum;
        
        %%%%% calculate the peak in frequency
        
        f_range= fq> (fq_min) & fq<(fq_max);
        baseline= min((m_pxx(f_range)));
        [pks,locs,w,p] = findpeaks((m_pxx(f_range)),fq(f_range));%%%% find the peaks frequency in the selected freq range
        [~,ind_sort]= sort(pks);                                           %%%% sort peaks and get an index 
        pks=pks(ind_sort); locs=locs(ind_sort);w=w(ind_sort);p=p(ind_sort);  %%% order all the variables with the same index
        Meta.Box(xx,yy).pks=pks; Meta.Box(xx,yy).locs=locs; Meta.Box(xx,yy).w=w; Meta.Box(xx,yy).p=p;     %%%%% load results in Rc
        Meta.Box(xx,yy).baseline=baseline;
        [~,where] = max(pks);       
        if numel(locs)==0 | good==0; 
            Meta.Box(xx,yy).fp1 =nan;   %%% fp1c is the frequency with higher peak. 
            F(xx,yy)=nan;
            P(xx,yy)=nan;
        else
            F(xx,yy)=locs(end);
            P(xx,yy)=phase(fq==locs(end)); 
        end
        
        Meta.Box(xx,yy).fp1 = F(xx,yy);
        Meta.Box(xx,yy).phase = P(xx,yy);
end
end

Meta.F=F;Meta.P=P; Meta.s=s;Meta.s_bin=s_bin;Meta.s_roi=s_roi;Meta.N_frames=N_frames;
Meta.box_size=box_size;
Meta.rect_out=rect_out;
Meta.median_freq= nanmedian(F(:));
Meta.std_freq= nanstd(F(:));
Meta.X=X;
Meta.Y=Y;
Meta.mask_tot=mask_tot;
end