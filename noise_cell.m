%%%%% noise measuring 
clear all;

d=dir('*noise*movie');
mo=moviereader(d(1).name);
FR=mo.FrameRate;
fs=mo.read();
N_frames= size(fs,3);
load('BW.mat');
temp_BW=logical(BW{1});
fs_roi= fs(repmat(temp_BW,[1,1,N_frames]));
fs_roi=reshape(fs_roi,[sum(temp_BW(:)),N_frames]);
roi=double(fs_roi)- mean(fs_roi,2);
%%% find mean freq
window = hann(floor(N_frames));
window= repmat(window,[1,size(roi,1)])';
n= floor(N_frames/2);
if mod(n,2)==0; n= n-1;end
pxx= abs(fft(double(roi).*window,n,2)).^2;
m_pxx= mean(pxx(:,1:floor(n/2)),1);
fq= (0:(FR./n):(FR./2-FR./n));    
figure();plot(fq,m_pxx);
prompt= 'what is the frequency guess ';
fq_guess = (input(prompt));
fq_max= fq_guess + (7);
fq_min  =fq_guess- (7);

f_range=fq> (fq_min) & fq<(fq_max);
[pks,locs,w,p] = findpeaks(m_pxx(f_range),fq(f_range));%%%% find the peaks frequency in the selected freq range
baseline= mean((m_pxx(f_range)));
[~,ind_sort]= sort(pks);                                           %%%% sort peaks and get an index 
pks=pks(ind_sort); locs=locs(ind_sort);w=w(ind_sort);p=p(ind_sort);  %%% order all the variables with the same index     %%%%% load results in Rc
[~,where] = max(pks);         
frequency_T = locs(end);

beat_rep=50;
frames_beat=floor(FR/ frequency_T);
dt= frames_beat*beat_rep;
N_int=floor(N_frames/dt);
frequency=zeros([1,N_int]);
%%
t_array= 1:frames_beat:(N_frames-dt);
cc=1;
parfor tt=t_array
BW_tt= repmat(temp_BW,[1,1,dt]);
%fs_tt=fs(:,:,1+ (tt-1)*dt:tt*dt);
fs_tt=fs(:,:,tt:dt+tt);
fs_roi=fs_tt(BW_tt);
fs_roi=reshape(fs_roi,[sum(temp_BW(:)),dt]);
roi=double(fs_roi)- mean(fs_roi,2);
%%% find mean freq
window = hann(floor(dt));
window= repmat(window,[1,size(roi,1)])';
n= floor(N_frames/2);
if mod(n,2)==0; n= n-1;end
pxx= abs(fft(double(roi).*window,n,2)).^2;
m_pxx= mean(pxx(:,1:floor(n/2)),1);
fq= (0:(FR./n):(FR./2-FR./n));    

f_range=fq>fq_min & fq< fq_max;
[pks,locs,w,p] = findpeaks(m_pxx(f_range),fq(f_range));%%%% find the peaks frequency in the selected freq range
baseline= mean((m_pxx(f_range)));
[~,ind_sort]= sort(pks);                                           %%%% sort peaks and get an index 
pks=pks(ind_sort); locs=locs(ind_sort);w=w(ind_sort);p=p(ind_sort);  %%% order all the variables with the same index     %%%%% load results in Rc
[~,where] = max(pks);         
frequency(cc) = locs(end);
disp(tt);
%%%%%%
cc=cc+1;
    
end

noise= beat_rep*var(frequency)/(mean(frequency)^2); 
%%

