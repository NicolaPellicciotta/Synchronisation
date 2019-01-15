%%% make some moviessss
%%

data_dir='/media/np451/Seagate Expansion Drive/ependymal/syn/29.7.18/P1';
cd(data_dir);

%% 

V_array=[2,4,5,6]
P_array=[0,58,0,58];
for kk=1:numel(V_array)
cd(data_dir);
V= V_array(kk);
P= P_array(kk);
lookfor= strcat('*V',num2str(V),'*P',num2str(P),'*movie');
d=dir(lookfor);
filename=d(1).name;
mo= moviereader(filename);
fs=mo.read();

fm=double(fs(:,:,1:200))-(movmean(fs(:,:,1:200),30,3));
fm=imadjustn(fm);
scroll_stack(fm);
cd(data_dir);
movie_fold= strcat('V',num2str(V),'_P',num2str(P),'_tiffs')
mkdir(movie_fold);
cd(movie_fold);

for hh=1:size(fm,3);
  
    imwrite(fm(:,:,hh),strcat(filename,'_N',num2str(hh),'.tiff'));
end
end