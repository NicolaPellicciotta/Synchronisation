%data_dir='/u/homes/np451/Desktop/ependymal/sync/28.6.18/P3/';
data_dir='/media/np451/Seagate Expansion Drive1/DATA/13.10.18/P5/';
cd(data_dir);


%%% select a mask from the movie without external flow %%%%
d=dir('*P0*movie');
mo=moviereader(d(end).name);
fs=mo.read();
N_frames= size(fs,3);
s=std(double(fs(:,:,1:100)),[],3);
s=mat2gray(s);

prompt= 'how many cells do you want to check? should be a number  ';
Ncell = (input(prompt));
for nc=1:Ncell
BW{nc}=roipoly(s);
end
%fs_roi= fs(repmat(BW,[1,1,N_frames]));

d=dir('*movie');
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
    Vind=strfind(filename,'V');
    R{nc,jj}.Volt= str2num(filename(Vind+1));
    Offind=strfind(filename,'O');
    R{nc,jj}.Offset= str2num(filename(Offind+1));
    Pind=strfind(filename,'P');
    Endind=strfind(R{nc,jj}.filename,'.10Jun');
    R{nc,jj}.Period= str2num(filename(Pind+1:Endind-1));
    
    R{nc,jj}.FR = FR;
    temp_BW=BW{nc};
    fs_roi= fs(repmat(temp_BW,[1,1,N_frames]));
    fs_roi=reshape(fs_roi,[sum(temp_BW(:)),N_frames]);
    %roi=double(fs_roi)- mean(fs_roi,2);
    roi=double(fs_roi)- movmean(fs_roi,30,2);

%%%% periodogram needs time on the row and pixel on the coloun; 
    N_seg=1;
    window = hann(floor(N_frames)/N_seg);

    c_roi=roi';
    [pxx, fq]=pwelch(double(c_roi), window,50,size(c_roi,1),FR);
    R{nc,jj}.m_pxx= mean(pxx,2);
    R{nc,jj}.fq= fq;
    f_range=fq>7 & fq<30;
    [pks,locs,w,p] = findpeaks(R{nc,jj}.m_pxx(f_range),fq(f_range));
    [~,ind_sort]= sort(p);
    pks=pks(ind_sort); locs=locs(ind_sort);w=w(ind_sort);p=p(ind_sort);
    R{nc,jj}.pks=pks; R{nc,jj}.locs=locs; R{nc,jj}.w=w; R{nc,jj}.p=p;
    [~,where] = max(pks);
    R{nc,jj}.fp1 = locs(end);
    R{nc,jj}.fp2 = locs(end-1);
    
    end
 
end
    
    save('fft_results.mat','R','BW','s');
    
    %% select cell to analyse
    
    nc=1;
    cc=1;
    for kk=nc:size(R,1):numel(R);    
    Rc{cc}=R{kk};
    cc=cc+1;
    end
    
    %% repopulate results in the right range of frq
 
  fq_guess=15;

  for jj=1:(numel(Rc)-1)
    
    Pind=strfind(Rc{jj}.filename,'P');
    Endind=strfind(Rc{jj}.filename,'.10Jul');
    Rc{jj}.Period= str2num(Rc{jj}.filename(Pind+1:Endind-1));      
      
    fq=Rc{jj}.fq;
    f_range=Rc{jj}.fq> (fq_guess-8) & Rc{jj}.fq<(fq_guess+8);
    [pks,locs,w,p] = findpeaks(Rc{jj}.m_pxx(f_range),fq(f_range));
    [~,ind_sort]= sort(p);
    pks=pks(ind_sort); locs=locs(ind_sort);w=w(ind_sort);p=p(ind_sort);
    Rc{jj}.pks=pks; Rc{jj}.locs=locs; Rc{jj}.w=w; Rc{jj}.p=p;
    [~,where] = max(pks);
    Rc{jj}.fp1c = locs(end);
    if numel(locs)~=1
        Rc{jj}.fp2c = locs(end-1);
    else
        Rc{jj}.fp2c = nan;
    end
  end
 
 f_rest=[];
 cc=0;
 for jj=1:(numel(Rc)-1)
      if Rc{jj}.Period==0; cc=cc+1; f_rest(cc)=Rc{jj}.fp1c;end
 end
 
 F_rest=mean(f_rest);  
 %% plots Arnaud tongue 
 
 
 
 Periods=40:100;
 Volts=2:8;
 Amplitude= [200,600,1000,1500,2500,3000,4500]./200;  %%% in um
 Mat1= zeros([numel(Periods),numel(Volts)]);
 Mat2= zeros([numel(Periods),numel(Volts)]);
 figure(3);
 for jj=1:(numel(Rc)-1);

 Pind=strfind(Rc{jj}.filename,'P');     
 Endind=strfind(Rc{jj}.filename,'.10Jul');
 Rc{jj}.Period= str2num(Rc{jj}.filename(Pind+1:Endind-1));    
 v=Rc{jj}.Volt;
 p=Rc{jj}.Period;
 
 fp1=Rc{jj}.fp1c;   %%% frequency peak corrected in a smaller range
 fp2=Rc{jj}.fp2c;
 Mat1(p==Periods,v==Volts)= fp1-(1000/p);
 a=Amplitude(v==Volts);
 ff= f_rest(v==Volts);
 
 if abs(fp1-(1000/p))<= 0.2  & ( abs(fp2-ff )> 1.2) 
     plot(1000./p - ff ,a,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[0,0,0]);hold on;
 else
     plot(1000./p - ff ,a,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[1,1,1]);hold on;
 end
 
 end
 xlabel('frequency [Hz]');ylabel('flow velocity [um/s]');
 title(strcat('mean freq',num2str(F_rest)));
% ylim([1.5,9]);
% figure(2);imagesc(Volts,1000./Periods,Mat1);
% figure();imagesc(Mat1);
 
 
  %% plot different spectra
%close all 
Hz_toplot=[12:16]; 
Volt_toplot=[5];
Period_toplot= floor(1000./Hz_toplot);

ind=[];
cc=1;
for jj=1:(numel(Rc)-1) 
if any((Rc{jj}.Period)== Period_toplot) & any((Rc{jj}.Volt)== Volt_toplot); ind(cc)=jj;cc=cc+1;end
end

figure(1);
cleg=1;
 for jj=ind;
 figure(1);plot(Rc{jj}.fq(5:end),Rc{jj}.m_pxx(5:end),'LineWidth',2); hold on;
 %figure(2); plot(str2num(Hz{jj}),Rc{jj}.fp-str2num(Hz{jj}),'o'); hold on;
 leg{cleg}=num2str(1000/Rc{jj}.Period);cleg=cleg+1;
 end
 figure(1); legend(leg)
 xlabel('frequecy [Hz]'),ylabel('periodogram');
 xlim([Hz_toplot(1)-3,Hz_toplot(end)+3]);
 
 