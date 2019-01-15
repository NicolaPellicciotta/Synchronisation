 %% plots Arnaud tongue 
%data_dir='/media/np451/Seagate Expansion Drive1/23.10.18/P1';
clear all;
close all;
Cal=[];
%data_dir='/media/np451/Seagate Expansion Drive/29.10.18/P';
data_dir='/media/np451/Seagate Backup Plus Drive/DATA/Cyto-D/Cyto_D.29.10.18/7hr/P10';
data_dir='/media/np451/Seagate Backup Plus Drive/DATA/Cyto-D/Cyto_D.29.10.18/24hr/P7';
%data_dir='/media/np451/Seagate Backup Plus Drive/DATA/1.11.18/P6';
data_dir='/media/np451/Seagate Expansion Drive/DATA/19.10.18/P5'

data_dir='/media/np451/Seagate Backup Plus Drive/DATA/26.11.18/P28';
cd(data_dir);

load('fft_results_fft.mat');
mkdir('plots');


%Rcen= R;
%clear R
%for i=53:numel(Rcen); R{i-52}=Rcen{i};end;
 
     
    %%% select cell to analyse
    Reds= cat(3,s(:,:,1),0*s(:,:,1),0*s(:,:,1));  
    imshow(imadjust(s));
    BWT=zeros(size(s));
    for nc=1:(size(R,1));
        BWT= BWT+BW{nc};
        [maskx,masky]=find(BW{nc});
        text(masky(1)-30,maskx(1)-30,num2str(nc),'Color','red','FontSize',24)
    end
    hold on
    h = imshow(Reds); % Save the handle; we'll need it later
    hold off        
    set(h, 'AlphaData', BWT); 
    saveas(gca,'plots/cell_analysed.jpg')
    
 %   %%%
   close all; 
   clear Rc;
    nc=7;
    cc=1;
    for kk= nc:size(R,1):numel(R);    
    Rc{cc}=R{kk};
    cc=cc+1;
    end
%  figure(5);imshow(BW{nc});   
%%%% get all the information about Period and Voltage right and then put in
%%%% the array 'ind' all the inex when P==0 and when Voltage Gain V<4 so
%%%% then to plot it
%%%
ind=[];
cc=1;
for jj=1:(numel(Rc)-1) 
        
    Pind=strfind(Rc{jj}.filename,'P');
%    Endind=strfind(Rc{jj}.filename,'.13Oct');
%    Rc{jj}.Period= str2num(Rc{jj}.filename(Pind+1:Endind-1));
    Rc{jj}.Period= str2num(Rc{jj}.filename(Pind+1:(end-25)));
    Vind=strfind(Rc{jj}.filename,'V');
    Vend=strfind(Rc{jj}.filename,'_O');
    Rc{jj}.Volt= str2num(Rc{jj}.filename(Vind+1:Vend-1));
%if (Rc{jj}.Period==0 & Rc{jj}.Volt <4); ind(cc)=jj;cc=cc+1;end    
if (Rc{jj}.Period==0 ); ind(cc)=jj;cc=cc+1;end    
    Rc{jj}.Freq= 1000./Rc{jj}.Period;     
end

%%%% plot all the spectra when P=0 e V<4 and decide the frequexy at rest of
%%%% the ciliated cell

figure(1);
cleg=1;
clear leg;
 for jj=ind;
 figure(1);plot(Rc{jj}.fq(5:end),smooth(Rc{jj}.m_pxx(5:end)),'LineWidth',2); hold on;
 %figure(2); plot(str2num(Hz{jj}),Rc{jj}.fp-str2num(Hz{jj}),'o'); hold on;
 leg{cleg}=num2str(Rc{jj}.Volt);cleg=cleg+1;
 end
 figure(1); legend(leg);
 xlabel('frequecy [Hz]'),ylabel('periodogram');
 xlim([12,35]);
  
prompt= 'what is the frequency guess ';
fq_guess = (input(prompt));
%prompt= 'what is estimate number of cilia? ';
%N_cilia = (input(prompt));
%Cell(nc).N_cilia=N_cilia;

%%%% start to add result to the class Cell, this will be useful later for
%%%% analysing multiple data

Cell(nc).s=s;
Cell(nc).BW=BW; 

%%% populate results in the right range of frq, use the frequecy at rest
%%% of the cilium to evaluate a good range where to look for
%%% synchronisation
  fq_max= fq_guess + (7);
  fq_min  =fq_guess- (7);

  for jj=1:(numel(Rc)-1)

    fq=Rc{jj}.fq;
    f_range=Rc{jj}.fq> (fq_min) & Rc{jj}.fq<(fq_max);
    baseline= min((Rc{jj}.m_pxx(f_range)));
    [pks,locs,w,p] = findpeaks((Rc{jj}.m_pxx(f_range)),fq(f_range));%%%% find the peaks frequency in the selected freq range
    baseline= 0%min((Rc{jj}.m_pxx(f_range)));
    [~,ind_sort]= sort(pks);                                           %%%% sort peaks and get an index 
    pks=pks(ind_sort); locs=locs(ind_sort);w=w(ind_sort);p=p(ind_sort);  %%% order all the variables with the same index
    Rc{jj}.pks=pks; Rc{jj}.locs=locs; Rc{jj}.w=w; Rc{jj}.p=p;     %%%%% load results in Rc
    Rc{jj}.baseline=baseline;
    [~,where] = max(pks);                                       
      
    
    if numel(locs)==0; Rc{jj}.fp1c =nan;   %%% fp1c is the frequency with higher peak. 
        else Rc{jj}.fp1 = locs(end);
    end
    if numel(locs)>1
        Rc{jj}.fp2 = locs(end-1);     %%% fp2c is the second frequency high second higher peak, it is a nan if it does not exixt
        Rc{jj}.pk2 = pks(end-1);
    else
        Rc{jj}.fp2 = 0; 
        Rc{jj}.pk2=0;
    end
  end
 
%%% f_rest is an array with all the frequency at rest of the cell, f_peak the corresponding peak value of the fft  
%%% F_rest is the mean frequency and F_peak the average peak
 f_rest=[];
 volt_exp=[];
 f_peak=[];
 cc=0;
 for jj=1:(numel(Rc)-1)
      if Rc{jj}.Period==0; cc=cc+1; f_rest(cc)=Rc{jj}.fp1;f_peak(cc)= Rc{jj}.pks(end); volt_exp(cc)= Rc{jj}.Volt;
      end
 end
 
 F_rest=mean(f_rest);  
 F_peak= mean(f_peak);
 F_pstd= std(f_peak);
 
 
 %%% load some data in Cell
 Cell(nc).F_rest=mean(f_rest);
 Cell(nc).f_rest=f_rest;
 Cell(nc).f_peak=f_peak;
 Cell(nc).F_peak= mean(f_peak);
 Cell(nc).F_pstd= std(f_peak);
 Cell(nc).f_range=f_range; 
 
 %%%% load the calibaration matrix 
 
%cd('/u/homes/np451/Documents/MATLAB/frequency/synchronization/');
Cal= load('calibration_matrix.mat');
%Cal=[];
cd(data_dir);
 
 %%% 
 Periods=51:120;
 Volts=2:0.5:6;
 Mat1= zeros([numel(Periods),numel(Volts)]);
 Mat2= zeros([numel(Periods),numel(Volts)]);

 
 for jj=1:(numel(Rc)-1);

 

 v=Rc{jj}.Volt;
 p=Rc{jj}.Period;    
 
 %%%%%% find the flow using the Calibration matrix
 if isempty(Cal) == 1
        flow= Rc{jj}.Volt;
 else   flow = Cal.Flow_int(Cal.volt_int==v , Cal.fq_int== 1000/p );  
 end

 fp1=Rc{jj}.fp1; %%% frequency peak corrected in a smaller range
 fpw=Rc{jj}.w(end);
 fp2=Rc{jj}.fp2;
 pk2= Rc{jj}.pk2;
 baseline= Rc{jj}.baseline
 Mat1(p==Periods,v==Volts)= fp1-(1000/p);
 %a=Amplitude(v==Volts);
 %ff= f_rest(v==Volts);
 ff= f_rest(1);
 
 
 threshold=0.66;
 Cell(nc).threshold=threshold;

 if p~=0%& p<=100
    peak1_cond=  abs(fp1-(1000/p))<= 0.5;% & fpw<1.5;  %%%(0.2 freq resolution)
%    peak2_cond=   abs(fp2-ff )> 3 |  abs( (F_peak-pk2)/F_peak )> threshold ;  %% 20% of the inital peak 
%
    peak2_cond=   abs(fp2-ff )> 3 | abs( (F_peak-pk2)/F_peak )> threshold ;  %% 20% of the inital peak 

    Rc{jj}.peak1_cond=peak1_cond;
    Rc{jj}.peak2_cond=peak2_cond;
 
    
    if  peak1_cond == 1 & fpw<2;
        if abs(fp2-ff )> 3; Sm=0; Rc{jj}.df= fp1-ff; 
        else Sm = abs( (pk2-baseline)/(F_peak-baseline)); 
            if Sm>0.40;
                Sm=1; Rc{jj}.df= fp2-ff; 
            else Sm=0; Rc{jj}.df= fp1-ff; 
            end
%        else;  Sm = abs( (pk2)/(F_peak)); if Sm>0.3; Sm=1; end
        end
    else Sm = 1; Rc{jj}.df= fp1-ff; 
    end
    Rc{jj}.Sm=Sm;

%    Rc{jj}.df= fp1-ff;

    Rc{jj}.df2= abs(Rc{jj}.df)/2; if Rc{jj}.df2 >1; Rc{jj}.df2=1; end; 
    
 %%%%% plot Arnould tonge with black and white depending on the above condition 
 figure(3);
 if peak1_cond & peak2_cond
     
 %hallo nicola!
 % hallo roberta and cornelius 
     plot(1000./p,flow,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[0,0,0]);hold on;
     Rc{jj}.syn=1;
 else
     plot(1000./p,flow,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[1,1,1]);hold on;
     Rc{jj}.syn=0;
 end
 
 %%%%% plot Arnould tonge with greyscale depending on the value Sm
 figure(4);
 %    plot(1000./p,flow,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',Rc{jj}.df2*[1,1,1]);hold on;
     plot(1000./p,flow,'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',Sm*[1,1,1]);hold on;
 
 end
 
 end
 
 cd('plots');
 
 figure(3); title(strcat('Cell ',num2str(nc) ,' mean freq',num2str(F_rest),'  B&W'));
 xlabel('External Frequency [Hz]')
 ylabel('Velocity [mm/s]');
 
 saveas(gca,strcat('Cell_',num2str(nc),'_BaW.jpg'))

 figure(4); title(strcat('Cell ',num2str(nc),' mean freq',num2str(F_rest), 'GreyScale'));
 xlabel('External Frequency [Hz]')
 ylabel('Velocity [mm/s]');
 saveas(gca,strcat('Cell_',num2str(nc),'_GreyScale.jpg'))
 
 cd(data_dir);
 % ylim([1.5,9]);
% figure(2);imagesc(Volts,1000./Periods,Mat1);
% figure();imagesc(Mat1);
Cell(nc).Rc=Rc; 

prompt= 'do you want to include in stat (1 for yes, 0 for no)? ';
good = (input(prompt));
Cell(nc).good= good;
%close all

%save('Cell.mat','Cell');

  %% plot different spectra

Hz_toplot=[10,12]; 
Volt_toplot=[4];
Period_toplot= unique(floor(1000./Hz_toplot));

ind=[];
cc=1;
for jj=1:(numel(Rc)-1) 
if any((Rc{jj}.Period)< max(Period_toplot)) & (Rc{jj}.Period> min(Period_toplot)) & any((Rc{jj}.Volt)== Volt_toplot); ind(cc)=jj;cc=cc+1;end
end

for jj=1:(numel(Rc)-1) 
if (Rc{jj}.Period) == 0 & (Rc{jj}.Volt)== Volt_toplot ; ind(cc)=jj;cc=cc+1;end
end


figure(1);
cleg=1;
 for jj=ind;
 figure(1);plot(Rc{jj}.fq(5:end),(Rc{jj}.m_pxx(5:end)),'LineWidth',2); hold on;
 %figure(2); plot(str2num(Hz{jj}),Rc{jj}.fp-str2num(Hz{jj}),'o'); hold on;
 if Rc{jj}.Period==0;
     leg{cleg}='at rest'
 else leg{cleg}=num2str(1000/Rc{jj}.Period);cleg=cleg+1;
 end
 end
 figure(1); legend(leg)
 xlabel('frequency [Hz]'),ylabel('fft over px intensity');
 xlim([Hz_toplot(1)-3,Hz_toplot(end)+3]);
 