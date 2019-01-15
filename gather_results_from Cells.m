%%% gather results from Cell array in plots

%data_dir='/media/np451/Seagate Expansion Drive/ependymal/syn/29.7.18/P1';
%data_dir='/media/np451/Seagate Expansion Drive1/DATA/13.10.18/P5/';
data_dir='/media/np451/Seagate Expansion Drive1/DATA/16.10.18/'
Positions=['P3','P4']

for ps=1:numel(Positions)
cd(Position(ps));
a=load('Cell.mat');

cd(data_dir);    
end
data_dir='/media/np451/Seagate Expansion Drive1/DATA/16.10.18/P3/';
cd(data_dir);
load('Cell.mat');
 
%% find  
Ncells= size(Cell,2);
Volts=[2:0.5:6];
%%% which frequency
cc=1;nc=1;
 for jj=1:(numel(Cell(nc).Rc)-1); if Cell(nc).Rc{jj}.Volt == 2; 
         Periods(cc)= Cell(nc).Rc{jj}.Period; cc=cc+1;
     end;end;
%%%%% occhio un po di casini con quando prende gli infiniti!!! da P=0;

 Freq=1000./double(Periods); Freq=sort(Freq);Freq=Freq(1:end-1);
         
%V=3.5;
Syn=zeros([numel(Volts),numel(Freq)]);
for nc=1:Ncells
    Cell(nc).Freq=Freq;
    for jj=1:numel(Cell(nc).Rc)-1;
        Cell(nc).Rc{jj}.Freq= 1000./double(Cell(nc).Rc{jj}.Period);
        p= Cell(nc).Rc{jj}.Period;
        if p~=0% & p<100;
        v=Cell(nc).Rc{jj}.Volt;
        f=Cell(nc).Rc{jj}.Freq;
        Syn(v==Volts,f==Freq)  =  Cell(nc).Rc{jj}.syn;
        end
     end
    Cell(nc).Syn=Syn;
    %%%%%% integral for find the synchonisation region epsilon
    for v=1:numel(Volts); Cell(nc).eps(v)=  trapz(Freq,Cell(nc).Syn(v,:)); end;
end

%% plot results with image
    whichVolt=5;
    v_eps= (Volts==whichVolt);
    Reds= cat(3,s(:,:,1),0*s(:,:,1),0*s(:,:,1));  
    imshow(imadjust(s));
    BWT=zeros(size(s));
    for nc=1:size(Cell,2);
        if Cell(nc).good
        BWT= BWT+BW{nc};
        [maskx,masky]=find(BW{nc});
        text(masky(1)-30,maskx(1)-30,num2str(Cell(nc).eps(v_eps),2),'Color','red','FontSize',24)
        end
        end
    hold on
    h = imshow(Reds); % Save the handle; we'll need it later
    hold off        
    set(h, 'AlphaData', BWT); 
    saveas(gca,'cell_analysed_epsilon.jpg')




%% epsilon vs frequency
figure();
whichVolt=5;
v_eps= (Volts==whichVolt);

cc=1;
for nc=1:Ncells;
    if Cell(nc).good
    plot(Cell(nc).F_rest,Cell(nc).eps(v_eps),'o');hold on;
    leg{cc}=num2str(nc); ;cc=cc+1;
    end
end
legend(leg)
xlabel('Frequency at rest[Hz]');
ylabel('Epsilon [Hz] synch region');
title(strcat('Epsilon vs Frequency V=',num2str(whichVolt)));
 saveas(gca,'epsilon_frequency.jpg');
%% epsilon vs Volts
M_epsilon=zeros(size(Volts));
for vv=1:numel(Volts);
    mean_epsilon_v=0; cc=1;
    for nc=1:Ncells;
    if Cell(nc).good
    mean_epsilon_v = mean_epsilon_v + Cell(nc).eps(vv); cc=cc+1;
    end
    end
    mean_epsilon_v/cc
    M_epsilon(vv)= mean_epsilon_v./cc;
end
figure()
plot(Volts,M_epsilon,'o');
xlabel('Voltage');
ylabel('epsilon[Hz]');
 saveas(gca,'epsilon_voltage.jpg');
 %% epsilon vs N_cilia
figure();
v_eps=3;
cc=1;
for nc=1:Ncells;
    if Cell(nc).good
    plot(Cell(nc).N_cilia,Cell(nc).eps(v_eps),'o','MarkerSize',8,'LineWidth',3,'MarkerFaceColor',[1,1,1]);hold on;
    leg{cc}=num2str(nc); ;cc=cc+1;
    end
end
legend(leg)
xlabel('number of cilia');
ylabel('Epsilon [Hz] synch region');

%% flagellar noise, tryng to measure flagellar noise from P0
for nc=1:Ncells
fq= Cell(nc).Rc{1}.fq;
F_rest= Cell(nc).F_rest;
fq_min= F_rest-5;
fq_max= F_rest+5;
f_range=fq> (fq_min) & fq<(fq_max);
m_pxx= Cell(nc).Rc{1}.m_pxx;
x=fq(f_range)
y=m_pxx(f_range)

f = fit(x.',y.','gauss1'); 
noise =f.c1/(f.b1);

figure(); plot(f,x,y);
xlabel('freq[Hz]'); ylabel('fft');
title(strcat('Cell ',num2str(nc),'  noise  ',num2str(noise)));
Cell(nc).noise=noise; 
Cell(nc).noise_fit=f;
end

%% %% epsilon vs flagellar noise
figure();

v_eps= (Volts==4);
cc=1;
for nc=1:Ncells;
    if Cell(nc).good
        Cell(nc).noise= var(1./Cell(nc).f_rest(1:4))./(1/Cell(nc).F_rest)^2;
    plot(Cell(nc).noise,Cell(nc).eps(v_eps),'o','MarkerSize',8,'LineWidth',3,'MarkerFaceColor',[1,1,1]);hold on;
    leg{cc}=num2str(nc); ;cc=cc+1;
    end
end
legend(leg)
xlabel('noise (std(f)/<f>2)');
ylabel('Epsilon [Hz] synch region');


%%

figure();
v_eps= (Volts==4);
cc=1;
for nc=1:Ncells;
    if Cell(nc).good
    plot(Cell(nc).N_cilia,Cell(nc).noise,'o','MarkerSize',8,'LineWidth',3,'MarkerFaceColor',[1,1,1]);hold on;
    leg{cc}=num2str(nc); ;cc=cc+1;
    end
end
legend(leg)
xlabel('noise ');
ylabel('Epsilon [Hz] synch region');