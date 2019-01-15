%%% gather results from Cell array in plots

 clear all;
% dir1='/media/np451/Seagate Backup Plus Drive/DATA/Cyto-D/Cyto_D.29.10.18/0hr/';
% dir2= '/media/np451/Seagate Backup Plus Drive/DATA/Cyto-D/Cyto_D.29.10.18/7hr/';
% dir3= '/media/np451/Seagate Backup Plus Drive/DATA/Cyto-D/Cyto_D.29.10.18/24hr/';
% dir_plots= '/media/np451/Seagate Backup Plus Drive/DATA/Cyto-D/Cyto_D.29.10.18/';
% Positions1={'P1','P2','P3','P4','P5','P6','P7','P8','P9','P10'};
% Positions2={'P2','P3','P4','P5','P6','P7','P9','P10'};
% Positions3={'P3','P4','P7'};

% PT1=strcat(dir1,Positions1);PT2=strcat(dir2,Positions2);PT3=strcat(dir3,Positions3);

dir{1}='/media/np451/Seagate Backup Plus Drive/DATA/15.11.18/0hr/';
dir{2}= '/media/np451/Seagate Backup Plus Drive/DATA/15.11.18/1hr/';
dir{3}= '/media/np451/Seagate Backup Plus Drive/DATA/15.11.18/6hr/';
%Positions1={'P1','P2','P3','P4','P5'}
%Positions2={'P1','P2','P3','P4','P5'};
%Positions3={'P4'};
Positions={'P1','P2','P3','P4','P5'}
%PT1=strcat(dir1,Positions1);PT2=strcat(dir2,Positions2);PT3=strcat(dir3,Positions3);

%cd(dir_plots)


for lag=1:3;
    pos_nc=0;
    cd(dir{lag}) 
%for ps=1:numel(eval(strcat('PT',num2str(lag)))) ;
for ps=1:numel(Positions) ;
    cd(dir{lag}) 
    pos_nc=pos_nc+1;  
%     pos=(eval(strcat('Positions',num2str(lag),'(',num2str(ps),')')));
%     pos=str2num(pos{1}(2:end));
%     temp_dir= eval(strcat('PT',num2str(lag),'{',num2str(ps),'}'));
    temp_dir=strcat('P',num2str(ps));
    if exist(temp_dir,'dir')==7
        disp(ps)
        cd(temp_dir);

        a=load('Cell.mat');
        temp_Ncells= (size(a.Cell,2));

        for nc=1:temp_Ncells

            pos_nc= pos_nc +nc-1;

            Cell(pos_nc,lag).ps=     ps;
            Cell(pos_nc,lag).nc=     nc;
            Cell(pos_nc,lag).Rc=     a.Cell(nc).Rc; 
            Cell(pos_nc,lag).good=   a.Cell(nc).good;
            Cell(pos_nc,lag).F_rest= a.Cell(nc).F_rest;
            Cell(pos_nc,lag).Pos=    temp_dir;
        end
        
        if lag==1; Cellxpos(ps)=nc;    end
    else   %%%%% if the directory do not exist or Cell do not exist
         for nc=1:Cellxpos 
            pos_nc=                 pos_nc +nc-1;
            Cell(pos_nc,lag).ps=    ps;
            Cell(pos_nc,lag).nc=    nc;
            Cell(pos_nc,lag).good=  0;
         end
    end
end

end
%% 
%%%%%% loading the variable Freq on Cell
Volts=[2:0.5:6];

%%% which frequency

for nc=1:numel(Cell)
    cc=1;Periods=[];clear Freq;
 for jj=1:(numel(Cell(nc).Rc)-1); if Cell(nc).Rc{jj}.Volt == 5; 
         Periods(cc)= Cell(nc).Rc{jj}.Period; cc=cc+1;
     end;end;
%%%%% occhio un po di casini con quando prende gli infiniti!!! da P=0;

 Freq=1000./double(Periods); Freq=sort(Freq);Freq=Freq(1:end-1);
 Cell(nc).Freq=Freq;
%%%%%
end



for nc=1:numel(Cell)
    Freq=Cell(nc).Freq;
    Syn=zeros([numel(Volts),numel(Freq)]);

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

%% treatment time vs epsilon
figure()
whichVolt=4;
v_eps= (Volts == whichVolt);
time=[0,1,6]
cc=1;

for nc=1:size(Cell,1)
    eps_cyto=[];
    goods=[];
    for lag=1:3;
        goods(lag)=Cell(nc,lag).good;
        if goods(lag)==1; eps_cyto=cat(1,eps_cyto,Cell(nc,lag).eps(v_eps));end
    end;
        leg{cc}=Positions{Cell(nc,1).ps};cc=cc+1;

        plot(time(logical(goods)),eps_cyto,'-o'); hold on;     
end
legend(leg)
xlabel('treatment time [hr]');
ylabel('Epsilon [Hz] synch region');


%%
time=[0,7]
cc=1;
whichVolt=4;
v_eps= (Volts == whichVolt);
for nc=1:size(Cell,1)
    if Cell(nc,1).good==1 & Cell(nc,2).good==1
        leg{cc}=Positions1{nc};cc=cc+1;
        freq_treatment=[Cell(nc,1).F_rest,Cell(nc,2).F_rest];
        eps_treatment=[Cell(nc,1).eps(v_eps),Cell(nc,2).eps(v_eps)];
        plot(diff(freq_treatment),diff(eps_treatment), '-o'); hold on;     
    end
end
legend(leg)
xlabel('freq_drop [Hz]');
ylabel('freq drop[Hz]');
%% treatment time vs epsilon
figure()
whichVolt=4;
v_eps= (Volts == whichVolt);
time=[0,7]
cc=1;
for nc=3:size(Cell,1)
    if Cell(nc,1).good==1 & Cell(nc,2).good==1
        leg{cc}=Positions1{nc};cc=cc+1;
        eps_treatment=[Cell(nc,1).eps(v_eps),Cell(nc,2).eps(v_eps)];
        plot(time,eps_treatment,'-o'); hold on;     
    end
end
legend(leg)
xlabel('treatment time [hr]');
ylabel('Epsilon [Hz] synch region');

