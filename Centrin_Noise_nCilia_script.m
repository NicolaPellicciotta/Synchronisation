clear all;
%data_dir='/media/np451/Seagate Backup Plus Drive/DATA/1.11.18/'  %%% experiment in DMEM P/S
data_dir='/media/np451/Seagate Expansion Drive/29.10.18/';
%Positions1={'P1','P2','P3','P4','P5','P6','P7','P8','P9','P10'}
data_dir='/media/np451/Seagate Backup Plus Drive/DATA/1.11.18/' ;
%Positions2={'P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12','P13'}
data_dir='/media/np451/Seagate Backup Plus Drive/DATA/26.11.18/';
%Positions3={'P1','P2','P3','P4','P5','P6','P7'};

cc=1;
for dd=[1:30];
directories{cc}=strcat('P',num2str(dd)); cc=cc+1;
end
cd(data_dir);
%% for counting number of cilia in a cell
for dd=1:numel(directories)
    dd=14;
    disp(dd);
    cd(strcat(data_dir,directories{dd}));
    load('Cell.mat');

    
    d=dir('*P0*movie');
    mo=moviereader(d(1).name);
    fs=mo.read();
    fsbw= double(fs(:,:,2:300))- movmean(fs(:,:,2:300),100,3);
    for kk=1:size(fsbw,3)
    fsbw2(:,:,kk)= wiener2(fsbw(:,:,kk),[5,5]);
    end
    k=0;
    scroll_stack(fsbw2)
    prompt= 'number of cilia';
    Ncilia = (input(prompt));
    if numel(Cell)>1;
        Cell(1).Ncilia=Ncilia;    
    else
        Cell(1).Ncilia=Ncilia;
    end
%    save('Cell_cilia.mat','Cell');
    close all;
    
end
%% for addin Centrin and Noise
for dd=1:numel(directories)
    disp(dd);
    cd(strcat(data_dir,directories{dd}));
    load('Cell.mat');
%     load('BW.mat');
%     Noise=noise_fft(Cell(1).F_rest,BW{1});
%     Cell(1).Noise=Noise;
    if numel(Cell)>1;
        Cell(1).Centrin= Centrin_fun();
        saveas(gcf,'plots/centrin.pdf');
        %%%% make images    
    else
        Cell.Centrin= Centrin_fun();
        saveas(gcf,'plots/centrin.pdf');
    end
        save('Cell.mat','Cell');
        close all;
end

%%  SPATIAL NOISE
for dd=1:numel(directories)

    disp(dd);
    cd(strcat(data_dir,directories{dd}));
    load('Cell.mat');
    box_size=6;
    Sp_noise=spatial_noise_fft(Cell(1).F_rest,box_size);
    Cell(1).Sp_noise=Sp_noise;
    
    %%% make figure spatial noise
    
    subplot(1,3,1);
    X= Cell(1).Sp_noise.X(1:end-1,1:end-1);
    Y= Cell(1).Sp_noise.Y(1:end-1,1:end-1);
    imagesc(Y(1:end,1),X(1,1:end),Cell(1).Sp_noise.F')
    axis image
    colorbar();
    subplot(1,3,2);
    imagesc(Cell(1).Sp_noise.s_roi);
    axis image
    subplot(1,3,3);
    imagesc(Cell(1).Sp_noise.s_bin);
    axis image
    saveas(gcf,'plots/sp_noise_freq_distribution.pdf');
    close all;
    
    %%%%% make another figure
    
    s_roi=Cell(1).Sp_noise.s_roi;
    mask_tot= Cell(1).Sp_noise.mask_tot;

    mask_overlay(s_roi,mask_tot,[1,0,0],0.3)
    hold on;
    for jj=1:numel(Cell(1).Sp_noise.Box)
        box=Cell(1).Sp_noise.Box(jj);
        pos= Cell(1).Sp_noise.Box(jj).vec_pos;
        if box.good ==1;
            rectangle('Position',[pos(2),pos(1),pos(4)-pos(2),pos(3)-pos(1)]);
        end
    end
 saveas(gcf,'plots/sp_noise_mask.pdf');
 
    
    
    
     save('Cell.mat','Cell');
     clear Cell;
     close all;
end