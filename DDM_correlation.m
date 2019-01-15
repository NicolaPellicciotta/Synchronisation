function [ Iqtau ] = DDM_correlation( frame_stack)
%DDM_core Computational core for DDM calculation
%   Detailed explanation goes here

%% input check



if size(frame_stack,1)~=size(frame_stack,2)
    disp('Warning: non-square frames, only a square region will be actually analysed');
    frame_size = min([size(frame_stack,1),size(frame_stack,2)]);
    frame_stack = frame_stack(1:frame_size,1:frame_size,:);
end



frame_size = min([size(frame_stack,1),size(frame_stack,2)]);
N_px_frame = frame_size*frame_size;
fft2_norm_factor = 1/N_px_frame;
N_frames = size(frame_stack,3);

%% general-purpose variables
if nargin<3
max_tau = floor(N_frames/2);
end

% max_tau = N_frames-1;%floor(N_frames/2);    %also equal to the number of lags
max_q = floor(frame_size/2);

%% distance map for fast radial average

jj = repmat((1:frame_size),frame_size,1);
ii=jj';
cc = max_q+1;
distance_map =fftshift(round(sqrt((ii-cc).*(ii-cc)+(jj-cc).*(jj-cc)))+1); %if we fftshit here is only 1 time instead of Nframes/2
distance_map = distance_map(:);

dist_counts = accumarray(distance_map,ones(frame_size*frame_size,1));

%% actual DDM calculation
%(difference of frames, then |FFT|^2, average and radial average)

Iqtau = zeros(max_q, max_tau, 'double');
% err_Iqtau = zeros(max_q, max_tau, 'double');

% Iqmattau = zeros(frame_size, frame_size, max_tau, 'double');

    accum_abs_FT_diff_image = zeros(frame_size,'single');
    
    cc=0;

for tau=1:max_tau


        % on-the-fly average of fft2 of difference images
        i=1;
        tempdiff = fft2( single(frame_stack(:,:,i)) - single(frame_stack(:,:,i+tau)) ) * fft2_norm_factor; %fft2 of difference image
        accum_abs_FT_diff_image = accum_abs_FT_diff_image + real( tempdiff ).^2 + imag( tempdiff ).^2;
        cc=cc+1;



    
    averaged_abs_FT_diff_image = accum_abs_FT_diff_image./cc; %average on initial times (ie on couple of frames at same lagtime)
        
    oneD_power_spectrum = accumarray(distance_map,averaged_abs_FT_diff_image(:))./dist_counts;	%radial average
    
    Iqtau(:,tau) = oneD_power_spectrum(2:max_q+1);	%fill each column of the output with the 1D power spectrum just calculated
end
% toc;    
end






