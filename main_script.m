%% Imaging for Neuroscience
% Homework - Group 2

% The goal of this homework is to investigate the relationship between 
% dynamic positron emission tomography and resting-state functional MRI 
% data. The following scans were simultaneously acquired with a hybrid
% PET/MRI system on a healthy subject. The tracer is irreversible.

% Aims of the study:
% 1. using a graphical approach, quantify the PET dynamic at the ROI level
%    and at the voxel level
% 2. compute the subject's functional connectivity matrix
% 3. investigate if there is a relationship between the PET estimates and
%    the computed connectivity measures at the ROI level.

clear 
close all
clc

%% PET analysis

%% loading PET data
load('DATA/PETdata.mat')

% The loaded files contains the following variables:
% Cbnew = plasmatic tracer concentration, corrected for metabolites [kBq/ml] {5439x1}
% mask = Hammers anatomical atlas resampled in the PET images space {128x128x95}
% PET = 4D dynamic PET data {128x128x95x24}
% time = PET hemi-scan times [min] {24x1}
% tnew = arterial sampling acquisition timing [min] {5439x1}

% We excluded Cpnew, delta, gamma since we didn't use them in our code

% adding the path to the NIfTI toolbox
addpath(fullfile('C:\NIfTI_20140122'))

% loading the grey matter tissue probability map {1x1 struct}
GM=load_untouch_nii('DATA/GM.nii.gz'); 


%% 1. Data preprocessing
% 1a. create a grey matter binary mask from the provided GM segmentation. 
%     Use a threshold equal to 0.3.

grey_matter=GM.img; % extracting the image from GM struct {128x128x95}

threshold=0.3;
grey_mask=grey_matter > threshold; % thresholding

% visualization of the grey matter and the respective mask
for i=1:size(grey_matter,3) 
    figure(1)
    sgtitle('Grey matter')
    subplot(121)
    imagesc(grey_matter(:,:,i)) 
    subplot(122)
    imagesc(grey_mask(:,:,i)) 
    pause(0.01)
end

% 1b. using the anatomical atlas and the grey matter mask, extract the average
%     tissue time activity curve (TAC) of the grey matter of each atlas region.

dim = size(PET); 
num_ROI = max(mask(:)); % number of ROIs in the atlas
tac = zeros(dim(4),num_ROI); % TAC initialization

for k=1:num_ROI
   idx = find(mask(:)==k); % indeces of the voxels beloging to the kth ROI
   for t=1:dim(4) 
       d = squeeze(PET(:,:,:,t)); % PET at instant t
       f = d.*grey_mask; % masking the slices of PET using grey_mask
       tac(t,k) = mean(f(idx)); % TAC value at instant t for the kth ROI
   end
end

% visualization of the TAC curves
figure('name','Time Activity Curve')
plot(time,tac)
title('TAC of each Atlas region')
xlabel('time [min]')
ylabel('kBeq/ml')
xlim([time(1) time(end)])

%% 2. Graphical methods
% 2a. based on the knowledge of the tracer kinetics, select the graphical 
%     method that better describes the PET data at the voxel level. 
%     If the Patlak method is used, estimate the net uptake (influx) rate 
%     constant (Ki). On the contrary, if the Logan method is used, estimate 
%     the tracer binding potential (BP) (right and left cerebellum as a 
%     single reference region). 
%     In both cases, use the unweighted linear least square estimator.

% Since the Logan Graphical Method can be used only for reversibile
% models and the tracer is said to be irreversible, we decided to
% describe the PET data at the voxel level using the Patlak method.
% Thus, in the following lines of code, we will estimate Ki.

% 2b. select a single t* value for all the ROIs: provide a justification 
%     for the choice of t* and report an example of Logan/Patlak plot 
%     (i.e. plot(X,Y) for a representative ROI).

% 2c. quantify PET data using the suitable graphical model (Platlak or Logan) 
%     both at the ROI and the voxel level.

%% Ki at the ROI level 

% interpolating the plasmatic tracer concentration, so that it is in
% the same time scale of the PET data
Cb_interp=interp1(tnew,Cbnew,time);

% computing X
integral_Cp=cumtrapz(time,Cb_interp);
X=integral_Cp./Cb_interp; % (the first three elements of X are NaN values)

% computing Y
Y=zeros(length(time),num_ROI);
for k=1:num_ROI

    C_measured=tac(:,k);

    Y(:,k)=C_measured./Cb_interp;
    
    % visualization
    figure(3)
    plot(X,Y(:,k),'o-')
    xlim([X(4) X(end)])
    xlabel('X')
    ylabel('Y')
    title('X,Y - at the ROI level')
    hold on
end

% visualization of a representative ROI
% we choose the 13th ROI, since the plot that shows Y as a function of X
% seems to respect the linear relationship that there should be between 
% the two after a certain time step.
figure('name','Representative ROI')
plot(X,Y(:,13),'o-') 
xlim([X(4) X(end)])
xlabel('X')
ylabel('Y')
title('X,Y - 13th ROI')

% by inspecting the plot, we can identify the 20th time point as the one
% after which the relationship between X and Y is linear.
idx_tstar=20;
tstar=time(idx_tstar); % t*

% keeping X and Y from idx_tstar to the end
X=X(idx_tstar:end);
Y=Y(idx_tstar:end,:);

% G matrix, used for the unweighted linear least squares
G=zeros(length(X),2);
G(:,1)=X;
G(:,2)=ones(length(X),1);

Ki_ROI=zeros(size(Y,2),1); % Ki initialization

for k=1:num_ROI

    pred_coeff=(G'*G)\(G'*Y(:,k)); % estimated coefficients

    Ki_ROI(k)=pred_coeff(1); % Ki in the kth ROI
end

% 2d. identify physiological estimates. How many ROIs show unreliable estimates?
idx_neg_ROI=find(Ki_ROI<0);
idx_nan_ROI=find(isnan(Ki_ROI));

% setting unreliable estimates to zero
Ki_ROI(idx_neg_ROI)=0;
Ki_ROI(idx_nan_ROI)=0;

disp(['Number of ROIs that show unreliable estimates: ', num2str(length(idx_neg_ROI)+length(idx_nan_ROI))])

%% Ki at the voxel level
% 2e. regarding the results at the voxel level: the physiological estimates 
%     must be organized in a 3D matrix (same size of mask variable). 
%     Replace the unreliable estimates with zeros.

% X and idx_tstar are the same as before

Ki_voxel=zeros(dim(1),dim(2),dim(3)); % Ki initialization

PET_masked=PET.*grey_mask; % masking PET data

for k=1:dim(3)
    for i=1:dim(1)
        for j=1:dim(2)
            
            % computing Y at the voxel level
            C_measured_voxel=squeeze(PET_masked(i,j,k,:));
            Y_voxel=C_measured_voxel./Cb_interp;

            % Y, from idx_tstar to the end
            Y_voxel=Y_voxel(idx_tstar:end);
            
            % G matrix, for the unweighted linear least squares
            G_voxel=zeros(length(X),2);
            G_voxel(:,1)=X;
            G_voxel(:,2)=ones(length(X),1);
            
            % estimated coefficients
            pred_coeff_voxel=(G_voxel'*G_voxel)\(G_voxel'*Y_voxel);
            
            % Ki evaluated in voxel (i,j,k)
            Ki_voxel(i,j,k)=pred_coeff_voxel(1);


        end
    end
end


% 2d. identify physiological estimates. How many voxels are discarded due 
%     to this selection?
idx_neg=find(Ki_voxel<0);
idx_nan=find(isnan(Ki_voxel));

% setting ureliable estimates to zero
Ki_voxel(idx_neg)=0;
Ki_voxel(idx_nan)=0;

disp(['Number of voxels that show unreliable estimates: ', num2str(length(idx_neg)+length(idx_nan))])


%% 2f. as for the results at the ROI level, are the physiological estimates 
%     of Ki are close to those obtainable by averaging the Ki values of the voxels 
%     contained within each ROI?
%     Provide a quantitative measure to evaluate the differences (at your choice, please motivate your choice)

Ki_average=zeros(num_ROI,1); % Ki initialization

for k=1:num_ROI

    idx=find(mask(:)==k); % voxels contained in the kth ROI

    % Ki as the mean of the voxels belonging to the kth ROI
    Ki_average(k)=mean(Ki_voxel(idx));
end

% scatterplot showing the correlation between Ki_ROI and Ki_average
figure('name','Scatterplot: Ki_ROI vs Ki_average')
scatter(Ki_ROI,Ki_average)
xlabel('Ki\_ROI')
ylabel('Ki\_average')
title('Scatterplot')

% quantitative measures
rmse_Ki=rmse(Ki_ROI,Ki_average); % Root Mean Square Error
r_Ki=corr(Ki_ROI,Ki_average);    % correlation
R2_Ki=r_Ki^2;                    % coefficient of determination

disp('Quantitative measures: ')
disp(['RMSE: ', num2str(rmse_Ki)])
disp(['Correlation: ', num2str(r_Ki)])
disp(['Coefficient of determination: ', num2str(R2_Ki)])

% Both the scatterplot and the evaluated coefficients show that averaging 
% the Ki estimates extracted from each voxel of a specific ROI provides a 
% result close to that obtained by evaluating Ki directly at the ROI level.


%% 2g. within the same figure, visualize the 45th slice of Ki matrix and
%     of the Hammers atlas after masking only the grey matter. 
%     Discuss the results: is the spatial pattern of the estimates homogeneous 
%     among the regions? 

Ki_45=Ki_voxel(:,:,45); % 45th slice of Ki matrix
slice_45=mask(:,:,45);  % 45th slice of the atlas
masked_slice_45=squeeze(slice_45.*grey_mask(:,:,45)); % keeping the gray matter

Ki_ROI_45=zeros(size(slice_45));
for k=1:num_ROI
   Ki_k=Ki_ROI(k);
   idx=find(masked_slice_45==k);
   Ki_ROI_45(idx)=Ki_k;
end

% visualization (at the voxel level)
figure('name','Comparison between Ki_45 and masked_slice_45')
subplot(121)
imagesc(Ki_45)
colorbar
colormap('hot')
title('Ki\_45')
subplot(122)
imagesc(masked_slice_45)
title('masked\_slice\_45')

% visualization (at the ROI level)
figure('name','Comparison between Ki_ROI_45 and masked_slice_45')
subplot(121)
imagesc(Ki_ROI_45)
colorbar
colormap('hot')
title('Ki\_ROI\_45')
subplot(122)
imagesc(masked_slice_45)
title('masked\_slice\_45')


% looking at the Ki at roi level it is possible to say that the activation is 
% not homogeneous among regions, infact it is possible to identify regions
% with a higher level of activity.

%% saving the required variables
save('PET_results_HW1_group2.mat','Ki_voxel','rmse_Ki','r_Ki','R2_Ki')
save('PET_data.mat', 'Ki_ROI','grey_mask')
clear all

%% fMRI analysis

addpath('DATA\')

%% 1. Data Loading 

WM_struct   = load_untouch_nii('WM.nii.gz'); % White Matter Data
WM_pmap   = double(WM_struct.img); % image data selection --> WM probability map
CSF_struct  = load_untouch_nii('CSF.nii.gz'); % CSF Data
CSF_pmap  = double(CSF_struct.img); % image data selection --> CSF propability map
PET  = load('PETdata.mat'); % for Hammers ATLAS
fMRI = load_untouch_nii('fMRI.nii.gz'); % fMRI Data
TR = 2.6; % time resolution for the fMRI
load('RegMov.mat') % contains reg_mov --> six column vector for the noise regression
load('PET_data.mat') %contain the grey_mask and the Ki_ROI from PET analysis

clear('WM_struct','CSF_struct')
%% 2. Binary Mask Creation

%White Matter

sl = 46; % Slice Number --> just to visualize something
figure('name','White Matter')
subplot(221)
imshow(WM_pmap(:,:,sl))
title('WM Original')
subplot(222)
histogram(WM_pmap(:,:,sl))
subplot(223)
WM_th = 0.75; % Threshold for WM --> keep white matter if probability>75%
WM = WM_pmap > WM_th; % Logical Matrix Extraction
imshow(WM(:,:,sl));
title('WM th')
subplot(224)
histogram(WM(:,:,sl))

%CSF
sl = 40;
figure('name','Cerebral Spinal Fluid')
subplot(221)
imshow(CSF_pmap(:,:,sl))
title('CSF Original')
subplot(222)
histogram(CSF_pmap(:,:,sl))
subplot(223)
CSF_th = 0.4;
CSF = CSF_pmap > CSF_th;
imshow(CSF(:,:,sl));
title('CSF th')
subplot(224)
histogram(CSF(:,:,sl))

clear("sl","WM_th","CSF_th","CSF_pmap","WM_pmap")
CSF = double(CSF); % Logical to Double Reconversion
WM  = double(WM);

%% 3. Image Erosion 2D

% WM
n_slices = size(WM,3); % total number of slices
WM_eroded = zeros(size(WM,1), size(WM,2),n_slices); % result array init.
structure = strel('disk',1); % define structuring element
for sl=1:n_slices % for every slice
    slice_data = WM(:,:,sl); % get nxm voxel image
    WM_eroded(:,:,sl) = imerode(slice_data,structure); % apply imerode and save
end
clear('sl',"n_slices",'slice_data',"WM")

%CSF
n_slices = size(CSF,3);
CSF_eroded = zeros(size(CSF,1), size(CSF,2),n_slices);
for sl=1:n_slices
    slice_data = CSF(:,:,sl);
    CSF_eroded(:,:,sl) = imerode(slice_data,structure);
end
clear('sl', 'structure',"n_slices",'slice_data','CSF')

sl = 46; % Slice Number --> just to visualize something
figure('name','Eroded Mask')
subplot(121)
imshow(WM_eroded(:,:,sl))
title('WM')

sl = 40;
subplot(122)
imshow(CSF_eroded(:,:,sl))
title('CSF')

%% 4. PCA Analysis

%WM - Masking Application

fMRI_Raw = fMRI.img; % Select 4D fmri raw data

nVolumes = size(fMRI_Raw,4); % Number of volumes acquired

WM_fMRI_Masked = zeros(size(fMRI_Raw));

for vol=1:nVolumes % for every volume
    WM_fMRI_Masked(:,:,:,vol) = fMRI_Raw(:,:,:,vol).*WM_eroded; % apply 3D mask
end


% WM - PCA 

nCols = size(WM_fMRI_Masked,1)*size(WM_fMRI_Masked,2)*size(WM_fMRI_Masked,3);

WM_pca.data = zeros(nVolumes,nCols);

counter = 0;
for slice=1:size(WM_fMRI_Masked,3) % for every slice
    for row=1:size(WM_fMRI_Masked,1) % row
        for col=1:size(WM_fMRI_Masked,2) % column
            counter = counter + 1;
            WM_pca.data(:,counter) = WM_fMRI_Masked(row,col,slice,:);
        end
    end
end
clear('counter', 'row', 'col', 'slice')

[WM_pca.coeff, WM_pca.score, ~,~,WM_pca.explained] = pca(WM_pca.data);

WM_first_score = WM_pca.score(:,1);
WM_exp_variance = WM_pca.explained(1);

clear("WM_pca","WM_fMRI_Masked","WM_eroded")

%CSF - Masking Application

CSF_fMRI_Masked = zeros(size(fMRI_Raw));

for vol=1:nVolumes % for every volume
    CSF_fMRI_Masked(:,:,:,vol) = fMRI_Raw(:,:,:,vol).*CSF_eroded; % apply 3D mask
end


% CSF - PCA 

nCols = size(CSF_fMRI_Masked,1)*size(CSF_fMRI_Masked,2)*size(CSF_fMRI_Masked,3);

CSF_pca.data = zeros(nVolumes,nCols);

counter = 0;
for slice=1:size(CSF_fMRI_Masked,3) % for every slice
    for row=1:size(CSF_fMRI_Masked,1) % row
        for col=1:size(CSF_fMRI_Masked,2) % column
            counter = counter + 1;
            CSF_pca.data(:,counter) = CSF_fMRI_Masked(row,col,slice,:);
        end
    end
end
clear('counter', 'row', 'col', 'slice','vol')

[CSF_pca.coeff, CSF_pca.score, ~,~,CSF_pca.explained] = pca(CSF_pca.data);

CSF_first_score = CSF_pca.score(:,1);
CSF_exp_variance = CSF_pca.explained(1);

clear("CSF_pca","fMRI_Raw","CSF_fMRI_Masked","CSF_eroded")

%Visualization of the first principal component
time_vector = [0:nVolumes-1].*TR;

figure('name','First Principal Component')
subplot(121)
plot(time_vector,WM_first_score,'b-o')
xlabel('time')
ylabel('activation level')
title('WM variance explained: '+string(WM_exp_variance))
subplot(122)
plot(time_vector,CSF_first_score,'b-o')
xlabel('time')
ylabel('activation level')
title('CSF variance explained: '+string(CSF_exp_variance))

%% 5. (B) ROI Time Activity Extraction and ROI 2d matrix cration

ATLAS = PET.mask; % Gets Hammer Atlas from file

GM_Mask = grey_mask; 

ATLAS_Masked = ATLAS.*GM_Mask; % apply mask

n_ROI = 73; % Number of ROIs (from PDF)

ROI_means = zeros(nVolumes,n_ROI); 
% each column will be the mean activity of a roi over time

ROI_2d = cell(1,n_ROI); %every cell will contain the 2D matrix of a roi
for roi = 1:n_ROI % for every roi
    roiMask = ATLAS_Masked == roi; % get Roi mask
    temp_matrix =[];
    for vol = 1:nVolumes
        squeezed_fMRI = squeeze(fMRI.img(:,:,:,vol));
        ROI_means(vol,roi) = mean(squeezed_fMRI.*(roiMask),'all');
        %for each volume extract the mean activity of the roi
        
        %creatrion of the 2D roi matrix as is needed the same for cycle the next step
        temp_matrix(vol,:) = squeezed_fMRI(roiMask);
        %there is no other way to do this because each 2d roi matrix has
        %different dimensions

    end
    ROI_2d{roi} = temp_matrix;
    clear('temp_matrix')
end

%visualization
figure('name','ROI fMRI mean activity')
plot(time_vector,ROI_means)
xlabel('time')
ylabel('ROIs')

%% 6. DENOISING
%for removing the non-neural undesired fluctuation is performed a
%linear regression on the rois 2D matrix (each row contains all the voxel 
%of the roi at a certain time). the regressor are about the mean activity
%of WM and CSF (their first pincipal component) and the movement artifact (reg_mov) 

reg_matrix = [reg_mov, WM_first_score, CSF_first_score]; %eigth regressors
reg_matrix = zscore(reg_matrix);

noisefree_ROI_2d = cell(1,n_ROI); %every cell will contain the noisefree 2d matrix of a roi
for roi = 1:n_ROI
    beta = (reg_matrix' * reg_matrix) \ reg_matrix' * ROI_2d{roi}; %LLS beta estimates
    noisefree_ROI_2d{roi} = ROI_2d{roi} - reg_matrix*beta;
end

clear('beta')

%visulaization of the regression matrix
figure('name','regression matrix for denoising')
imagesc(reg_matrix)
title('Regression Matrix')
colormap gray
colorbar

%is necessary an high pass filtering in order to minimize the contribution
%of non neural low frequency activity like respiration, cardiac pulsation
%and base drift caused by the instrumentation. for this aim is used a high
%pass filter with the cutoff frequency set at 1/(1.5*stimulation period)
%stimulation period is not provided so is used the SPM12 default 1/128 Hz

f_co = 1/128; %cutoff frequency
fs = 1/TR; %sample frequency
order = 3; %slow filter
[b,a] = butter(order,f_co/(fs/2),'high'); %butterworth filter 
figure('name','high pass Butterwoth filer')
freqz(b,a,512,fs)

filterd_ROI_2d = cell(1,n_ROI);
for roi = 1:n_ROI
    filterd_ROI_2d{roi} = filtfilt(b,a,noisefree_ROI_2d{roi});
    %filt all the column without change the signal phase
    %(each clumn in noisefreee_ROI_2d{roi} is the time series of a pixel
    %of that roi)
end

clear('a','b','order','roi','f_co','fs')


%% 7. check of preprocessing step

% hippocampus has ROI  = 1
% in the first cell of ROI_2d is present the time_temporal description
% (along columns) of each voxel that belongs to the first ROI

hippo_ROI = ROI_2d{1}; % every column is a voxel of the first ROI

% mediating on second dimension all elements in the same istant of time

mean_hippo_ROI = nanmean(hippo_ROI, 2);

figure('name','Hippocampus mean time series')
subplot(121)
plot(time_vector, mean_hippo_ROI)
title('pre-denoising')
xlabel('time')


hippo_ROI_filtered = filterd_ROI_2d{1};

% mediating on second dimension all elements in the same istant of time
% for the filtered roi

mean_hippo_ROI_filtered = nanmean(hippo_ROI_filtered, 2);

subplot(122)
plot(time_vector, mean_hippo_ROI_filtered)
title('after-denoising')
xlabel('time')


% plotting 

% mean_hippo_ROI subtracted with 'mean(mean_hippo_ROI)' to have 
% mean_hippo_ROI in the same level of mean_hippo_ROI_filtered

% mean_hippo_ROI_filtered has already mean=0 so I don't 
% have to subtract it for 'mean(mean_hippo_ROI_filtered)'

figure('name','Comparison')
plot(time_vector, mean_hippo_ROI - mean(mean_hippo_ROI),...
    time_vector , mean_hippo_ROI_filtered,'LineWidth',2)
title('Hippo time course after performing the denoising steps')
legend('Original signal','Preprocessed signal')
xlabel('time [min]')

clear('mean_hippo_ROI_filtered','mean_hippo_ROI','hippo_ROI')


%% 8. Computing functional connectivity matrix
% each pair of time series of the ROI

% extracting each mediated ROIs in a matrix
mean_ROIs = zeros(nVolumes, n_ROI);
for roi = 1:n_ROI % for every roi
    mean_ROIs(:,roi) = nanmean(filterd_ROI_2d{roi}, 2);
end

% correlation matrix between every mediated ROI
[FC, p] =   corr(mean_ROIs)  ; %compute FC matrix
figure('name','Correlation')
imagesc(FC)
title('correlation matrix')
colormap jet
colorbar

% appling Fisher z-transform to the coefficients
zFC = atanh(FC); 

figure('name','Fisher-transform')
subplot(121)
imagesc(zFC)
axis square
colormap jet
colorbar
title('z-Fisher transformed FC ')

subplot(122)
imagesc(p)
colormap jet
axis square
colorbar
title('p-values matrix')


%% 9. Clustering Coefficients

% 1st step: Binarizing the FC matrix
% We have to consider only the statistically significant functional
% connections, in order to find these we set the alpha-value at 0.05 and we
% also correct the alpha value by dividing it by the number of ROIs

%binarized_FC = p<=0.05/(n_ROI*(n_ROI-1)/2); % original Bonferroni correction
binarized_FC = p<=0.05/n_ROI;

figure('name','Binarized Correlation matrix')
imagesc(binarized_FC)
axis square
colorbar
title('Binarized FC matrix')

% 2nd step: Computing the clustering coefficient for each ROI
% The clustering_coef_bu function returns a vector containing the
% clustering coefficient for each node in the graph. Each element of this
% vector represents the clustering coefficient of the corresponding node

clustering_coeffs = clustering_coef_bu(binarized_FC);

% Creating a 3D map with Atlas masked dimensions, then assigning to every
% ROI in the matrix its clustering coefficient
cc_map = ATLAS_Masked;
for roi=1:n_ROI
    roiMask = ATLAS_Masked == roi;
    cc_map(roiMask) = clustering_coeffs(roi);
end

% Clustering coefficients visualization across different slices acquired
% over time with jet colormap
handle_player = implay(cc_map);
handle_player.Visual.ColorMap.MapExpression = 'jet';

%% saving the required variables
save('fMRI_results_HW1_group2.mat','zFC','clustering_coeffs')

%% CONCLUSION

[cf_ki_corr, cf_ki_p] = corr(clustering_coeffs,Ki_ROI);

figure('name','scatter KI-CF')
scatter(clustering_coeffs,Ki_ROI)
title('Correlation: ' + string(cf_ki_corr) + ' p value: ' + string(cf_ki_p))
xlabel('clustering coefficent')
ylabel("Roi's Ki")

disp("Roi's Ki and clustering coefficent correlation: " + num2str(cf_ki_corr))
disp("p value: " + num2str(cf_ki_p) + "--> p>0.05 --> not statisticaly correlated")


%looking at the graphs we can say that there's not correlation between the
%clustering coefficent and the Ki. our though is that they represent
%somenthing different: Ki is related to the activation/consumption of a
%certain area, while the clustering coefficent is related to how much the activation
%of thet area is linked to the activation of the near areas. 
%Therefore some areas can have an high Ki and a low clustering coeffincent 
%for example in case of tasks that demand the activation of segregated areas instead of
%diffuse activation of the brain the activated area will show an high Ki but
%a low clustering coefficent.
