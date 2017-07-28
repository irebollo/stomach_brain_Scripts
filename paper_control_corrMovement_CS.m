%{
 
paper_control_corrMovement_CS

Perform correlation between the framewise displacement betas and coupling
strength in all voxels of the gastric network, reported in section
Gastric-coupling in gastric network is not related to head micromovements


Input:
Betas FWD regression
    Y:\Subjects\Subject13\Timeseries\Other\BetasInstantaneousMovementEstimates_S_13

Coupling strenght x voxel
    Y:\Subjects\Subject13\Timeseries\GlobalSignal\couplingStrenght_csfr_S_13_kw3_fir2_fspread_015_fOrder_5_tw_15

Output
%}

%% Filenames and parametars

plotDir = strcat (global_path2root,'Figures',filesep,'Global',filesep);
plotFilename = strcat(plotDir,'EffectMovementonGSandCS');

insideBrain = tools_getIndexBrain('inside');
[gasnet] = global_getGastricNetwork;
gasnet = gasnet(insideBrain);
gasnet = gasnet(:);

%% Do


% initialize matrices

BetasMovement = zeros(length(subjects),length(insideBrain));
cs_allsubjects = zeros(length(subjects),length(insideBrain));

Mov_cs_r_wholeBrain = zeros(length(subjects),1);
Mov_cs_r_gasnet = zeros(length(subjects),1);


% iterate through subjects, 
for iSubj = 1:length(subjects)
subj_idx = subjects(iSubj)


% load betas of framewise displacement
filename_BetasInstantaneousMovementEstimates = global_filename(subj_idx,cfgMain,'filename_BetasInstantaneousMovementEstimates');
load(filename_BetasInstantaneousMovementEstimates)
BetasMovement(iSubj,:) = betas_instantaneousMovement;

% load coupling strenght
couplingStrenghtFilename_csfr = global_filename(subj_idx,cfgMain,'couplingStrenghtFilename_csfr');
load(couplingStrenghtFilename_csfr);
cs_allsubjects(iSubj,:) = couplingstrenght;


% identify and exclude NaNs
ind_nan_mov = isnan(BetasMovement(iSubj,:));
ind_nan_cs = isnan(cs_allsubjects(iSubj,:));
not_nan = ~ind_nan_mov & ~ind_nan_cs ;

[r p ] = corrcoef(BetasMovement(iSubj,not_nan),cs_allsubjects(iSubj,not_nan));
Mov_cs_r_wholeBrain(iSubj) = r(3)

% we are only interested in the correlation inside the gastric network
BetasMovement_gasnet = BetasMovement(gasnet);
couplingstrenght_gasnet = couplingstrenght(gasnet);

ind_nan_mov_gasnet = isnan(BetasMovement_gasnet);
ind_nan_cs_gasnet = isnan(couplingstrenght_gasnet);
not_nan = ~ind_nan_mov_gasnet & ~ind_nan_cs_gasnet ;

[r p ] = corrcoef(BetasMovement_gasnet(not_nan),couplingstrenght_gasnet(not_nan));
Mov_cs_r_gasnet(iSubj) = r(3)

end

%average
BetasMovement_mean = nanmean(BetasMovement);
cs_allsubjects_mean = nanmean(cs_allsubjects);
    

% group average
figure
    [r,p]=corrcoef(cs_allsubjects_mean(gasnet),BetasMovement_mean(gasnet));
    plot(cs_allsubjects_mean(gasnet),BetasMovement_mean(gasnet),'ok')
    title([ 'Group r movement and CS gasnet ' num2str(r(3)) ' p ' num2str(p(3))],'fontsize',20)
    xlabel('Coupling Strenght','fontsize',20)
    ylabel('Instantaneous Movement Beta','fontsize',20)
    h=gca
    set(h,'fontsize',16)
    h=lsline
    set(h,'LineWidth',4)
    set(h,'Color',[0.5 0.5 0.5])

% individual subjects
    figure
    Mov_cs_r_gasnet_fisherZ = fisherz(Mov_cs_r_gasnet)
    violin(Mov_cs_r_gasnet_fisherZ)
    h=gca
    set(h,'fontsize',16)
    [h p ci stats] = ttest(Mov_cs_r_gasnet_fisherZ)
    title(['Individual r values Gastric Network p ' num2str(p)])
    
% values reported
mean (Mov_cs_r_gasnet_fisherZ)
std (Mov_cs_r_gasnet_fisherZ)
        