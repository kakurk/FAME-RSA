% The master script for the FAME RSA analysis
%
% Written by Kyle Kurkela, kyleakurkela@gmail.com
% July-August, 2017

run('/storage/home/kak53/Downloads/set_kyle_mail_prefrences.m')
addpath(genpath('/gpfs/group/nad12/default/nad12/spm12'))

runSpec = true;
runEst  = true;
runST   = true;
runRSA  = true;

models  = {'001_visualcategory' '002_trialtype' '003_memory'};

for m = 1:lenth(models)

    root = fullfile('/gpfs/group/nad12/default/nad12/FAME8/RSA/models', models{m});
    
    %% Step 1: Specify the Model
    
    if runSpec

        sendmail('kyleakurkela@gmail.com', 'Specifying...')

        pause off

        try
            if m == 1
                SpecifyModel_visualcategory;
            elseif m == 2
                SpecifyModel_trialtype;
            elseif m == 3
                SpecifyModel_memory;
            else
                error('No Model Specified')
            end
        catch
            sendmail('kyleakurkela@gmail.com', 'Problem with Specifying')
        end

        pause on

    end

    %% Step 2: Estimate Model

    if runEst

        sendmail('kyleakurkela@gmail.com', 'Estimating...')
        
        try
            EstimateModel(root);
        catch
            sendmail('kyleakurkela@gmail.com', 'Problem with Estimating')
        end

    end

    %% Step 3: Generate Single Trial Models

    if runST

        sendmail('kyleakurkela@gmail.com', 'Single Trial...')

        try
            generate_spm_singletrial(root);
        catch
            sendmail('kyleakurkela@gmail.com', 'Problem Generating Single Trial')
        end

    end

    %% Step 4: run RSA

    if runRSA

        % copy over the ROIs
        copyfile('/gpfs/group/nad12/default/nad12/FAME8/Analysis_ret/FAMEret8RSA_hrf/*.nii', '/gpfs/group/nad12/default/nad12/FAME8/RSA/models/SingleTrialModel/')

        sendmail('kyleakurkela@gmail.com', 'RSA...')

        try
            run_CANLAB_rsa;
        catch ME
            sendmail('kyleakurkela@gmail.com', 'Problem with RSA')
            rethrow(ME)
        end

    end

end
