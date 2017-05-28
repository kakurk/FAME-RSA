% The master script for the FAME RSA analsysis
%
% Written by Kyle Kurkela, kyleakurkela@gmail.com
% May-June, 2017

run('/storage/home/kak53/Downloads/set_kyle_mail_prefrences.m')

runSpec = true;
runEst  = true;
runST   = true;
runRSA  = true;

%% Step 1: Specify the Model

if runSpec

sendmail('kyleakurkela@gmail.com', 'Specifying...')

pause off

try
    SpecifyModel;
catch
    sendmail('kyleakurkela@gmail.com', 'Problem with Specifying')
end

pause on

end

%% Step 2: Estimate Model

if runEst

sendmail('kyleakurkela@gmail.com', 'Estimating...')

try
    EstimateModel;
catch
    sendmail('kyleakurkela@gmail.com', 'Problem with Estimating')
end

end

%% Step 3: Generate Single Trial Models

if runST

sendmail('kyleakurkela@gmail.com', 'Single Trial...')

try
    generate_spm_singletrial;
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
catch
    sendmail('kyleakurkela@gmail.com', 'Problem with RSA')
end

end
