function DCM = spm_dcm_load(P)
    % Loads a cell array of DCM filenames into an array for DCM
    
    DCM = cell(size(P));
    
    for s = 1:size(P,1)
        for m = 1:size(P,2)
            try
                model = load(P{s,m});
                DCM{s,m} = model.DCM;
            catch
                error('Failed to load model for subject %s model %m', s, m);
            end
        end
    end
end