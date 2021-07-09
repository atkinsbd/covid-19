%PURPOSE:
% Open MAT file and resave variables
%--------------------------------------------------------------------------

% Declare function
function resave_MAT_file(filedir,job_ID)

    % find all the name types in the folder (will pick up anything with a
    % file that ends #1.mat)
    filetypes = dir([filedir,'/*#',num2str(job_ID),'.mat']);
    
    % iterate over those name types
    for filetypes_it = 1:length(filetypes)
        clearvars -global -except filetypes*
        load([filedir,filetypes(filetypes_it).name])
        save([filedir,filetypes(filetypes_it).name],'-regexp', '^(?!(filetypes_it|filetypes|filedir|job_ID)$).') % Save all variables except filename
    end
end