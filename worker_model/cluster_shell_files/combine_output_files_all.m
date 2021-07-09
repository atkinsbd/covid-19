
% Declare function
function combine_output_files_all

    % find all the name types in the folder (will pick up anything with a
    % file that ends #1.mat)
    filetypes = dir('../../Results/worker_model/*#1.mat');

    % iterate over those name types
    for filetypes_it = 1:length(filetypes)
        load_and_save_filetype(filetypes(filetypes_it))
    end
end

function load_and_save_filetype(filetype)

    % find the actual filenames (i.e. filename_#1.mat, filename_#2.mat
    % etc)
    filenames = dir(['../../Results/worker_model/',filetype.name(1:end-5),'*']);

    % find what variables are in those files (assumes all files have
    % the same format)
    vars = whos('-file',['../../Results/worker_model/',filenames(1).name])

    % create new storage variables for the concatenated variables
    % assumes the concatenation happens along the second axis
    % i.e. newinf is of size (timesteps,runs_per_job,num_scenarios)
    % so newinf_combined will be of size (timesteps, total_runs, num_scenarios)
    runs = vars(1).size(2);
    total_runs = runs*length(filenames);
    for vars_it = 1:length(vars)
        size_var = vars(vars_it).size;
        if (length(size_var)>1)  % don't try to combine scalars,
            if (size_var(1) == runs) % Vectors of form runs by n_configs
                size_var(1) = total_runs;
            else
                size_var(2) = total_runs;
            end
            eval([vars(vars_it).name,'_combined = zeros(size_var);'])
        end
    end

    % initialise index
    start_idx = 1;

    % iterate over the actual filenames
    for filenames_it = 1:length(filenames)
        % Get final index for current replicate
        end_idx = filenames_it*runs;

        % load the variables in the file
        load(['../../Results/worker_model/',filenames(filenames_it).name])

        % iterate over the variables
        for vars_it = 1:length(vars)
            size_var = vars(vars_it).size;

            % enter the values into the storage variable
            if (size_var(1) == runs) % Vectors of form runs by n_configs
                eval_string = [vars(vars_it).name,'_combined(start_idx:end_idx,:'];
            else
                eval_string = [vars(vars_it).name,'_combined(:,start_idx:end_idx'];
            end
            if (length(size_var)>1)&&(size_var(1)~=1)&&(size_var(2)~=1) % don't try to combine scalars
                for i=3:length(size_var)
                    eval_string = [eval_string,',:'];
                end
                eval_string = [eval_string,') = ',vars(vars_it).name,';'];
                eval(eval_string)
            end
        end

        % Increment initial access index
        start_idx = end_idx + 1;
    end

    % Save combined storage arrays to output file
    save_filename = ['../../Results/worker_model/',filetype.name(1:end-7),'_combined.mat'];
    %save(save_filename,'*_combined');
    save(save_filename,'*_combined','-v7.3');
end
