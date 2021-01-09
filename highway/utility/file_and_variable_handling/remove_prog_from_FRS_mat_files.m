%% description
% Given a directory, this goes through each .mat file in the directory and
% loads it, removes the "prog" field from the "out" object, and re-saves it

%% user parameters
file_directory = '/Users/hannahlarson/Desktop/RSS19_data/dyn_10_stat_10'; %'/Users/shreyas/MATLAB/obstacle_avoidance/reach_sets/segway/combined_190116/dyn_8_stat_10/'

%% automted from here ;
% add a slash at the end of the directory just in case
file_directory = [file_directory,'/'] ;

% get files from the dir
files = dir(file_directory) ;

% go through each file and remove prob and re-save
for idx = 1:length(files)
    disp('---')
    
    clearvars -except files file_directory idx
    
    current_filename = [file_directory,files(idx).name] ;
    disp(['Current file: ',files(idx).name])
    
    % for each .mat file...
    if strcmp(current_filename(end-3:end),'.mat')
        disp(['Loading file ',num2str(idx)])
        tic
        load(current_filename)
        toc
        
        
        if exist('out','var') 
            if isfield(out,'prog')
                disp('Removing prog from out')
                out = rmfield(out,'prog') ;
            else
                disp('No prog field in out')
            end
        else
            disp('No out var')
        end
        
        if exist('out_hat','var') 
            if isfield(out_hat,'prog')
                disp('Removing prog from out_hat')
                out_hat = rmfield(out_hat,'prog') ;
            else
                disp('No prog field in out_hat')
            end
        else
            disp('No out_hat var')
        end
        
        disp('Re-saving file')
        save(current_filename)
    end
end

clear

disp('---')