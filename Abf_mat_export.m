% This script loads abf files and writes them to .mat files
% general initializing procedure
clear all

% acquiring all data from all files in the data folder

FullPathName = uigetdir; %(is the full name of the path)
cd(FullPathName);
FolderList = regexp([(genpath(FullPathName)) ';'],'(.*?);','tokens'); %(here are stored all the directories in which the files are contained)

% (in this cycle i'm creating a list of all the files contained in the folder;
% the variable "List" returns a list of all the files that have to be analysed;

for i = 1:size(FolderList,2)-1
    List{i} = dir([FolderList{1,i}{:} '\*.abf']);
end

% (After these lines you have a variable called list in which are stored
% information about the names and the paths of all the .abf files contained in the selected folder)

% I'm clearing now all the variables that I don't need from now on

clear ix;
clear i;
clear FolderList;
clear FullPathName;

% In this first part all the empty records from the variable List will be
% discarded

i=1;
for F = 1:numel(List)
    if(~isempty(List{F}))
        tempName = ['file_',List{F}.name(1:end-4)];
        tempName = replace(tempName," ","_");
        Database_Files.(tempName) = List{F};
        disp(fullfile(Database_Files.(tempName).folder,Database_Files.(tempName).name))
        
        % read abf file, two cannels.
        
        % in this section all the data relative to each file is stored into the
        % database
        
        answer = 1;
        
        the_file = [Database_Files.(tempName).folder '\' Database_Files.(tempName).name];
        clc;
        
        [Database_Files.(tempName).traces, Database_Files.(tempName).si, Database_Files.(tempName).h] = abfload(fullfile(the_file), 'start', 0, 'stop', 'e');
        Database_Files.(tempName).Time = (0:(Database_Files.(tempName).si/1000):((Database_Files.(tempName).h.sweepLengthInPts-1)/10))';
        
        Database_Files.(tempName).Voltages = squeeze ((Database_Files.(tempName).traces(:,2,1:end)));
        Database_Files.(tempName).Currents = squeeze ((Database_Files.(tempName).traces(:,1,1:end)));
        
        % with the next two lines I'm cancelling the field traces from the
        % database because is already stored in other two fields
        
        field = 'traces';
        Database_Files.(tempName) = rmfield(Database_Files.(tempName),field);
        
        % in these lines I'm deleting all the variables and saving all the analysed data by file
        
        File_Name = char(sprintf('%s.mat',Database_Files.(tempName).name (1:end-4)));
        field = {'bytes','isdir','datenum','h','si','name','folder','date'};
        Database_Files.(tempName) = rmfield(Database_Files.(tempName),field);
        File = Database_Files.(tempName);
        save (File_Name, 'File');
        close all;
        clearvars -except Database_Files F i List
        i=i+1;
        
    end
end

% saving the Database file in the directory for consequent analysis
save Database_File Database_Files