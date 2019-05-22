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
List = squeeze(List);
for F = 1:numel(List)
    if(~isempty(List{F}))
        for i=1:size(List{F})
        tempName = ['file_',List{F}(i).name(1:end-4)];
        tempName = replace(tempName," ","_");
        Database_Files.(tempName) = List{F}(i);
        disp(fullfile(Database_Files.(tempName).folder,Database_Files.(tempName).name))
        
        % read abf file, two cannels.
        
        % in this section all the data relative to each file is stored into the
        % database
        
        answer = 1;
        
        the_file = [Database_Files.(tempName).folder '\' Database_Files.(tempName).name];
        clc;
        
        [Database_Files.(tempName).traces, Database_Files.(tempName).si, Database_Files.(tempName).h] = abfload(fullfile(the_file), 'start', 0, 'stop', 'e');
        Database_Files.(tempName).Time = (0:(Database_Files.(tempName).si/1000):((Database_Files.(tempName).h.sweepLengthInPts-1)*(Database_Files.(tempName).si/1000)))';
        
        Database_Files.(tempName).Voltages = squeeze ((Database_Files.(tempName).traces(:,2,1:end)));
        Database_Files.(tempName).Currents = squeeze ((Database_Files.(tempName).traces(:,1,1:end)));
        
        % with the next two lines I'm cancelling the field traces from the
        % database because is already stored in other two fields
        
        field = 'traces';
        Database_Files.(tempName) = rmfield(Database_Files.(tempName),field);
        
        % in these lines I'm deleting all the variables and saving all the 
        % analysed data by file
        
        File_Name = char(sprintf('%s.mat',Database_Files.(tempName).name (1:end-4)));
        field = {'bytes','isdir','datenum','h','si'};
        Database_Files.(tempName) = rmfield(Database_Files.(tempName),field);        
        j=1;
        
        %------------- from here, to the end of this cycle ---------------%
        % All the informations are stored in a cell array and in each 
        % columns there will be a different type of data( 1=Time,
        % 2=Currents, 3=Voltages, 4=sweep number, 5=file name, 6=folder tag
        
        for y=1:size(Database_Files.(tempName).Currents,2)
            Database_Files.(tempName).CompatibleSheet (j:j+(size(Database_Files.(tempName).Currents,1)-1),1) =  Database_Files.(tempName).Time(:,1);
            Database_Files.(tempName).CompatibleSheet (j:j+(size(Database_Files.(tempName).Currents,1)-1),4) = y;
            j = j+size(Database_Files.(tempName).Currents,1);
        end
            
        Database_Files.(tempName).CompatibleSheet (:,2) = reshape(Database_Files.(tempName).Currents,[],1);
        Database_Files.(tempName).CompatibleSheet (:,3) = reshape(Database_Files.(tempName).Voltages,[],1);
        
        Database_Files.(tempName).CompatibleSheet = num2cell(Database_Files.(tempName).CompatibleSheet);
        
        for y=1:numel(Database_Files.(tempName).Currents)
            Database_Files.(tempName).CompatibleSheet {y,5} =  Database_Files.(tempName).name;
            Database_Files.(tempName).CompatibleSheet {y,6} = Database_Files.(tempName).folder;
        end
        
        %-----------------------------------------------------------------%
        
        File = Database_Files.(tempName).CompatibleSheet;
        save (File_Name, 'File');
        close all;
        clearvars -except Database_Files F i List File
        i=i+1;
        end  
    end
end

% saving the Database file in the directory for consequent analysis
% save Database_File Database_Files

