%%%%%%%%%%%%% Get all files in the directory and subdirectories %%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dirName: Input folder 
% search_sac: term of research 
% fileList: list of the files/paths
% fileFolder: list of the folders

function [fileList,fileFolder] = getAllFile(dirName,search_sac)

  dirData = dir(strcat(dirName,'*/',search_sac));%Get the data for the current directory
  dirIndex = [dirData.isdir];% Find the index for directories
  fileList = {dirData(~dirIndex).name}';%Get a list of the files
  fileFolder = {dirData(~dirIndex).folder};%Get a list of the subfolders
  if ~isempty(fileList)
% loop through the list of file
      for k=1:length(fileList)
          fileList(k)=strcat(fileFolder(k),'/',fileList(k));%Append path to k-th file
      end
  end  
end