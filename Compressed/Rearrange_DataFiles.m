function [Corrected_DataFile] = Rearrange_DataFiles(tone_number)
% This function rearranges the sensors using Pablo's Rearrange_Data
% function for one tone at a time

    %The MyFolder variable in line 8 will need to be changed on a
    %computer-by-computer basis.
    myFolder = 'C:\Users\MaxSlezak\Google Drive\Classes\2018-2019\Senior Design\Measurement_Data';
    myFolder = fullfile(myFolder,string(tone_number));
    % Check to make sure that folder actually exists.  Warn user if it doesn't.
    if ~isfolder(myFolder)
      errorMessage = sprintf('Error: This folder does not exist:\n%s', myFolder);
      uiwait(warndlg(errorMessage));
      return;
    end
    % Get a list of all files in the folder with the desired file name pattern.
    filePattern = fullfile(myFolder, '*.mat'); 
    theFiles = dir(filePattern);
    %
    correct_order = [16	32	15	31	14	30	13	29	12	28	11	27	10	26	9   25 ...
                     8	24	7	23	6	22	5	21	4	20	3	19	2	18	1	17 ...	
                     33	49	34	50	35	51	36	52	37	53	38	54	39	55	40	56 ...
                     41	57	42	58	43	59	44	60	45	61	46	62	47	63	48	64];

    Corrected_DataFile = cell(length(theFiles),1);
    for k = 1:length(theFiles)
      baseFileName = theFiles(k).name;
      fullFileName = fullfile(myFolder, baseFileName);
      totalData = load(fullFileName);
      Corrected_DataFile{k} = totalData(:,correct_order););
    end
end