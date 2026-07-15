function [eeg, sFreq, ADBitVolts] = readCRTsd(fname, nBlockMin, nBlockMax, readParent)
% Reads a CSC file (new NT format) and returns a tsd
%
% function tsd = ReadCR_tsd(cr)
%
%
% INPUT:
%       fname ... full filename of Cheetah_NT CSC*.dat file     
%  
% OUTPUT:
%
%       tsd of the csc data.
%
% ADR 1999 called CR2EEG
% status PROMOTED
% version 1.0
% cowen Sat Jul  3 14:59:47 1999
% lipa  modified for NT   Jul 18 1999

% o Got rid of the diplay progress and rounded the timestamps
% o Fixed the dT to be in timestamps 
% o Added a dummy value to the end of cr.ts
% o Made the for look go to nBlocks and not nBlocks -1 

%       ReadCR_nt returns 2 arrays and 1 double of the form...
%       ts = nrec x 1 array of the timestamps that start each block.
%       cr = nrec x 512 array of the data
%       sFreq = sampling frequency


fss = regexp(fname, filesep, 'once');
if isempty(fss) % assume file is located in current directory
	fname = fullfile(pwd, fname);
end

try
	if nargin == 1
		[ts,cr,sFreq] = readCRcore(fname);  %  timestams ts are in 0.1 milliseconds units!!!!!
		readParent = true;
	elseif nargin == 3
		 [ts,cr,sFreq] = readCRcore(fname, nBlockMin, nBlockMax);  %  timestams ts are in 0.1 milliseconds units!!!!!
		 readParent = true;
	elseif nargin == 4
		[ts,cr,sFreq] = readCRcore(fname, nBlockMin, nBlockMax);  %  timestams ts are in 0.1 milliseconds units!!!!!
	end
catch err
	if strcmp(err.identifier, 'FPB:UtilsIO:EmptyNCSFile')
		error('readCRTsd:FileEmpty', 'Please provide another file name. This file is empty!')
	end
	if strcmp(err.identifier, 'FPB:UtilsIO:readCRcore:ErrorReadingNCSFile')
		error('readCRTsd:ErrorReadingNCSFile', 'There was an error openning the file. Perhaps non-existent.')
	end
end
dd=reshape(cr',1,length(cr(:)));
blockSize = 512;
nBlocks = size(cr,1);
dT = 10000/sFreq; % in tstamps

clear cr;
TIME = zeros(size(dd));
ts = [ts;ts(end) + 512*dT];     
for iBlock = 1:(nBlocks)
  %DisplayProgress(iBlock, nBlocks-1);
  TIME((blockSize * (iBlock-1) + 1):(blockSize * iBlock)) = ...
      linspace(ts(iBlock), ts(iBlock+1) - dT, blockSize);
end

if readParent
	H = getHeader(fname);
	if isempty(H)
		parentName = findParent(fname);
		if ~isempty(parentName)
			H = getHeader(parentName);
		else % parent doesn't exist; assume bitvolt is 6.1e-8
			warning('Parent not found. Assuming ADBitVolts = 6.1e-8')
			ADBitVolts = 6.1e-8;
			dd = dd*ADBitVolts;
			eeg = tsd.tsd(TIME', dd', 'ts');
			return
		end
	end
	i = 1;
	while isempty(strfind(H{i}, '-ADBitVolts'))
		i = i + 1;
	end
	bvInd = strfind(H{i}, '-ADBitVolts');
	ADBitVolts = str2double(H{i}(11+bvInd:end));
	if ~isempty(strfind(fname, 'rat100'))
		ADBitVolts = 3.0500e-08; % I know rat100 has this value
	end

	if isnan(ADBitVolts) % several values for '-ADBitVolts' in header
		% first find the first space
		for j = 12:length(H{i})
			if strcmp(H{i}(j), ' ')
				break;
			end
		end
		% now convert
		ADBitVolts = str2double(H{i}(12:j-1));
	end
	dd = dd*ADBitVolts;
end

eeg = tsd(TIME', dd', 'ts');