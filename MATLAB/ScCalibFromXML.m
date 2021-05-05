function ScCalibXML = ScCalibFromXML(datafolder, ScanXML)
%% Load Calibration.xml for Scatter Correction
% Reference: Improved scatter correction using adaptive scatter kernel superposition
% Date: 2021-05-05
% Author: Yi Du (yi.du@hotmail.com)

foldername = [datafolder filesep 'Calibrations' filesep 'SC-' ScanXML.Acquisitions.Voltage.Text 'KV-*' filesep 'Factory'];
srcdir = dir([foldername filesep 'Calibration.xml']);
srcfilename = [srcdir.folder filesep 'Calibration.xml'];
tmpfilename = [srcdir.folder filesep 'tmp.xml'];
copyfile(srcfilename,tmpfilename);

%% Delete redundent comment line
fidsrc = fopen(srcfilename);
fidtmp = fopen(tmpfilename, 'wt');
count = 0;
while(true)
    tmp = fgetl(fidsrc);
    % End-of-File
    if(tmp == -1)
        break;
    end
    
    if(contains(tmp, '!--'))
        continue;
    end
    
    fprintf(fidtmp, '%s\n', tmp);
    count = count+1;
end
fclose(fidsrc);
fclose(fidtmp);

%% Export as struct
tmp = xml2struct(tmpfilename);
ScCalibXML = tmp.Calibration;
% delete temperary file
delete(tmpfilename);

end
