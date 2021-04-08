
clear,clc;

%%
foldername = 'D:\MATLABWorkplace\TIGRE-master\MATLAB\2020-06-12_231835\Calibrations\SC-100KV-Bowtie-Full\Factory';
srcfilename = [foldername filesep 'Calibration.xml'];
tmpfilename = [foldername filesep 'tmp.xml'];
copyfile(srcfilename,tmpfilename);

%%
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
fclose(fidsrc)
fclose(fidtmp)


%%
calib = xml2struct(tmpfilename);
delete(tmpfilename);

