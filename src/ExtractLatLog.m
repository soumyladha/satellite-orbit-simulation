file=('latlog.txt'); % File Path
    format=('%f %f'); % specifiing input format
    fid=fopen(file); % FID syntax
    celltemp = textscan(fid,format,'headerlines',1); % importing as cell
    latlog = cell2mat(celltemp); 
fclose(fid);
