function [vineMat, familyMat, parMat, par2Mat, USample] = vinecopulawrapper( U, nSamples, vineType, indTest, pLevel, nWorkers, familySet)
% Defaults from VineCopula documentation
if(nargin < 3); vineType = 0; end
if(nargin < 4); indTest = 0; end
if(nargin < 5); pLevel = 0.05; end 
if(nargin < 6); nWorkers = 1; end
if(nargin < 7); familySet = []; end

%% Save files
argsFile = tempname;
dataFile = tempname;
paramsFile = tempname;
sampleFile = tempname;

fid = fopen(argsFile,'w+');
fprintf( fid, '%d,%d,%d,%g,%d', nSamples, vineType, indTest, pLevel, nWorkers );
if(~isempty(familySet))
    fprintf( fid, ',%d', familySet );
end
fprintf( fid, '\n' );
fclose(fid);

fid = fopen(dataFile,'w+');
fprintf( fid, ['%d', repmat(',%d', 1, full(size(U,1))-1), '\n'], full(U(:)));
fclose(fid);
fprintf('Finished writing files for vinecopulawrapper.m\n');

%% Call R script
[status,result] = system('hostname -f');
if(status == 0 && ~isempty(strfind(result,'cs.utexas.edu')))
    Rcmd = '/lusr/opt/R-3.2.2/bin/Rscript';
else
    Rcmd = 'Rscript';
end
cmd = sprintf('cd comparison && %s --vanilla vine.copula.wrapper.R %s %s %s %s', ...
    Rcmd, argsFile, dataFile, paramsFile, sampleFile);
fprintf('VineCopulaWrapper.m: Executing the following system command:\n%s\n\n', cmd);
for attempt = 1:10
    status = system(cmd);
    if(status == 0)
        fprintf('Correctly executed command on attempt = %d\n',attempt);
        break;
    end
end

%% Extract parameters
params = csvread(paramsFile);
p = size(U,2);
assert(size(params,1) == 4*p, 'params have 4*p rows');
assert(size(params,2) == p, 'params have p cols');
vineMat = params(1:p,:);
familyMat = params((p+1):(2*p),:);
parMat = params((2*p+1):(3*p),:);
par2Mat = params((3*p+1):(4*p),:);

%% Extract samples
USample = csvread(sampleFile);

end
