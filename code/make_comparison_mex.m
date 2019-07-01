function make_comparison_mex
%% Attempt to make FLPGM (a.k.a. LPMRF) model code
if(~exist(sprintf('+mrfs/+samplers/ais_mex.%s',mexext()),'file'))
    fprintf('Attempting to make the mex files for the FLPGM model\n');
    cd('mexcode');
    try
        makemex();
    catch e
        warning('Error making the mex files for the FLPGM (a.k.a. LPMRF) model (proceeding without these files) \n');
        disp(e.message);
    end
    cd('..');
end
%% Attempt to amke QUIC code
if(~exist(sprintf('QUIC/QUIC.%s',mexext()),'file'))
    cd('QUIC');
    try
        mex -llapack QUIC.C QUIC-mex.C -output QUIC
    catch
        warning('Error making the mex files for the Gaussian GM (a.k.a. QUIC files) model (proceeding without these files)\n');
    end
    cd('..');
end