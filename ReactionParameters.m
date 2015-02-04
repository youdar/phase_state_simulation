function [klpf,klmf,klpm,klmm,kap,kam,kn,cnmRatio,kbf,kbm,y00,p00,rs,rfm,n,m,vScale] = ReactionParameters(DataArray)

    klpf = DataArray(1); 
    klmf = DataArray(2);
    klpm = DataArray(3);
    klmm = DataArray(4);
    kap = DataArray(5);
    kam = DataArray(6);
    kn = DataArray(7);
    cnmRatio = DataArray(8);        % cnmRatio = (cnm equivalent)/cn = cnm*p0/cn
    kbf = DataArray(9);
    kbm = DataArray(10);
    y00 = DataArray(11);
    p00 = DataArray(12);
    rs = DataArray(13);
    rfm = DataArray(14);
    n = DataArray(15);
    m = DataArray(16);
    vScale = DataArray(17);              % system scaling factor - used to reduce particle number

end

