sf.name = 'gray_scott';
if ~isfield(sf,'caption')
    sf.caption='TODO';
end
curr_dir = cd;
try
    texfig(sf,cd,'gray_scott',gcf,p)
catch me
    cd(curr_dir);
    texfig(sf,cd,'gray_scott',gcf)
end
cd(curr_dir);