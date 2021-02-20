%% Main file to run examples %%
prompt = 'What example do you want to run? (for example 1, enter 1):';
inp = input(prompt);

switch inp
    case 1
        ibipp('Example pics\half_mbb.png',150,0.4)
    case 2
        ibipp('Example pics\hammerhead.png',200,0.5,[2 1 1 2],...
            [180 180 180 180],'optimization','BESO','filterRadius',3)
    case 3
        ibipp('Example pics\2point.png',300,0.45,[1 1],[0 0],'optimization',...
            'levelset','tau',5e-5,'preserveLoad',1,'preserveSupport',1)
    case 4
        ibipp('Example pics\half_spanner.png',250,0.5,'pressure',[1,1],...
            'preserveload',2,'preservesupport',1,'modelname','spanner.stl',...
            'symmetry','right','modeltype','extrude','extrudelength',0.2)
    otherwise 
        error('Error: Invalid input selection (select 1, 2, 3, or 4)')
end 