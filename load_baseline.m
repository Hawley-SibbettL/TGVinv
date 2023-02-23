% function loads the original image used in the synthetic inversion, and
% plots the resulting image
function [] = load_baseline(filename, plotflag)

folder = ['C:\Users\', getenv('username'), '\OneDrive - The University of Nottingham\UoN Box Migration\Data\New TGV test models'];

load(fullfile(folder, [filename,'.mat']))


end