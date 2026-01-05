@echo off
matlab -nodesktop -r "currentDir = pwd; repo_path = fullfile(currentDir,'./NetworkGen/'); addpath(genpath(repo_path)); savepath;"