function [ S_old, step, VPI, oilrecovery, cumulateoil, watercut, dt_ref, countime, time ] = ...
         condinicial( S_old, step, VPI, oilrecovery, cumulateoil, watercut, dt_ref, countime, time, metsat )
%
global caso malha formethod

if exist(sprintf('%s\\%s\\%s\\%s\\Results\\Saturation.mat',caso,malha,formethod,metsat),'file') ~= 0
    load(sprintf('%s\\%s\\%s\\%s\\Results\\Saturation',caso,malha,formethod,metsat));
end

if exist(sprintf('%s\\%s\\%s\\%s\\Results\\VPI.mat',caso,malha,formethod,metsat),'file') ~= 0
    load(sprintf('%s\\%s\\%s\\%s\\Results\\VPI',caso,malha,formethod,metsat));
end

if exist(sprintf('%s\\%s\\%s\\%s\\Results\\Countime.mat',caso,malha,formethod,metsat),'file') ~= 0
    load(sprintf('%s\\%s\\%s\\%s\\Results\\Countime',caso,malha,formethod,metsat));
end

if exist(sprintf('%s\\%s\\%s\\%s\\Results\\Oilrecovery.mat',caso,malha,formethod,metsat),'file') ~= 0
    load(sprintf('%s\\%s\\%s\\%s\\Results\\Oilrecovery',caso,malha,formethod,metsat));
end

if exist(sprintf('%s\\%s\\%s\\%s\\Results\\Cumulateoil.mat',caso,malha,formethod,metsat),'file') ~= 0
    load(sprintf('%s\\%s\\%s\\%s\\Results\\Cumulateoil',caso,malha,formethod,metsat));
end

if exist(sprintf('%s\\%s\\%s\\%s\\Results\\Watercut.mat',caso,malha,formethod,metsat),'file') ~= 0
    load(sprintf('%s\\%s\\%s\\%s\\Results\\Watercut',caso,malha,formethod,metsat));
end

if exist(sprintf('%s\\%s\\%s\\%s\\Results\\TimeStep.mat',caso,malha,formethod,metsat),'file') ~= 0
    load(sprintf('%s\\%s\\%s\\%s\\Results\\TimeStep',caso,malha,formethod,metsat));
end

if exist(sprintf('%s\\%s\\%s\\%s\\Results\\DT.mat',caso,malha,formethod,metsat),'file') ~= 0
    load(sprintf('%s\\%s\\%s\\%s\\Results\\DT',caso,malha,formethod,metsat));
end

if exist(sprintf('%s\\%s\\%s\\%s\\Results\\TIME.mat',caso,malha,formethod,metsat),'file') ~= 0
    load(sprintf('%s\\%s\\%s\\%s\\Results\\TIME',caso,malha,formethod,metsat));
end

end

