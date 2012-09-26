function report = wellCalculateProduction(rSol,wSol,W,fluid,time)
%calculate well data for each well
%
% SYNOPSIS:
% report = wellCalculateProduction(rSol,wSol,W,fluid,time)
% 
% PARAMETERS:
%   rSol   - stuct with cell solutions
%   wSol   - stuct with well solutions
%   W      - well stuct
%   fluid  - fluid struct
%
% RETURNS:
%   report - return quantities of each well with eclipse keyword
%
    report.TIME = time;
    report.WVPT = cellfun(@(f) sum(f(:)),{wSol.flux})';
    report.WBHP = [wSol.pressure]';    
    fw = cellfun(@(c) fluid.fw(struct('s', rSol.s(c))), { W.cells }, 'UniformOutput', false);
    report.WWPR = cellfun(@(f,w) sum(f(:).*w),{wSol.flux}, fw)';
    report.WOPR = cellfun(@(f,w) sum(f(:).*(1-w)),{wSol.flux}, fw)';
    %report.WOPR = cellfun(@(f,w) sum(f(:).*(1-fluid.fw(w))),{wSol.flux},wSat)';
    report.WWCT = report.WWPR./report.WVPT;
    %report.WWPT = cellfun(@(f,cells) sum(f(:).*rSol.s(cells(:))),{wSol.flux},{W.cells});
%     for i=1:numel(wsat)
%         wSat{end+1} = struct('s',wsat{i});
%     end
%     wsat = cellfun(@(cells) rSol.s(cells),{W.cells},'UniformOutput',false);
%     wSat={};