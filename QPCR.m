% QPCR analysis
% This class manages references and interpolating the references to get the concentrations of unknowns on the same plate
% It also provides error estimates using bootstrapping of the interpolation process
classdef QPCR < handle
  properties
    options;	% Option settings
    wellnames;	% Cell array of well names (e.g. wellnames(8,12) )
    ctgrid;	% Grid of plate ct values (e.g ctgrid(8,12) )
    primers;	% Cell array of primer pairs used
    refs;	% Map of primer->ref entry; each entry is a struct containing the reference information for that primer pair
  end

  methods(Static)
    function m=makemap(grid,samplewells)
      if nargin<2
        samplewells={};
      end
      m=containers.Map();
      for i=1:size(grid,1)
        for j=1:size(grid,2)
          wellname=sprintf('%c%d','A'+i-1,j);
          if isempty(samplewells) || ismember(wellname,samplewells)
            m(wellname)=grid(i,j);
          end
        end
      end
    end
  end
  
  methods
    function obj=QPCR(ctgrid,varargin)
      defaults=struct('extrapolate',true,'ctnoise',0.3,'nboot',100,'ci',80,'interpmethod','linear','fitrange',[8,20]);
      obj.options=processargs(defaults,varargin);

      obj.ctgrid=ctgrid;
      % Create a map to go from wellname to ct
      obj.primers=cell(size(ctgrid));
      obj.refs=containers.Map();
      obj.wellnames=cell(size(ctgrid));
      for i=1:size(ctgrid,1)
        for j=1:size(ctgrid,2)
          obj.wellnames{i,j}=sprintf('%c%d','A'+i-1,j);
        end
      end
    end
    
    function wlist=parsewells(obj,wells)
    % Parse a well input and return a cell array of individual well names
    % Well string is of form 'A1,C2,A3:B4,C5' or a cell array of well strings
    % Ranges increase row number most frequently (i.e. G3:B4 is G3,H3,A4,B4 )
      if isnumeric(wells)
        % Already a well list
        if any(wells<1 | wells>length(obj.ctgrid(:)))
          error('Bad well number in list [%s]\n',sprintf('%d ',wells));
        end
        wlist=wells;
      elseif iscell(wells)
        wlist=[];
        for i=1:length(wells)
          w=obj.parsewells(wells{i});
          wlist=[wlist,w];
        end
      elseif any(wells==',')
        % List
        comma=find(wells==',');
        wlist=[];
        prev=1;
        comma=[comma,length(wells)+1];
        for i=1:length(comma)
          w=obj.parsewells(wells(prev:comma(i)-1));
          wlist=[wlist,w];
          prev=comma(i)+1;
        end
      elseif any(wells==':')
        % Range
        colon=find(wells==':',1);
        first=obj.parsewells(wells(1:colon-1));
        last=obj.parsewells(wells(colon+1:end));
        if first>last || last>length(obj.ctgrid(:)) || first<1
          error('Bad well range %s\n', wells);
        end
        wlist=first:last;
      else
        if length(wells)<2 || wells(1)<'A' || wells(1)>'A'-1+size(obj.ctgrid,1) || str2double(wells(2:end))>size(obj.ctgrid,2)
          error('Bad well string: %s\n', wells);
        end
        row=wells(1)-'A'+1;
        col=str2double(wells(2:end));
        wlist=(col-1)*size(obj.ctgrid,1)+row;
      end
    end

    function ct=getct(obj,well)
      wlist=obj.parsewells(well);
      ct=obj.ctgrid(wlist);
    end
    
    function addref(obj,primer,refwells,refconcs,varargin)
      defaults=struct('units','nM');
      args=processargs(defaults,varargin);

      if isKey(obj.refs,primer)
        error('addref: Already have a reference defined for %s\n',primer);
      end
      refwlist=obj.parsewells(refwells);
      ct=obj.getct(refwells);
      if length(ct)~=length(refconcs)
        error('addref: wells and refconcs must have same length\n');
      end
      interpct=5:.05:30;
      interpdata=obj.bootcompute(ct,refconcs,interpct);
      % Fit points with ct in fitrange to a straight line
      ctsel=ct>=obj.options.fitrange(1) & ct<=obj.options.fitrange(2);
      minconc=min(refconcs(ctsel));
      maxconc=max(refconcs(ctsel));
      concsel=refconcs>=minconc & refconcs<=maxconc & isfinite(ct);
      obj.refs(primer)=struct('name',primer,'wells',refwlist,'welldescr',refwells,'ct',ct,'concs',refconcs,'interpdata',interpdata,'units',args.units);
      if sum(concsel)>=2
        fit=polyfit(log(refconcs(concsel)),ct(concsel),1);
        if any(~isfinite(fit))
          keyboard;
        end
        fitconcrange=[minconc,maxconc];
        eff=exp(-1/fit(1));
        ct0=exp(-fit(2)/fit(1));
        ct10=ct0*eff^-10;
        predict=polyval(fit,log(refconcs));
        ctnoise=sqrt(sum((ct(concsel)-predict(concsel)).^2)/(sum(concsel)-2));  % d.f.=N-2 since fit has 2 parameters
        deviation=nan(size(refconcs));
        deviation(concsel)=ct(concsel)-predict(concsel);
        fprintf('Primer %s model:  efficiency=%.2f, Conc(Ct=10)=%.1f%s, Conc(Ct=0)=%.2g%s\n', primer, eff, ct10, args.units,ct0, args.units);
        tmp=obj.refs(primer);
        tmp.mdl=struct('fit',fit,'concrange',fitconcrange,'eff',eff,'ct0',ct0,'ct10',ct10,'ctnoise',ctnoise,'deviation',deviation);
        obj.refs(primer)=tmp;
      end
      obj.setref(primer,refwells);
    end
    
    function dupref(obj,existing,new)
      obj.refs(new)=obj.refs(existing);
    end
    
    function ctnoise=estimateNoise(obj)
    % Check if we have reference replicates; use to compute a ct noise level
      ctnoise=nan;
      if length(unique(refconcs))< length(refconcs)
        stdlist=[];
        uy=unique(refconcs);
        for i=1:length(uy)
          sel=isfinite(x)&refconcs==uy(i);
          if sum(sel)>1
            stdlist(end+1)=std(x(sel));
          end
        end
        if length(stdlist)>0
          fprintf('Ct stdev = %.2f over %d replicated concentrations\n',mean(stdlist), length(stdlist));
        end
        if length(stdlist)>=3 && isempty(obj.options.ctnoise)
          ctnoise=mean(stdlist);
        end
      end
    end

    function setref(obj,ref,wells)
    % Set the reference used for the given wells
      if ~isKey(obj.refs,ref)
        error('setref: Reference %s has not been defined using addref()\n', ref);
      end
      wlist=obj.parsewells(wells);
      for i=1:length(wlist)
        if ~isempty(obj.primers{wlist(i)}) && ~strcmp(obj.primers{wlist(i)},ref)
          error('setref: Attempted to set reference for well %s to %s, but was previously set to reference %s\n', obj.wellnames{wlist(i)}, ref, obj.primers{wlist(i)});
        end
        obj.primers{wlist(i)}=ref;
      end
    end
    
    function [conc,cilow,cihigh]=getconc(obj,ref,rep1,rep2,rep3)
      w1=obj.parsewells(rep1);
      wells=w1;
      if nargin>3
        w2=obj.parsewells(rep2);
        if length(w2)~=length(w1)
          error('Replicate lists must all have the same length');
        end
        wells=[wells,w2];
      end
      if nargin>4
        w3=obj.parsewells(rep3);
        if length(w3)~=length(w1)
          error('Replicate lists must all have the same length');
        end
        wells=[wells,w3];
      end
      
      obj.setref(ref,wells); 	% Set it for error-checking and subsequent plotting
      refdata=obj.refs(ref);
      ct=obj.getct(wells);
      conc=interp1([refdata.interpdata.ct],[refdata.interpdata.conc],ct);
      cilow=interp1([refdata.interpdata.ct],[refdata.interpdata.cilow],ct);
      cihigh=interp1([refdata.interpdata.ct],[refdata.interpdata.cihigh],ct);
      conc=reshape(conc,[],length(w1));
      cilow=reshape(cilow,[],length(w1));
      cihigh=reshape(cihigh,[],length(w1));
      conc=mean(conc,1);
      cilow=mean(cilow,1);
      cihigh=mean(cihigh,1);
      range=max(conc)./min(conc);
      if any(range)>2
        fprintf('Warning: Some concentrations vary by %.1fx between replicates over wells %s\n',max(range),sprintf('%s ',wells));
      end
    end

    function plot(obj)
      refnames=sort(keys(obj.refs));
      nx=ceil(length(refnames)/4);
      ny=ceil(length(refnames)/nx);
      clf;
      bnds=containers.Map();
      for i=1:length(refnames)
        r=obj.refs(refnames{i});
        subplot(ny,nx,i);
        obj.plotref(refnames{i});
        c=axis;
        if ~isKey(bnds,r.units)
          bnds(r.units)=c(1:2);
        else
          b=bnds(r.units);
          bnds(r.units)=[min(b(1),c(1)),max(b(2),c(2))];
        end
      end
      % Align all the x axes
      for i=1:length(refnames)
        r=obj.refs(refnames{i});
        subplot(ny,nx,i);
        c=axis;
        c(1:2)=bnds(r.units);
        axis(c);
      end
    end
    
    function plotref(obj,ref,sampdata)
    % Plot reference information on conc vs. ct plot
      ti=ref;	% Title
      r=obj.refs(ref);
      h=[]; leg={};
      h(end+1)=semilogx(r.concs,r.ct,'or');	% References
      leg{end+1}='Reference';
      
      % Draw lines showing best point, confidence intervals
      hold on;
      plot(smooth([r.interpdata.conc]),[r.interpdata.ct],'g');
      h(end+1)=plot(smooth([r.interpdata.cilow]),[r.interpdata.ct],'g:');
      leg{end+1}='Reference CI';
      plot(smooth([r.interpdata.cihigh]),[r.interpdata.ct],'g:');

      % Draw fit
      if isfield(r,'mdl')
        fitconc=[r.interpdata.conc];
        sel=fitconc>=r.mdl.concrange(1) & fitconc<=r.mdl.concrange(2);
        fitconc=fitconc(sel);
        fitct=polyval(r.mdl.fit,log(fitconc));
        h(end+1)=plot(fitconc,fitct,'r');
        leg{end+1}='Linear Fit';
        if isfield(r.mdl,'eff') && isfield(r.mdl,'ct10')
          ti=[ti,sprintf(' eff=%.2f, Conc(Ct=10)=%.2f%s, ctnoise=%.2f', r.mdl.eff, r.mdl.ct10,r.units, r.mdl.ctnoise)];
        end
      end
      
      % Overlay interpolated data for samples that use this reference
      sampwells=find(strcmp(obj.primers(:),ref));
      isref=ismember(sampwells,r.wells);	% Don't replot references
      sampwells=sampwells(~isref);
      sampcts=obj.getct(sampwells);

      if isempty(sampcts)
        fprintf('No non-reference data points for %s to plot\n', ref);
      else
        [mid,low,high]=obj.getconc(ref,sampwells);
        sampdata=obj.bootcompute(r.ct,r.concs,sampcts);
        h(end+1)=plot(mid,sampcts,'xb');
        leg{end+1}='Data';
        for i=1:length(mid)
          hh=plot([low(i),high(i)],sampcts(i)*[1,1],'-b');
          text(high(i),sampcts(i),obj.wellnames{sampwells(i)});
          if i==1
            h(end+1)=hh;
            leg{end+1}='Data CI';
          end
        end
      end

      c=axis;
      c(1)=max(c(1),min(r.concs)/10);
      c(2)=min(c(2),max(r.concs)*10);
      axis(c);

      legend(h,leg);
      title(ti);
      xlabel(sprintf('Concentration (%s)',r.units));
      ylabel('Ct');
    end
    
    function plotdeviations(obj)
    % Plot model deviations to diagnose poor dilution implementation
      clf;
      refnames=sort(keys(obj.refs));
      devs=[];
      for i=1:length(refnames)
        d=obj.refs(refnames{i}).mdl.deviation;
        if i==1
          devs=d;
        else
          devs(i,1:length(d))=obj.refs(refnames{i}).mdl.deviation;
        end
      end
      bar(devs');
      legend(refnames);
      xlabel('Dilution step');
      ylabel('Ct difference from linear fit');
    end
    
    function data=bootcompute(obj,x,y,xv)
    % Bootstrap estimate
      bootstat=bootstrp(obj.options.nboot,@(xa,ya)(obj.compute(xa+randn(size(xa))*obj.options.ctnoise,ya,xv+randn(size(xv))*obj.options.ctnoise)),x,y);
      fracfinite=mean(isfinite(bootstat));
      %bootstat(:,fracfinite<0.5)=nan;		% Make sure that points that are often unusable are not given a value
      med=nanmedian(bootstat);
      med(fracfinite<0.5)=nan;
      data=struct('ct',num2cell(xv(:)'),'conc',num2cell(med),'std',num2cell(nanstd(bootstat)),'cilow',num2cell(prctile(bootstat,(100-obj.options.ci)/2)),'cihigh',num2cell(prctile(bootstat,100-(100-obj.options.ci)/2)),'ci',obj.options.ci,'bootstat',num2cell(bootstat,1),'fracfinite',num2cell(fracfinite));
      data=reshape(data,size(xv));
    end
    
    function conc=compute(obj,x,y,xv)
      sel=isfinite(x) & isfinite(y);
      x=x(sel);y=y(sel);
      uy=sort(unique(y(:)'));
      ux=[];
      for i=1:length(uy)
        sel=(y==uy(i) & isfinite(x));
        ux(i)=nanmean(x(sel));
      end
      while any(ux(2:end)>ux(1:end-1))
        for i=1:length(ux)-1
          if ux(i+1)>ux(i)
            %fprintf('Fixing inverted: ');
            ux=[ux(1:i-1),mean(ux(i:i+1)),ux(i+2:end)];
            uy=[uy(1:i-1),mean(uy(i:i+1)),uy(i+2:end)];
            break;
          end
        end
      end

      if length(ux)<2
        fprintf('compute: Only %d interpolation points\n', sum(sel));
        conc=nan;
        return;
      end

      if obj.options.extrapolate
        conc=exp(interp1(ux,log(uy),xv(:),obj.options.interpmethod,'extrap'));
      else
        conc=exp(interp1(ux,log(uy),xv(:),obj.options.interpmethod,NaN));
      end

      conc=reshape(conc,size(xv));
    end

  end % methods
end % classdef

function y=smooth(x)
  y=medfilt1(x,5);
end

