% QPCR analysis
% This class manages references and interpolating the references to get the concentrations of unknowns on the same plate
classdef QPCR < handle
  properties
    options;	% Option settings
    wellnames;	% Cell array of well names (e.g. wellnames(8,12) )
    ctgrid;	% Grid of plate ct values (e.g ctgrid(8,12) )
    dilgrid;	% Dilutions from samples
    lengrid;	% Length of products
    strandgrid;	% Number of strands in samples (1 or 2)
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

    function [c,clow,chigh]=concsub(c1,c1low,c1high,c2,c2low,c2high)
    % Subtract concentrations including confidence intervals
      c=max(0,c1-c2);
      c(~isfinite(c1) | ~isfinite(c2))=nan;
      clow=max(0,c-sqrt((c1low-c1).^2+(c2low-c2).^2));
      chigh=c+sqrt((c1high-c1).^2+(c2high-c2).^2);
    end
    
    function [c,clow,chigh]=concadd(c1,c1low,c1high,c2,c2low,c2high)
    % Add concentrations including confidence intervals
      c=c1+c2;
      clow=c-sqrt((c1low-c1).^2+(c2low-c2).^2);
      chigh=c+sqrt((c1high-c1).^2+(c2high-c2).^2);
    end
    

  end
  
  methods
    function obj=QPCR(ctgrid,varargin)
      defaults=struct('extrapolate',true,'ci',80,'minct',7,'dilgrid',[],'lengrid',[],'strandgrid',[],'plate2',[]);
      args=processargs(defaults,varargin);

      obj.options=rmfield(args,{'dilgrid','lengrid','strandgrid'});
      if ~isempty(args.plate2)
        % Append a second plate
        ctgrid=[ctgrid,args.plate2];
      end
      obj.ctgrid=ctgrid;
      if isempty(args.dilgrid)
        obj.dilgrid=ones(size(ctgrid));
      else
        assert(all(size(ctgrid)==size(args.dilgrid)));
        obj.dilgrid=args.dilgrid;
      end
      if isempty(args.lengrid)
        obj.lengrid=ones(size(ctgrid));
      else
        assert(all(size(ctgrid)==size(args.lengrid)));
        obj.lengrid=args.lengrid;
      end
      if isempty(args.strandgrid)
        obj.strandgrid=ones(size(ctgrid));
      else
        assert(all(size(ctgrid)==size(args.strandgrid)));
        assert(all(args.strandgrid(:)==1 | args.strandgrid(:)==2 | isnan(args.strandgrid(:))));
        obj.strandgrid=args.strandgrid;
      end
      
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
        for i=1:length(wells(:))
          w=obj.parsewells(wells{i});
          wlist=[wlist,w];
        end
        wlist=reshape(wlist,size(wells));
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

    function nms=getwellnames(obj,wlist) 
    % Convert a list of well numbers to names
      nms=reshape(obj.wellnames(wlist(:)),size(wlist));
    end
    
    function ct=getct(obj,well)
      wlist=obj.parsewells(well);
      ct=obj.ctgrid(wlist);
    end
    
    function addmanualref(obj,primer,efficiency,ct10,varargin)
      defaults=struct('units','nM','length',1,'strands',2,'ctlod',20,'dilution',1);
      args=processargs(defaults,varargin);
      obj.refs(primer)=struct('name',primer,'wells',[],'welldescr',{{}},'ct',10,'concs',ct10,'units',args.units,'ctwater',[],'samples',containers.Map(),'len',args.length,'dilution',args.dilution,'strands',args.strands);

      ct0=ct10*efficiency^10;
      fit(1)=-1/log(efficiency);
      fit(2)=-log(ct0)*fit(1);

      lod=ct0*efficiency^-args.ctlod;
      fprintf('Primer %4s model:  efficiency=%4.2f, Conc(Ct=10)=%.3g%s, Det Limit=%6.2g%s (@Ct=%.1f)\n', primer, efficiency, ct10, args.units,lod,args.units,args.ctlod);
      tmp=obj.refs(primer);
      tmp.mdl=struct('N',0, 'fit',fit,'conc',[],'ct',[],'concrange',[nan,nan],'eff',efficiency,'ct0',ct0,'ct10',ct10,'ctnoise',nan,'deviation',[],'ctlod',args.ctlod,'lod',lod,'sxloglod',nan);
      obj.refs(primer)=tmp;
    end

    function addref(obj,primer,refwells,refconcs,varargin)
      defaults=struct('units','nM','length',[],'strands',[],'efficiency',[]);
      args=processargs(defaults,varargin);

      if isKey(obj.refs,primer)
        error('addref: Already have a reference defined for %s\n',primer);
      end
      refwlist=obj.parsewells(refwells);
      dilution=unique(obj.dilgrid(refwlist));
      if length(dilution)>1
        error('Reference %s has multiple different dilutions specified by dilgrid\n',primer);
      end
      if ~isempty(args.strands)
        obj.strandgrid(refwlist)=args.strands;
      end
      strands=unique(obj.strandgrid(refwlist));
      if length(strands)>1
        error('Reference %s has multiple different strand counts specified by strandgrid\n',primer);
      end
      if ~isempty(args.length)
        obj.lengrid(refwlist)=args.length;
      end
      len=unique(obj.lengrid(refwlist));
      if length(len)>1
        error('Reference %s has multiple different lengths specified by lengrid\n',primer);
      end
      
      ct=obj.getct(refwells);
      ct=ct(:);
      refconcs=refconcs(:);
      if length(ct)~=length(refconcs)
        error('addref: wells and refconcs must have same length\n');
      end

      uconcs=unique(refconcs);
      wnames=obj.getwellnames(refwlist);
      for i=1:length(uconcs)
        ct1=ct(refconcs==uconcs(i));
        rw1=wnames(refconcs==uconcs(i));
        if range(ct1)>1.0
          fprintf('Warning: Cts vary by %.1f between reference replicates over wells: ',range(ct1));
          for j=1:length(ct1)
            fprintf('%s [%.1f] ',rw1{j},ct1(j));
          end
          fprintf('\n');
        end
      end

      % Fit points with ct in fitrange to a straight line
      ctwater=ct(refconcs==0);

      % Compute limit of detection
      if isempty(ctwater)
        ctlod=max([ct;19])+1;
        fprintf('Warning: Missing "water" sample(s) for primer %s -- assuming maximum CT of %.1f\n',primer,ctlod);
      elseif all(isnan(ctwater))
        ctlod=max([ct;19])+1;
        fprintf('Warning: No Ct returned for "water" sample(s) for primer %s -- assuming maximum CT of %.1f\n',primer,ctlod);
      else
        ctlod=min(ctwater);
      end

      ctsel=ct>=obj.options.minct & ct<=(ctlod-2.5);    % This much away from the noise floor will add 0.27 Ct of "noise" of the measurement, which is on the order of the usual noise
      minconc=min(refconcs(ctsel&refconcs>0));
      maxconc=max(refconcs(ctsel));
      if isempty(minconc)
        fprintf('No Ct values for any concentrations of reference for primer %s\n', primer);
        return;
      end
      concsel=refconcs>=minconc & refconcs<=maxconc & isfinite(ct);
      obj.refs(primer)=struct('name',primer,'wells',refwlist(refconcs>0),'welldescr',{obj.getwellnames(refwlist(refconcs>0))},'ct',ct(refconcs>0),'concs',refconcs(refconcs>0),'units',args.units,'ctwater',ctwater,'samples',containers.Map(),'len',len,'dilution',dilution,'strands',strands);
      % samples will hold individual sample data
      if sum(concsel)>=1
        if ~isempty(args.efficiency)
          fprintf('Forced efficiency to %.2f\n', args.efficiency);
          fit(1)=-1/log(args.efficiency);
          fit(2)=mean(ct(concsel)-fit(1)*log(refconcs(concsel)));
        elseif sum(concsel)==1
          fprintf('Only 1 concentration point, assuming perfect efficiency\n');
          fit=polyfit(log(refconcs(concsel)*[1,0.5]),ct(concsel)+[0,1],1);	% Assume 2x
        else
          fit=polyfit(log(refconcs(concsel)),ct(concsel),1);
        end
        if any(~isfinite(fit))
          keyboard;
        end
        fitconcrange=[minconc,maxconc];
        eff=exp(-1/fit(1));
        ct0=exp(-fit(2)/fit(1));
        ct10=ct0*eff^-10;
        predict=polyval(fit,log(refconcs));
        ctnoise=sqrt(median((ct(concsel)-predict(concsel)).^2)*sum(concsel)/(sum(concsel)-2));  % d.f.=N-2 since fit has 2 parameters  (this is also known as Sy|x )
        % ctnoise=sqrt(sum((ct(concsel)-predict(concsel)).^2)/(sum(concsel)-2));  % d.f.=N-2 since fit has 2 parameters  (this is also known as Sy|x )
        deviation=nan(size(refconcs));
        deviation(concsel)=ct(concsel)-predict(concsel);
        lod=ct0*eff^-ctlod;
        sxloglod=ctnoise/eff*sqrt(1+1/sum(concsel)+(ctlod-mean(ct(concsel))).^2/(eff^2*sum((log(refconcs(concsel))-mean(log(refconcs(concsel)))).^2)));
        fprintf('Primer %s model:  efficiency=%.2f, Conc(Ct=10)=%.1f%s, Conc(Ct=0)=%.2g%s Det Limit=%.2g%s (Ct=%.1f)\n', primer, eff, ct10, args.units,ct0, args.units,lod,args.units,ctlod);
        tmp=obj.refs(primer);
        tmp.mdl=struct('N',sum(concsel), 'fit',fit,'conc',refconcs(concsel),'ct',ct(concsel),'concrange',fitconcrange,'eff',eff,'ct0',ct0,'ct10',ct10,'ctnoise',ctnoise,'deviation',deviation,'ctlod',ctlod,'lod',lod,'sxloglod',sxloglod);
        obj.refs(primer)=tmp;
      end
      obj.setref(primer,refwells);
    end
    
    function f=fitref(obj,primer)
    % Fit a model of the form C+Cbg = C0 * eff^-Ct
      ref=obj.refs(primer);
      %ft=fittype('(C0-log(x+CB))/log(EFF)','coefficients',{'EFF','C0','CB'});
      %      ft=fittype('C0*EFF.^-x-CB','coefficients',{'EFF','C0','CB'});
      ft=fittype('log(abs(exp(C0)*EFF.^-x-exp(CB)))','coefficients',{'EFF','C0','CB'});
      options=fitoptions(ft);
      %      options.StartPoint=[2,nanmean(ref.concs),nanmin(ref.concs)];
      options.StartPoint=[2,log(nanmean(ref.concs)),log(nanmin(ref.concs))];
      options.Lower=[1.5,-Inf,-Inf];
      options.Upper=[2.5,Inf,Inf];
      options.Display='iter';
      options.Algorithm='Trust-Region';
      options.TolFun=1e-50;
      options.TolX=1e-10;
      if false
        f=fit(ref.concs,ref.ct,ft,options);
        setfig(['fitref-',primer]);clf;
        plot(f,ref.concs,ref.ct);
        hold on;
        plot(f,'predfun');
        set(gca,'XScale','log');
        xlabel('Concentration');
        ylabel('Ct');
      else
        f=fit(ref.ct,log(ref.concs),ft,options);
        setfig(['fitref-',primer]);clf;
        plot(f,ref.ct,log(ref.concs));
        hold on;
        plot(f,'predfun');
        ylabel('Concentration');
        xlabel('Ct');
      end
    end
    
    function dupref(obj,existing,new)
      obj.refs(new)=obj.refs(existing);
    end
    
    function setref(obj,ref,wells)
    % Set the reference used for the given wells
      if ~isKey(obj.refs,ref)
        error('setref: Reference %s has not been defined using addref()\n', ref);
      end
      wlist=obj.parsewells(wells(:));
      for i=1:length(wlist)
        if ~isempty(obj.primers{wlist(i)}) && ~strcmp(obj.primers{wlist(i)},ref)
          error('setref: Attempted to set reference for well %s to %s, but was previously set to reference %s\n', obj.wellnames{wlist(i)}, ref, obj.primers{wlist(i)});
        end
        obj.primers{wlist(i)}=ref;
      end
    end
    
    function cmat=getconcmat(obj,ref,rep1,rep2,rep3,rep4,varargin)
      defaults=struct('length',nan,'dilution',nan,'strands',[]);
      args=processargs(defaults,varargin);
      if nargin<4; rep2=[]; end;
      if nargin<5; rep3=[]; end;
      if nargin<6; rep4=[]; end;
      [c,l,h]=getconc(obj,ref,rep1,rep2,rep3,rep4,'length',args.length,'dilution',args.dilution,'strands',args.strands);
      cmat=[c(:),l(:),h(:)];
    end
    
    function [conc,cilow,cihigh]=getconc(obj,ref,rep1,rep2,rep3,rep4,varargin)
      defaults=struct('length',nan,'dilution',[],'strands',[]);
      args=processargs(defaults,varargin);
      w1=obj.parsewells(rep1);
      wells=w1(:);
      if nargin>3 && ~isempty(rep2)
        w2=obj.parsewells(rep2);
        if any(size(w2)~=size(w1))
          error('Replicate lists must all have the same length');
        end
        wells=[wells,w2(:)];
      end
      if nargin>4 && ~isempty(rep3)
        w3=obj.parsewells(rep3);
        if any(size(w3)~=size(w1))
          error('Replicate lists must all have the same length');
        end
        wells=[wells,w3(:)];
      end
      if nargin>5 && ~isempty(rep4)
        w4=obj.parsewells(rep4);
        if any(size(w4)~=size(w1))
          error('Replicate lists must all have the same length');
        end
        wells=[wells,w4(:)];
      end
      
      obj.setref(ref,wells); 	% Set it for error-checking and subsequent plotting
      ct=obj.getct(wells);

      for i=1:size(ct,1)
        if range(ct(i,:))>1.0
          fprintf('Warning: Cts vary by %.1f between replicates over wells: ',range(ct(i,:)));
          for j=1:size(ct,2)
            fprintf('%s [%.1f] ',obj.wellnames{wells(i,j)},ct(i,j));
          end
          fprintf('\n');
        end
      end
      [conc,cilow,cihigh]=obj.getconcfromct(ref,ct);

      s=obj.refs(ref).samples;
      for i=1:length(conc)
        sel=isfinite(ct(i,:));
        if sum(sel)==0
          continue;
        end
        wellnms=obj.wellnames(wells(i,sel));
        sampname=wellnms{1};
        for j=2:length(wellnms)
          sampname=[sampname,'+',wellnms{j}];
        end
        s(sampname)=struct('name',sampname,'wells',{wellnms},'ct',ct(i,sel),'conc',conc(i),'cilow',cilow(i),'cihigh',cihigh(i));
      end
      % Correct for length, dilution, strands
      r=obj.refs(ref);
      w1=w1(:);
      if ~isempty(args.dilution)
        obj.dilgrid(w1)=args.dilution;
      end
      if ~isnan(args.length)
        obj.lengrid(w1)=args.length;
      end
      conc=conc.*(obj.dilgrid(w1)/r.dilution)./(obj.lengrid(w1)/r.len)./(obj.strandgrid(w1)/r.strands);
      cilow=cilow.*(obj.dilgrid(w1)/r.dilution)./(obj.lengrid(w1)/r.len)./(obj.strandgrid(w1)/r.strands);
      cihigh=cihigh.*(obj.dilgrid(w1)/r.dilution)./(obj.lengrid(w1)/r.len)./(obj.strandgrid(w1)/r.strands);
      conc=reshape(conc,size(w1));
      cilow=reshape(cilow,size(w1));
      cihigh=reshape(cihigh,size(w1));
    end
    
    function [conc,cilow,cihigh]=getconcfromct(obj,ref,ct)
      ctmean=nanmean(ct,2);
      ctcnt=sum(isfinite(ct),2);
      m=obj.refs(ref).mdl;

      % From www.chem.utoronto.ca/coursenotes/analsic/LinRegr2b.pdf
      % sxy=syx/b*sqrt(1/m+1/n+(y0-ybar)^2/( b^2 *sum(xi-xbar)^2 )
      % Modified to handle different measurement noise from model noise
      ctnoise=nanmax(m.ctnoise,std(ct,[],2));
      sx2measure=(ctnoise/m.fit(1)).^2./ctcnt;
      sx2model=(m.ctnoise/m.fit(1)).^2*(1/m.N+(ctmean-mean(m.ct)).^2/(m.fit(1)^2*sum((log(m.conc)-mean(log(m.conc))).^2)));
      sx=sqrt(sx2measure+sx2model);

      % [xx,xxlo,xxhi]=invpred(log(m.conc),m.ct,ctmean,'predopt','observation');
      tfact=tinv(1-(1-obj.options.ci/100)/2,m.N-2);	% T-statistic to provide desired CI
      logconc=(ctmean-m.fit(2))/m.fit(1);
      logconcerr=sx*tfact;

      conc=exp(logconc);
      cilow=exp(logconc-logconcerr);
      cihigh=exp(logconc+logconcerr);

      % Subtract water background concentration from estimated concentrations
      logloderr=m.sxloglod*tfact;
      lod=m.lod;
      lodlow=exp(log(lod)-logloderr);
      lodhigh=exp(log(lod)+logloderr);
      [conc,cilow,cihigh]=obj.concsub(conc,cilow,cihigh,lod,lodlow,lodhigh);
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
        c(4)=30;
        c(3)=5;
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
      hold on;
      if ~isempty(r.wells)
        wellnms=obj.getwellnames(r.wells);
        for i=1:length(r.concs)
          hh=text(r.concs(i)*1.05,r.ct(i),wellnms{i});
          set(hh,'Color','red');
        end
      end
      
      % Draw lines showing best point, confidence intervals
      hold on;

      % Draw fit
      if isfield(r,'mdl')
        ctplot=5:0.1:r.mdl.ctlod;
        [conc,cilow,cihigh]=obj.getconcfromct(ref,ctplot');
        plot(conc,ctplot,'g');
        h(end+1)=plot(cilow,ctplot,'g:');
        leg{end+1}='Reference CI';
        plot(cihigh,ctplot,'g:');

        fitconc=r.mdl.concrange(1):0.1:r.mdl.concrange(2);
        fitct=polyval(r.mdl.fit,log(fitconc));
        h(end+1)=plot(fitconc,fitct,'r');
        leg{end+1}='Linear Fit';
        if isfield(r.mdl,'eff') && isfield(r.mdl,'ct10')
          ti=[ti,sprintf(' eff=%.2f, Conc(Ct=10)=%.3g%s, ctnoise=%.2f, det limit=%.2g%s', r.mdl.eff, r.mdl.ct10,r.units, r.mdl.ctnoise,r.mdl.lod,r.units)];
        end
        c=axis;
        plot(r.mdl.lod*[1,1],c(3:4),':r');
      end
      
      % Overlay data for samples that use this reference
      sampkeys=r.samples.keys;
      allconcs=r.concs;
      for i=1:length(sampkeys)
        s=r.samples(sampkeys{i});
        hh=plot(s.conc,nanmean(s.ct),'xb');
        plot([s.cilow,s.cihigh],nanmean(s.ct)*[1,1],'-b');
        plot(s.conc*ones(1,length(s.ct)),s.ct,'ob');
        plot(s.conc*[1,1],[nanmin(s.ct),nanmax(s.ct)],'-b');
        text(s.cihigh*1.05,nanmean(s.ct),s.name);
        allconcs(end+1)=s.conc;
      end
      if length(sampkeys)>0
        leg{end+1}='Sample';
        h(end+1)=hh;
      end

      c=axis;
      c(1)=max(c(1),min(allconcs(allconcs>0))/10);
      c(2)=min(c(2),max(allconcs(allconcs>0))*10);
      axis(c);

      legend(h,leg);
      title(ti);
      xlabel(sprintf('Equiv. Ref. Concentration (%s)',r.units));
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
    
  end % methods
end % classdef

function y=smooth(x)
  y=medfilt1(x,5);
end

