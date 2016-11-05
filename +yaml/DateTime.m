classdef DateTime
% 		Copyright (c) 2011
% 		This program is a result of a joined cooperation of Energocentrum
% 		PLUS, s.r.o. and Czech Technical University (CTU) in Prague.
%         The program is maintained by Energocentrum PLUS, s.r.o. and 
%         licensed under the terms of MIT license. Full text of the license 
%         is included in the program release.  		
%         Author(s): 
% 		Jiri Cigler, Dept. of Control Engineering, CTU Prague & Automatic Control Laboratory, ETH Zurich		
% 		Jan  Siroky, Energocentrum PLUS s.r.o.		
%         Implementation and Revisions:
%         Auth  Date        Description of change
%         ----  ---------   -------------------------------------------------
%         jc    01-Mar-11   First implementation
%         jc    30-Sep-11   Added function colon
%         jc    07-Jan-12   Added functions addtodate,datevec,weekday

    properties
        serialDate
    end
    methods

function this = DateTime(varargin)    
    import yaml.*;
    if numel(varargin)==1 && isa(varargin{1},'java.util.Date')
                    sec = varargin{1}.getTime/1000;                   
                    this.serialDate=datenum(1970,1,1,0,0,sec);
            else
                this.serialDate=datenum(varargin{:});
            end
        end

function this = plus(this,val)      
    import yaml.*;
    o =@plus;
            this = doFun(this,o,val);
        end

function this = minus(this,val)    
    import yaml.*;
    o =@minus;
            this = doFun(this,o,val);
        end

function this = times(this,val)    
    import yaml.*;o =@times;
            this = doFun(this,o,val);
        end

function this = mtimes(this,val)        
    import yaml.*;o =@mtimes;
            this = doFun(this,o,val);
        end

function this = mrdivide(this,val)        
    import yaml.*;o =@mrdivide;
            this = doFun(this,o,val);
        end

function this = rdivide(this,val)      
    import yaml.*;o =@rdivide;
            this = doFun(this,o,val);
        end

function this = horzcat(this,varargin)   
    import yaml.*;for i=1:numel(varargin)
                this.serialDate = [this.serialDate, varargin{i}.serialDate];
            end
        end

function out = colon(this,step,to)       
    import yaml.*;vect = [double(this):double(step):double(to)]';
            out =DateTime(vect);
        end

function this = vertcat(this,varargin)     
    import yaml.*;for i=1:numel(varargin)
                this.serialDate = [this.serialDate; varargin{i}.serialDate];
            end
        end

function this = ctranspose(this)       
    import yaml.*;this.serialDate = this.serialDate';
        end

function this = transpose(this)     
    import yaml.*;this.serialDate = this.serialDate';
        end

function  disp(this)      
    import yaml.*;disp([this.serialDate])
        end

function out = double(this)  
    import yaml.*;out = this.serialDate;
        end

function out = length(this)    
    import yaml.*;out = length(this.serialDate);
        end

function out = size(this,varargin)      
    import yaml.*;out = size(this.serialDate,varargin{:});
        end

function out = numel(this)   
    import yaml.*;out = numel(this.serialDate);
        end

function out = isreal(this)       
    import yaml.*;out = isreal(this.serialDate);
        end

function out = isnan(this) 
    import yaml.*;out = isnan(this.serialDate);
        end

function out = isfinite(this)     
    import yaml.*;out = isfinite(this.serialDate);
        end

function out = le(this,B)    
    import yaml.*;if isa(B,'DateTime')
                out = le(this.serialDate,B.serialDate);
            else
                out = le(this.serialDate,B);
            end
        end

function out = lt(this,B)  
    import yaml.*;fun=@lt;
            if isa(B,'DateTime')
                out = fun(this.serialDate,B.serialDate);
            else
                out = fun(this.serialDate,B);
            end
        end

function out = gt(this,B)    
    import yaml.*;fun=@gt;
            if isa(B,'DateTime')
                out = fun(this.serialDate,B.serialDate);
            else
                out = fun(this.serialDate,B);
            end
        end

function out = eq(this,B)       
    import yaml.*;fun=@eq;
            if isa(B,'DateTime')
                out = fun(this.serialDate,B.serialDate);
            else
                out = fun(this.serialDate,B);
            end
        end

function out = diff(this)      
    import yaml.*;out = diff(this.serialDate);
        end

function out = norm(this,varargin)    
    import yaml.*;out = norm(this.serialDate,varargin{:});
        end

function [this k] = sort(this,varargin)      
    import yaml.*;[this.serialDate k] = sort(this.serialDate,varargin{:});
        end

function this = subsref(this,S)     
    import yaml.*;if isa(S.subs{1},'DateTime')
                S.subs{1}=double(S.subs{1});
            end
            this.serialDate =  subsref(this.serialDate,S);
        end

function idx = subsindex(this)     
    import yaml.*;idx = double(this)-1;
        end

function endidx = end(this,k,n)      
    import yaml.*;if size(this.serialDate,1)==1 || size(this.serialDate,2)==1
                endidx=numel(this.serialDate);
            else
                endidx = size(this.serialDate,k);
            end
        end

function this = subsasgn(this, S, B)     
    import yaml.*;if not(isa(B,'DateTime'))
                B=DateTime(B);
            end
            this.serialDate =subsasgn(this.serialDate, S, B);
        end

function res = bsxfun(fun,A,B)       
    import yaml.*;res = fun(A,B);
        end

function out =superiorfloat (x,y,xi)    
    import yaml.*;if isa(x,'DateTime') && isa(xi,'DateTime')
                out = superiorfloat(x.serialDate,y,xi.serialDate);
            elseif isa(x,'DateTime') && not(isa(xi,'DateTime'))
                out = superiorfloat(x.serialDate,y,xi);
            elseif not(isa(x,'DateTime')) && isa(xi,'DateTime')
                out = superiorfloat(x,y,xi.serialDate);
            else
                out = superiorfloat(x,y,xi);
            end
        end

function this = floor(this)      
    import yaml.*;this.serialDate = floor(this.serialDate);
        end

function this = max(this,varargin)   
    import yaml.*;this.serialDate = max(this.serialDate,varargin{:});
        end

function this = min(this,varargin)    
    import yaml.*;this.serialDate = min(this.serialDate,varargin{:});
        end

function out = datestr(this,varargin)  
    import yaml.*;out = datestr(this.serialDate,varargin{:});
        end

function out = addtodate(this,varargin)     
    import yaml.*;out = addtodate(this.serialDate,varargin{:});
        end

function varargout= datevec(this,varargin)      
    import yaml.*;nout = nargout;
            if nout <=1
                varargout{1} = datevec(this.serialDate,varargin{:});
            elseif nout ==2
                [varargout{1} varargout{2}] = datevec(this.serialDate,varargin{:});
            elseif nout ==3
                [varargout{1} varargout{2} varargout{3}] = datevec(this.serialDate,varargin{:});
            elseif nout ==4
                [varargout{1} varargout{2} varargout{3} varargout{4}] = datevec(this.serialDate,varargin{:});
            elseif nout ==5
                [varargout{1} varargout{2} varargout{3} varargout{4} varargout{5} ] = datevec(this.serialDate,varargin{:});
            elseif nout ==6
                [varargout{1} varargout{2} varargout{3} varargout{4} varargout{5} varargout{6} ] = datevec(this.serialDate,varargin{:});
            else 
                error('Unknown function call');
            end
        end
    end
    methods (Access = private)

function this = doFun (this,o, val)
end
    end
end