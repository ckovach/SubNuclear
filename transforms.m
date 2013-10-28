% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

classdef transforms   
    %Class definition for a volumeview transform class.
    properties
        label = ''; 
        notes= '';
        file = '';
        units = 'units';
        type = 'linear'; %Transformation type
        objectid = -1;
%         tr; %Transform function
%         itr; %Inverse transform function
        nldfun = @(x)0; %Non-linear diplacement
        inldfun = @(x)0; %Non-linear displacement for the inverse 
       
        
    end
    properties (SetAccess=private)
        mytr = [];
       chain = [];
    end
    properties (Dependent=true)
        trmat ;  % Linear transformation
%         DispNonlin; % Non-linear displacement function (for non-linear transformations)
    end
        
   methods
        function me = transforms(varargin)
            
             for i = 1:2:length(varargin)
                me.(varargin{i})=varargin{i+1};
            end
        end
        
        function  a = invtrmat(me)
            if ~isempty(me.trmat)
                a = me.trmat^-1;
            else
                a = [];
            end
        end
        
        function me = set.trmat(me,T)
            if (T(end)~=1 || sum(abs(T(:,end)))~=1) % Check whether matrix include translation
                T(end,end+1) = 1;
            end
            me.mytr = double(T);
            me.chain = me;
%             me.tr = @(pt)me.applytr(pt);
%             me.itr = @(pt)me.applyitr(pt);
       
        end
       function T = get.trmat(me)
             T=me.mytr ;
       end
%        %%%
%         function dfun = get.DispNonlin(me)
%              dfun=me.nldfun;
%         end
        %%%
%         function me = set.DispNonlin(me,a)
%              if ~isa(a,'function_handle')
%                  error('Input must be a function handle with one input argument')          
%              end
%             me.nldfun=a;
%         end
     
        function b = mtimes(me,a)
            b = transforms('label',sprintf('%s*%s',me.label,a.label),'trmat',me.trmat*a.trmat,'chain',[me.chain,a.chain]);
        end
        function b = eq(me,a)
           if length(me)==1
               x1 = me.trmat;
               c = a;               
           elseif length(a)==1
              x1 = a.trmat;
               c = me;
           else
               b = arrayfun(@(d,e) numel(d.trmat)==numel(e.trmat) && all(d.trmat(:)==e.trmat(:)),me,a);
               return
           end
           b = arrayfun(@(d) numel(d.trmat)==numel(x1) && all(d.trmat(:)==x1(:)),c);

        end
        %%%
        function me = compute(me,fromx,tox)
            
            nd = size(fromx,2);
            switch nd
                case 3
                    fromx(:,4) = 1;
                case 2
                     tox(:,3) = 1;
            end
            me.trmat = tox\fromx;
            
        end
       %%%
       function a = tr(me,pt)
             
                % apply transform
                nd = size(pt,2);

                if nd < size(me.trmat,1)
                    pt(:,end+1) = 1;
                end
                a = pt*me.mytr(:,1:end-1) + me.nldfun(pt(:,1:end-1));
                
%                 a = a(:,1:end-1);
       end
       function a = itr(me,pt)
         % apply inverse  transform 
         
                nd = size(pt,2);
                if nd < size(me.trmat,1)
                    pt(:,end+1) = 1;
                end
                T = me.mytr^-1;
                a = pt*T(:,1:end-1) + me.inldfun(pt(:,1:end-1));
%                 a = a(:,1:end-1);
       end
   end
end