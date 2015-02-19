% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

classdef plotobj < imobj  
    properties
        ploth = [];
        trmat = [];
        plotcolor =[];
        plotargs = {};
        %show = false;
    end
    
    properties (Dependent = true)
        linestyle;
        marker = 'none';
        show;
    end
    
    properties (Hidden = true)
        lnstyle = 'none';
        mrkstyle = '+';
        visible = false;
    end
    
    methods
        
        function set.linestyle(me,a)
            
            ism = ismember(a,'-:');
            lnst = a(ism);
            mrkst = a(~ism);
            
            if isempty(lnst)
                lnst = 'none';
            end
            if isempty(mrkst)
                mrkst = 'none';
            end
            me.lnstyle = lnst;
            me.marker = mrkst;
        end
        
        function a = get.linestyle(me)
            a = regexprep({me.lnstyle,me.mrkstyle},'none','');
            a = [a{:}];
        end
        
        function set.marker(me,a)
           
              if isempty(a) || ~ismember(a,{'+','o','*','.','x','square','diamond','v','^','>','<','pentagram','hexagram','none'})
                   a = 'none';
              end
              me.mrkstyle = a;
        end
        function a = get.marker(me)
           a = me.mrkstyle; 
        end
        
        function set.show(me,a)
           
            if a == me.visible
                return
            end
            me.visible = a;
            ish = ishandle(me.ploth);
            if any(ish) 
                
                if a
                    set(me.ploth(ish),'visible','on')
                else
                    set(me.ploth(ish),'visible','off')                    
                end
            end
            
        end
            
        function a = get.show(me)
           
            a = me.visible;
                        
        end
        
    end
end


