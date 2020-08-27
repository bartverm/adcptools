function csid=cross_section_selector(vmadcp)
% Interactively define cross sections for vessel mounted adcp data
%
%   csid = cross_section_selector(vmadcp) will allow to select
%   cross-sections to process adcp data by drawing polygons around the
%   track. vmadcp must be a VMADCP object.
%
%   see also: VMADCP
    vmadcp.plot_track;
    title('Press Enter to end selection')
    [xa,ya]=vmadcp.xy;
    hold on
    uiwait(msgbox('Click on the figure to draw a polygon around the repeat transect track and press Enter to end'))
    get_next=true;
    csid=[];
    inpol=false(size(xa));
    polygon=patch();
    polygon.FaceColor='none';
    selection=plot(1,1,'ro');
    leg=legend(selection,'Current selection');
    while (get_next)
        polygon.XData=[];
        polygon.YData=[];
        selection.XData=[];
        selection.YData=[];
        while (true)
            try
                [xt,yt]=ginput(1);
            catch err
                if strcmp(err.identifier,'MATLAB:ginput:FigureDeletionPause')
                    uiwait(msgbox('You closed the figure, selection ends here!','Figure closed'))
                    return
                else
                    rethrow(err)
                end
            end
            if isempty(xt)
                break
            end
            set(polygon,'XData',[get(polygon,'XData'); xt],'YData',[get(polygon,'YData'); yt])
            if numel(get(polygon,'XData'))>2
                inpol=inpolygon(xa,ya,get(polygon,'XData'),get(polygon,'YData'));
                set(selection,'XData',xa(inpol),'YData',ya(inpol));
            end
        end
        if ~any(inpol)
            uiwait(msgbox('The polygon did not contain any adcp positions!','Empty polygon!'))
        else
            if strcmp('Yes',questdlg('Are you happy with this selction','Selection finished!','Yes','No','Yes'))
                csid=[csid; zeros(1,numel(inpol))]; %#ok<AGROW>
                csid(end,inpol)=1;
            end
        end
        if strcmp('No',questdlg('Do you want to select another repeat transect?','Another one','Yes','No','Yes'))
            break
        else
            if any(inpol)
                plot(xa(csid(end,:)==1),ya(csid(end,:)==1),'o')
                leg.String{end}=['Cross section ',num2str(size(csid,1))];
            end
        end
    end
    close(gcf)
end