function [ef, xs]=cross_section_selector(vmadcp)
% Interactively define cross sections for vessel mounted adcp data
%
%   [ef, xs] = cross_section_selector(vmadcp) will allow to select
%   cross-sections to process adcp data by drawing polygons around the
%   track. vmadcp must be a VMADCP object. It will return a vector with
%   EnsembleFilter objects and XSection objects defining the different
%   cross-sections
%
%   see also: VMADCP, EnsembleFilter, XSection
    figure
    vmadcp.plot_track('-','color',[.5 .5 .5]);
    title('Press Enter to end selection')
    hpos=[vmadcp.horizontal_position];
    [xa,ya]=deal(hpos(1,:),hpos(2,:));
    hold on
    uiwait(msgbox('Click on the figure to draw a polygon around the repeat transect track and press Enter to end'))
    get_next=true;
    inpol=false(size(xa));
    polygon=patch();
    polygon.FaceColor='none';
    tang_quiv=quiver(mean(xa,'omitnan'),mean(ya,'omitnan'),0,0,'color','r');
    norm_quiv=quiver(mean(xa,'omitnan'),mean(ya,'omitnan'),0,0,'color','g');
    selection=plot(1,1,'r.');
    leg=legend([selection tang_quiv, norm_quiv],'Current selection', 'Tangential dir','Normal dir');
    leg.AutoUpdate='off';
    ef=EnsembleFilter.empty(1,0);
    xs=XSection.empty(1,0);
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
                    % when figure is close just end selecting and return
                    % result so far
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
            if strcmp('Yes',questdlg('Are you happy with this selction?','Selection finished!','Yes','No','Yes'))
                ef=[ef EnsembleFilter(~inpol)]; %#ok<AGROW>
                xs=[xs XSection(vmadcp,ef(end))]; %#ok<AGROW>
                leg.AutoUpdate='on';
                plot(xa(~ef(end).bad_ensembles),ya(~ef(end).bad_ensembles),'.')
                leg.String{end}=['Cross section ',num2str(numel(ef))];
                leg.AutoUpdate='off';
                quivs=xs(end).plot();
                if strcmp('Yes', questdlg('Should the direction of the cross-section be reversed?', 'Cross section generated','No','Yes','No'))
                    xs(end).revert();
                    delete(quivs)
                    xs(end).plot();
                end
            end
        end
        if strcmp('No',questdlg('Do you want to select another repeat transect?','Another one','Yes','No','Yes'))
            break
        end
    end
    leg.AutoUpdate='on';
    title('Cross sections')
    delete(selection);
    delete(polygon);
end