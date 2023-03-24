classdef ModelParameters < handle
    %MODELPARAMETERS Class for capturing model parameters as output of
    %get_parameters

    %   Detailed explanation goes here

    properties
        M (:,:) double % matrix M such that b = Mp

        b (:,1) double; %rhs of system of eqs (b = Mp)

        p (:,:) double; % model parameters (b = Mp)

        reg (1,1) Regularization

        opts (1,1) SolverOptions
    end
    properties(Dependent)
        vel_cmap
        phi_cmap
        s_cmap
    end

    methods
        function obj = ModelParameters(varargin)
            % Overwrite default options
            for ia = 1:2:nargin
                obj.(varargin{ia}) = varargin{ia+1};
            end
        end

        function plot_solution(obj, names_selection, par_idx)
            
            if nargin < 2
                names_selection = [obj.reg.model.names{:}];
            end

            P = obj.p(:, par_idx);
            nreg = size(P, 2);
            nn = length(names_selection);
            nc = obj.reg.mesh.ncells;
            np = sum(obj.reg.model.npars); % Number of parameters in each cell
            Np = size(P,1); %= nc*np;

            if nreg>1 % Compare different vectors
                t = tiledlayout(nn, nreg, TileSpacing = "tight", Padding = "tight", TileIndexing = "columnmajor");
            else
                t = tiledlayout('flow', TileSpacing="tight", Padding="tight");
            end

            t.XLabel.String = 'y [m]';
            t.YLabel.String = 'z [m]';
            t.XLabel.Interpreter = 'latex';
            t.YLabel.Interpreter = 'latex';

            par_idx = obj.get_par_idx(names_selection);
            titles = obj.modify_names(names_selection);

            for col = 1:nreg
                for row = 1:nn
                    nexttile;
                    var = P(par_idx(row):np:Np, col);
                    obj.reg.mesh.plot(var)
                    %     loc_tit = str,old,new)
                    title(['$', titles{row}, '$'], 'interpreter', 'latex')

                    %         set(c,'TickLabelInterpreter','latex')
                    if ~contains(titles{row}, 'phi')
                        amax = max(abs(var), [], 'omitnan') + 1e-5;
                        caxis([-amax, amax])
                        c = colorbar;
                        colormap(gca, obj.vel_cmap)
                        if ~contains(titles{row}, '\partial')
                            ylabel(c, '$m/s$','Rotation',270, 'interpreter', 'latex');
                        elseif contains(titles{row}, '\sigma')
                            ylabel(c, '$m/s$','Rotation',270, 'interpreter', 'latex');
                        else
                            ylabel(c, '$m/s^2$','Rotation',270, 'interpreter', 'latex');
                        end

                    else
                        caxis([-180, 180])
                        temp = get(gca, 'Children');
                        temp(2).CData =  temp(2).CData*180/pi;
                        c = colorbar;
                        colormap(gca, obj.phi_cmap)
                        ylabel(c, 'deg','Rotation',270, 'interpreter', 'latex');

                    end
                    pos = get(c,'Position');
                    if row == 1
                        pos1 = pos;
                    end
                    %         disp(pos)
                    c.Label.Position(1) = pos1(1)+.5/col; % to change its position
                    %         c.Label.Position(2) = c.Label.Position(2) + .2; % to change its position
                    c.Label.HorizontalAlignment = 'center'; % to change its position
                    c.TickLabelInterpreter = 'latex';
                    %         c.Label.Rotation = 270; % to rotate the text
                    axis tight
                    %     hAxes.TickLabelInterpreter = 'latex';
                    %     title(sprintf('%s, %s, %s', names{i}, names{i+np(1)}, names{i + np(1) + np(2)}))
                    set(gca, 'XDir','reverse') % Very important
                    %         xlabel('y [m]', 'interpreter', 'latex')
                    %         ylabel('z [m]', 'interpreter', 'latex')
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                end
            end
            % Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
            % %[Left Bottom Right Top] spacing
            % NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
            % set(gca, 'Position', NewPos);
        end

        function par_idx = get_par_idx(obj, names_selection)
            par_idx = nan([1,numel(names_selection)]);
            for name_idx = 1:numel(names_selection)
                par_idx(name_idx) = find(strcmp([obj.reg.model.names{:}], names_selection{name_idx}));
            end
        end

        function mod_names = modify_names(obj, names_selection)
            mod_names = names_selection;
            for idx = 1:numel(names_selection)
                mod_names{idx} = strrep(mod_names{idx}, 'sig', '\sigma');
                mod_names{idx} = strrep(mod_names{idx}, '^1', '');
                mod_names{idx} = strrep(mod_names{idx}, 'u0', 'u_0');
                mod_names{idx} = strrep(mod_names{idx}, 'v0', 'v_0');
                mod_names{idx} = strrep(mod_names{idx}, 'w0', 'w_0');
                mod_names{idx} = strrep(mod_names{idx}, 'd', '\partial ');
            end
        end


    
        function vmap = get.vel_cmap(obj)
            
            vmap = brewermap(20, 'RdBu');
        end

        function smap = get.s_cmap(obj)
            smap = brewermap(15, 'YlOrBr');
        end

        function pmap = get.phi_cmap(obj)
            pmap = [obj.vel_cmap ;flipud(obj.vel_cmap)];
        end

    end  
end

