classdef SigmaZetaMesh < Mesh & helpers.ArraySupport & matlab.mixin.Copyable
    % Defines a sigma-z mesh
    %
    %   The SigmaZetaMesh should be generated with a SigmaZetaMeshGenerator.
    %
    %   The mesh consists of verticals. Each vertical has a certain number of
    %   cells. In each vertical the cells follow the bed. The name Sigma-Zeta
    %   refers to the fact that this mesh is a hyrbid between a z-mesh and a
    %   sigma-mesh. It tries to combine the best of both worlds: it follows
    %   nicely the bed just like a sigma mesh would, but the vertical spacing
    %   of the cells is constant, just like a z-mesh would. This means that
    %   each cell will hold approximately the same number of adcp velocity
    %   estimates.
    %
    %   Each mesh cell consists of six vertices, called:
    %           *         ---> top-middle
    %         /   \
    %       /       \
    %     /          *    ---> top-right
    %   *             |   ---> top-left
    %   |       *     |   ---> bottom-middle
    %   |     /   \   |
    %   |   /       \ |
    %   | /           *   ---> bottom-right
    %   *                 ---> bottom-left
    %
    %   Note that the left, middle and right vertices are always vertically
    %   stacked and share the same n-coordinate.
    %
    %   Data can be stored either in a vector with a value for each cell or
    %   in a matrix which has a toplogy similar to the mesh. The indexing
    %   properties (see below) help map between these formats.
    %
    %   SigmaZetaMesh properties (read only):
    %   * General
    %   xs - defines the cross-section of the mesh
    %   water_level - the water level for the mesh
    %   time - time of the mesh
    %   ncells - number of cells in the mesh
    %   nverticals - number of columns in mesh
    %   max_ncells_vertical - maximum number of cells in one vertical
    %   area_cells - area of cells in the mesh
    %   dn_cells - lateral size of cells
    %
    %   * Vertex Positions
    %   z_bottom_left - z-coordinate of the bottom-left vertex
    %   z_top_left - z-coordinate of the top-left vertex
    %   z_bottom_mid - z-coordinate of the bottom-middle vertex
    %   z_top_mid - z-coordinate of the top-middle vertex
    %   z_bottom_right - z-coordinate of the bottom-right vertex
    %   z_top_right - z-coordinate of the top-right vertex
    %   sig_bottom_left - sigma-coordinate of the bottom-left vertex
    %   sig_top_left - sigma-coordinate of the top-left vertex
    %   sig_bottom_mid - sigma-coordinate of the bottom-middle vertex
    %   sig_top_mid - sigma-coordinate of the top-middle vertex
    %   sig_bottom_right - sigma-coordinate of the bottom-right vertex
    %   sig_top_right - sigma-coordinate of the top-right vertex
    %   n_left - n-coordinate of the left vertices
    %   n_middle - n-coordinate of the middle vertices
    %   n_right - n-coordinate of the right vertices
    %   x_left - x-coordinate of the left vertices
    %   x_middle - x-coordinate of the middle vertices
    %   x_right - x-coordinate of the right vertices
    %   y_left - y-coordinate of the left vertices
    %   y_middle - y-coordinate of the middle vertices
    %   y_right - y-coordinate of the right vertices
    %
    %   * Bed position
    %   zb_left - z-coordinate of bed at left vertices
    %   zb_middle - z-coordinate of bed at middle vertices
    %   zb_right - z-coordinate of bed at right vertices
    %   zb_all - z-coordinate of bed at all vertices from left to right
    %   nb_all - n-coordinate of bed at all vertices from left to right
    %   xb_all - x-coordinate of bed at all vertices from left to right
    %   yb_all - y-coordinate of bed at all vertices from left to right

    %
    %   * Water surface position
    %   nw - n-coordinates of water surface boundaries
    %   xw - x-coordinates of water surface boundaries
    %   yw - y-coordinates of water surface boundaries
    %
    %   * Patch coordinates (for use with patch plotting function)
    %   n_patch - n-coordinate of vertices for use with patch function
    %   x_patch - x-coordinate of vertices for use with patch function
    %   y_patch - y-coordinate of vertices for use with patch function
    %   z_patch - z-coordinate of vertices for use with patch function
    %
    %   * Indexing
    %   col_to_mat - map column based data to matrix layout
    %   row_to_mat - map row based data to matrix layout
    %   mat_to_cell - map matrix to cell vector layout
    %   cell_to_mat - map cell vector layout to matrix layout
    %   row_to_cell - map row based data to cell layout
    %   col_to_cell - map column based data to cell layout
    %
    %   * Mesh topology
    %   neighbors - holds the index of neighboring cells
    %   domains - domain cells belong to
    %   jacobian - jacobian matrix for i,j to sigma,n coordinates
    %
    %   SigmaZetaMesh methods:
    %   index - returns mesh cell indices for given positions
    %   plot - plot the mesh optionally coloring with a given variable
    %   plot3 - 3D plot the mesh optionally coloring with a given variable
    %   plot_neighbors - plot the mesh connectivity and numbering
    %
    %   see also: Mesh, Solver

    properties
        % SigmaZetaMesh/time
        %
        %   scalar datetime object holding the time of the mesh. This is
        %   mainly relevant when the waterlevel is changing in time
        %
        %   see also: SigmaZetaMesh
        time (1,1) datetime
    end
    properties (SetAccess=?SigmaZetaMeshGenerator)
        % SigmaZetaMesh/nverticals
        %
        %   number of verticals in the mesh. In matrix representation this
        %   is the number of columns.
        %
        % see also: SigmaZetaMesh
        nverticals

        % SigmaZetaMesh/max_ncells_vertical
        %
        %   Maximum number of cells in a vertical. This is also the number
        %   of rows in the matrix representation of the mesh
        %
        % see also: SigmaZetaMesh
        max_ncells_vertical

        % SigmaZetaMesh/xs
        %
        %   cross section on which the mesh is generated. Scalar XSection
        %   object
        %
        %   see also:SigmaZetaMesh
        xs (1,1) XSection

        % SigmaZetaMesh/water_level
        %
        %   scalar double value holding the water level for the mesh.
        %
        % see also: SigmaZetaMesh
        water_level (1,1) double

        % SigmaZetaMesh/z_bottom_left
        %
        %   z-coordinates of bottom left vertex of mesh cells.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        z_bottom_left (:,1) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/z_top_left
        %
        %   z-coordinates of top left vertex of mesh cells.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        z_top_left (:,1) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/z_bottom_mid
        %
        %   z-coordinates of bottom center vertex of mesh cells.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        z_bottom_mid (:,1) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/z_top_mid
        %
        %   z-coordinates of top central vertex of mesh cells.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        z_top_mid (:,1) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/z_bottom_right
        %
        %   z-coordinates of bottom right vertex of mesh cells.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        z_bottom_right (:,1) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/z_top_right
        %
        %   z-coordinates of top right vertex of mesh cells.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        z_top_right (:,1) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/n_left
        %
        %   n-coordinates of left vertices of mesh cells.
        %   size: 1 x nverticals
        %
        % see also: SigmaZetaMesh
        n_left (1,:) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/n_middle
        %
        %   n-coordinates of the central vertices of mesh cells.
        %   size: 1 x nverticals
        %
        % see also: SigmaZetaMesh
        n_middle (1,:) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/n_left
        %
        %   n-coordinates of right vertices of mesh cells.
        %   size: 1 x nverticals
        %
        % see also: SigmaZetaMesh
        n_right (1,:) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/zb_left
        %
        %   bed elevation below the right vertices of mesh cells.
        %   size: 1 x nverticals
        %
        % see also: SigmaZetaMesh
        zb_left (1,:) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/zb_middle
        %
        %   bed elevation below the central vertices of mesh cells.
        %   size: 1 x nverticals
        %
        % see also: SigmaZetaMesh
        zb_middle (1,:) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/zb_right
        %
        %   bed elevation below the right vertices of mesh cells.
        %   size: 1 x nverticals
        %
        % see also: SigmaZetaMesh
        zb_right (1,:) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/zb_all
        %
        %   bed elevation below the left, center, and right verticels of
        %   mesh cells. Usefull to plot bed elevation
        %   size: 1 x (1 + 2 * nverticals)
        %
        % see also: SigmaZetaMesh, nb_all
        zb_all (1,:) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/nb_all
        %
        %   n-coordinate below the left, center, and right verticels of
        %   mesh cells. Usefull to plot bed elevation
        %   size: 1 x (1 + 2 * nverticals)
        %
        % see also: SigmaZetaMesh, zb_all
        nb_all (1,:) double {mustBeFinite, mustBeReal}


        % SigmaZetaMesh/nw
        %
        %   n-coordinates of water surface boundaries
        %
        % see also: SigmaZetaMesh, zb_all
        nw (2,:) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/col_to_mat
        %
        % map column based data to matrix layout. Typical use: you want to
        % have matrix with one value for each cell. To get the n-coordinate
        % of the left vertices you can e.g. use:
        %   mesh.n_left(mesh.col_to_mat)
        %
        %   see also: SigmaZetaMesh
        col_to_mat (:,:) double {mustBeInteger, mustBeFinite mustBeReal}

        % SigmaZetaMesh/row_to_mat
        %
        % map row based data to matrix layout.
        %
        %   see also: SigmaZetaMesh
        row_to_mat (:,:) double {mustBeInteger, mustBeFinite mustBeReal}

        % SigmaZetaMesh/mat_to_cell
        %
        % map matrix to cell vector layout. Typical use, you want to fill a
        % matrix representing the mesh, with mesh values, such as the
        % z-coordinate of the top_left vertex:
        %
        %   ZTopLeft = nan(mesh.max_ncells_vertica, mesh.nverticals);
        %   ZTopLeft(mesh.mat_to_cell) = mesh.z_top_left;
        %
        %   see also: SigmaZetaMesh
        mat_to_cell (:,1) logical

        % SigmaZetaMesh/cell_to_mat
        %
        %  map cell vector layout to matrix layout. Typical use: I have a
        %  matrix of size nverticals x max_ncells_vertical:
        %
        %   MyVel = zeros(mesh.nverticals, mesh.max_ncells_vertical);
        %
        %   This can be mapped to a vector with a value for each cell as
        %   follows:
        %
        %   vels_vec = MyVel(mesh.cell_to_mat);
        %
        %   see also: SigmaZetaMesh
        cell_to_mat (:,1) double {mustBeInteger, mustBeFinite mustBeReal}

        % SigmaZetaMesh/row_to_cell
        %
        %   map row based data to cell layout
        %
        %   see also: SigmaZetaMesh
        row_to_cell (:,1) double {mustBeInteger, mustBeFinite mustBeReal}

        % SigmaZetaMesh/col_to_cell
        %
        %   map column based data to cell layout. Typical use: if you want
        %   to e.g. get the n-coordinate of the left vertices for each of
        %   the cells you can use:
        %
        %   n_left_cells = mesh.n_left(obj.col_to_cell)
        %
        %   see also: SigmaZetaMesh
        col_to_cell (:,1) double {mustBeInteger, mustBeFinite mustBeReal}

        % SigmaZetaMesh.neighbors
        %
        %   array with size 4 x ncells holding the index of the
        %   neighbors to the right, top left and bottom in the first to
        %   fourth row respectively. The figure below shows the row number
        %   where the indices of the specific neighbors end up:
        %
        %                     |   row 2 -> top  |
        %           __________|_________________|___________
        %                     |                 |
        %       row 3 -> left |       cell      | row 1 -> right
        %           __________|_________________|___________
        %                     |                 |
        %                     | row 4 -> bottom |
        %
        %   If the neighbor does not exist, e.g. at the corners or sides of
        %   the mesh, a NaN is returned.
        %
        %   see also: SigmaZetaMesh
        neighbors (:,4) double {mustBeReal}

        % SigmaZetaMesh/domains
        %
        % The domain array hold an index for each cell depending on
        % where the cell is located in the mesh:
        %
        %  4: top left corner   |    3: top side    | 2: top right corner
        %  _____________________|___________________|_____________________
        %                       |                   |
        %     5: left side      |    0: interior    | 1: right side
        %  _____________________|___________________|_____________________
        %                       |                   |
        %  6: bottom left corner|  7: bottom side   | 8: bottom right corn
        %
        %   9: degenerate cell, i.e. having outer border on more than
        %   two sides, or having it on the bottom and top or left and
        %   right
        %
        %   see also: SigmaZetaMesh
        domains (:,1) double {mustBeInteger, mustBeFinite, mustBeReal,...
            mustBeNonnegative, mustBeLessThanOrEqual(domains,9)}

        % SigmaZetaMesh/jacobian
        %
        %   jacobian of transformation from i,j to sigma,n coordinates,
        %   with i being the cell number in the vertical (numbering from
        %   surface to bed and j the column number from left to right.
        %
        %   jacobian is an ncells x 2 x 2 with in the trailing dimensions
        %   the jacobian matrix:
        %        __          __
        %       |   dn    dn   |
        %       |   --    --   |
        %       |   di    dj   |
        %       |              |
        %       |  dsig  dsig  |
        %       |  ----  ----  |
        %       |__ di    dj __|
        %
        %   see also: SigmaZetaMesh
        jacobian (:,2,2) double {mustBeReal}
    end
    properties (Dependent, SetAccess=protected, GetAccess=public)
        % SigmaZetaMesh/x_left
        %
        %   x-coordinates of left vertices of mesh cells.
        %   size: 1 x nverticals
        %
        % see also: SigmaZetaMesh
        x_left

        % SigmaZetaMesh/x_middle
        %
        %   x-coordinates of central vertices of mesh cells.
        %   size: 1 x nverticals
        %
        % see also: SigmaZetaMesh
        x_middle

        % SigmaZetaMesh/x_right
        %
        %   x-coordinates of right vertices of mesh cells.
        %   size: 1 x nverticals
        %
        % see also: SigmaZetaMesh
        x_right

        % SigmaZetaMesh/y_left
        %
        %   y-coordinates of left vertices of mesh cells.
        %   size: 1 x nverticals
        %
        % see also: SigmaZetaMesh
        y_left

        % SigmaZetaMesh/y_middle
        %
        %   y-coordinates of central vertices of mesh cells.
        %   size: 1 x nverticals
        %
        % see also: SigmaZetaMesh
        y_middle

        % SigmaZetaMesh/y_right
        %
        %   y-coordinates of right vertices of mesh cells.
        %   size: 1 x nverticals
        %
        % see also: SigmaZetaMesh
        y_right

        % SigmaZetaMesh/z_center
        %
        %   z-coordinates of mesh cells center.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        z_center

        % SigmaZetaMesh/xb_all
        %
        %   x-coordinate of bed points below the left, center, and right
        %   vertices of mesh cells. Usefull to plot bed elevation
        %   size: 1 x (1 + 2 * nverticals)
        %
        % see also: SigmaZetaMesh, nb_all
        xb_all

        % SigmaZetaMesh/yb_all
        %
        %   y-coordinate of bed points below the left, center, and right
        %   vertices of mesh cells. Usefull to plot bed elevation
        %   size: 1 x (1 + 2 * nverticals)
        %
        % see also: SigmaZetaMesh, nb_all
        yb_all

        % SigmaZetaMesh/xw
        %
        %   x-coordinates of water surface boundaries
        %
        % see also: SigmaZetaMesh, nw
        xw

        % SigmaZetaMesh/yw
        %
        %   y-coordinates of water surface boundaries
        %
        % see also: SigmaZetaMesh, nw
        yw

        % SigmaZetaMesh/n_patch
        %
        %   n-coordinates of all vertices of mesh cells to be used for
        %   plotting with the patch function
        %
        %   size: 7 x ncells, with 7 being the coordinates of the 6
        %       vertices with the first being repeated to close the polygon
        %
        % see also: SigmaZetaMesh
        n_patch

        % SigmaZetaMesh/x_patch
        %
        %   x-coordinates of all vertices of mesh cells to be used for
        %   plotting with the patch function
        %
        %   size: 7 x ncells, with 7 being the coordinates of the 6
        %       vertices with the first being repeated to close the polygon
        %
        % see also: SigmaZetaMesh
        x_patch

        % SigmaZetaMesh/y_patch
        %
        %   y-coordinates of all vertices of mesh cells to be used for
        %   plotting with the patch function
        %
        %   size: 7 x ncells, with 7 being the coordinates of the 6
        %       vertices with the first being repeated to close the polygon
        %
        % see also: SigmaZetaMesh
        y_patch

        % SigmaZetaMesh/z_patch
        %
        %   z-coordinates of all vertices of mesh cells to be used for
        %   plotting with the patch function
        %
        %   size: 7 x ncells, with 7 being the coordinates of the 6
        %       vertices with the first being repeated to close the polygon
        %
        % see also: SigmaZetaMesh
        z_patch

        % SigmaZetaMesh/sig_patch
        %
        %   sigma-coordinates of all vertices of mesh cells to be used for
        %   plotting with the patch function
        %
        %   size: 7 x ncells, with 7 being the coordinates of the 6
        %       vertices with the first being repeated to close the polygon
        %
        % see also: SigmaZetaMesh
        sig_patch


        % SigmaZetaMesh/sig_bottom_left
        %
        %   sigma-coordinates of bottom left vertex of mesh cells.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        sig_bottom_left (:,1) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/sig_top_left
        %
        %   sigma-coordinates of top left vertex of mesh cells.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        sig_top_left (:,1) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/sig_bottom_mid
        %
        %   sigma-coordinates of bottom center vertex of mesh cells.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        sig_bottom_mid (:,1) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/sig_top_mid
        %
        %   sigma-coordinates of top central vertex of mesh cells.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        sig_top_mid (:,1) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/sig_bottom_right
        %
        %   sigma-coordinates of bottom right vertex of mesh cells.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        sig_bottom_right (:,1) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/sig_top_right
        %
        %   sigma-coordinates of top right vertex of mesh cells.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        sig_top_right (:,1) double {mustBeFinite, mustBeReal}

        % SigmaZetaMesh/sig_center
        %
        %   sigma-coordinates of mesh cell centers.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        sig_center

        % SigmaZetaMesh/area_cells
        %
        %   area of mesh cells.
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        area_cells

        % SigmaZetaMesh/dn_cells
        %
        %   delta n or lateral size of mesh cells
        %   size: ncells x 1
        %
        % see also: SigmaZetaMesh
        dn_cells
    end
    methods
        function val = get.nverticals(obj)
            val = size(obj.col_to_mat,2);
        end
        function val = get.max_ncells_vertical(obj)
            val = size(obj.col_to_mat,1);
        end
        function val=get.x_left(obj)
            [val,~]=obj.xs.sn2xy(obj.n_left*0, obj.n_left);
        end
        function val=get.x_middle(obj)
            [val,~]=obj.xs.sn2xy(obj.n_middle*0, obj.n_middle);
        end
        function val=get.x_right(obj)
            [val,~]=obj.xs.sn2xy(obj.n_right*0, obj.n_right);
        end
        function val=get.y_left(obj)
            [~,val]=obj.xs.sn2xy(obj.n_left*0, obj.n_left);
        end
        function val=get.y_middle(obj)
            [~,val]=obj.xs.sn2xy(obj.n_middle*0, obj.n_middle);
        end
        function val=get.y_right(obj)
            [~,val]=obj.xs.sn2xy(obj.n_right*0, obj.n_right);
        end
        function val=get.z_center(obj)
            val=(obj.z_bottom_mid+obj.z_top_mid)/2;
        end
        function val=get.xb_all(obj)
            nall=obj.nb_all;
            [val,~]=obj.xs.sn2xy(nall*0, nall);
        end
        function val=get.yb_all(obj)
            nall=obj.nb_all;
            [~,val]=obj.xs.sn2xy(nall*0, nall);
        end
        function val=get.xw(obj)
            nwl=obj.nw;
            [val,~]=obj.xs.sn2xy(nwl*0, nwl);
        end
        function val=get.yw(obj)
            nwl=obj.nw;
            [~,val]=obj.xs.sn2xy(nwl*0, nwl);
        end
        function val=get.x_patch(obj)
            npatch=obj.n_patch;
            [val,~]=obj.xs.sn2xy(npatch*0, npatch);
        end
        function val=get.y_patch(obj)
            npatch=obj.n_patch;
            [~,val]=obj.xs.sn2xy(npatch*0, npatch);
        end
        function val=get.sig_bottom_left(obj)
            val=obj.z_to_sigma(obj.z_bottom_left, reshape(obj.zb_left(obj.col_to_cell),[],1));
        end
        function val=get.sig_top_left(obj)
            val=obj.z_to_sigma(obj.z_top_left, reshape(obj.zb_left(obj.col_to_cell),[],1));
        end
        function val=get.sig_bottom_mid(obj)
            val=obj.z_to_sigma(obj.z_bottom_mid, reshape(obj.zb_middle(obj.col_to_cell),[],1));
        end
        function val=get.sig_top_mid(obj)
            val=obj.z_to_sigma(obj.z_top_mid, reshape(obj.zb_middle(obj.col_to_cell),[],1));
        end
        function val=get.sig_bottom_right(obj)
            val=obj.z_to_sigma(obj.z_bottom_right, reshape(obj.zb_right(obj.col_to_cell),[],1));
        end
        function val=get.sig_top_right(obj)
            val=obj.z_to_sigma(obj.z_top_right, reshape(obj.zb_right(obj.col_to_cell),[],1));
        end
        function val=get.sig_center(obj)
            val=(obj.sig_bottom_mid+obj.sig_top_mid)/2;
        end
        function s=z_to_sigma(obj,z, zb)
            s=(z-zb)./(obj.water_level-zb);
        end
        function z=sigma_to_z(obj,sig, zb, water_level)
            if nargin < 4
                water_level=obj.water_level;
            end
            z=sig.*(water_level-zb)+zb;
        end
        function val=get.n_patch(obj)
            val=[obj.n_left(obj.col_to_cell)
                obj.n_middle(obj.col_to_cell)
                obj.n_right(obj.col_to_cell)
                obj.n_right(obj.col_to_cell)
                obj.n_middle(obj.col_to_cell)
                obj.n_left(obj.col_to_cell)
                obj.n_left(obj.col_to_cell)];
        end
        function val=get.z_patch(obj)
            val=[obj.z_bottom_left,...
                obj.z_bottom_mid,...
                obj.z_bottom_right,...
                obj.z_top_right,...
                obj.z_top_mid,...
                obj.z_top_left,...
                obj.z_bottom_left]';
        end
        function val=get.sig_patch(obj)
            val=[obj.sig_bottom_left,...
                obj.sig_bottom_mid,...
                obj.sig_bottom_right,...
                obj.sig_top_right,...
                obj.sig_top_mid,...
                obj.sig_top_left,...
                obj.sig_bottom_left]';
        end
        function val=get.area_cells(obj)
            val = (1/2.*(obj.n_patch(2,:)-obj.n_patch(1,:)).*(2.*(obj.z_patch(5,:) - obj.z_patch(2,:)) +...
                obj.z_patch(4,:) + obj.z_patch(6,:) - obj.z_patch(3,:) - obj.z_patch(7,:)))';
        end
        function val=get.dn_cells(obj)
            val = (obj.n_patch(3,:) - obj.n_patch(1,:))';
        end
        function mesh=mesh_at_water_level(obj,target_wl, constant_z)
            mesh(numel(target_wl))=SigmaZetaMesh;
            for cm=1:numel(mesh)
                mesh(cm)=obj.copy();
            end
            if nargin<3
                constant_z=false;
            end

            % adapt z_coordinates, unlees z is to be constant

            for cm=1:numel(mesh)
                if ~constant_z
                    mesh(cm).z_bottom_left = obj.sigma_to_z(...
                        obj.sig_bottom_left,...
                        reshape(obj.zb_left(obj.col_to_cell),[],1),...
                        target_wl(cm));
                    mesh(cm).z_top_left=obj.sigma_to_z(...
                        obj.sig_top_left,...
                        reshape(obj.zb_left(obj.col_to_cell),[],1),...
                        target_wl(cm));
                    mesh(cm).z_bottom_mid=obj.sigma_to_z(...
                        obj.sig_bottom_mid,...
                        reshape(obj.zb_middle(obj.col_to_cell),[],1),...
                        target_wl(cm));
                    mesh(cm).z_top_mid=obj.sigma_to_z(...
                        obj.sig_top_mid,...
                        reshape(obj.zb_middle(obj.col_to_cell),[],1),...
                        target_wl(cm));
                    mesh(cm).z_bottom_right=obj.sigma_to_z(...
                        obj.sig_bottom_right,...
                        reshape(obj.zb_right(obj.col_to_cell),[],1),...
                        target_wl(cm));
                    mesh(cm).z_top_right=obj.sigma_to_z(...
                        obj.sig_top_right,...
                        reshape(obj.zb_right(obj.col_to_cell),[],1),...
                        target_wl(cm));

                    % set the target water level in the mesh
                    mesh(cm).water_level=target_wl(cm);

                    % regenerate neighbors
                    mesh(cm).neighbors =...
                        SigmaZetaMeshGenerator.get_neighbors_and_domain(...
                        mesh(cm));

                    % regenerate jacobian
                    mesh(cm).jacobian =...
                        SigmaZetaMeshGenerator.get_jacobian(mesh(cm));
                end
            end

        end

        function idx=index(obj,n,sigma)
            % Indices of mesh cells for given positions
            %
            %   idx = index(obj, n, sigma) returns the indices of the mesh cells that
            %   hold the points given in (n,sigma) coordinates
            %
            %  see also: SigmaZetaMesh
            idx=nan(size(sigma));
            nl=obj.n_left(obj.col_to_cell);
            nm=obj.n_middle(obj.col_to_cell);
            nr=obj.n_right(obj.col_to_cell);
            sbl=obj.sig_bottom_left;
            sbm=obj.sig_bottom_mid;
            sbr=obj.sig_bottom_right;
            stl=obj.sig_top_left;
            stm=obj.sig_top_mid;
            str=obj.sig_top_right;
            % TODO: vectorize code below!
            for cc=1:obj.ncells
                % left side
                fleft = n >= nl(cc) & n < nm(cc);
                fright = n >= nm(cc) & n < nr(cc);
                fsigleft = sigma >  obj.fit_sig(n,nl(cc),nm(cc),sbl(cc),sbm(cc)) &...
                    sigma <= obj.fit_sig(n,nl(cc),nm(cc),stl(cc),stm(cc)) &...
                    fleft;
                fsigright = sigma > obj.fit_sig(n,nm(cc),nr(cc),sbm(cc),sbr(cc)) &...
                    sigma <= obj.fit_sig(n,nm(cc),nr(cc),stm(cc),str(cc)) &...
                    fright;
                fincell = fsigleft | fsigright;
                idx(fincell) = cc;
            end
        end

        function varargout = plot(obj,varargin)
            % Plot the mesh optionally colored with a variable
            %
            %   plot(obj) plots the mesh with the bed and water surface
            %   plot(obj,var) with numel(var) = ncells plots the mesh and colors the cells with the variable var
            %   plot(obj,var) with numel(var) = 2*ncells or 3*ncells plots the mesh and
            %   colors the cells with the variable var and superimposes a quiver plot
            %   consisting of the 2nd and 3rd columns of var in n and z directions.
            %   plot(obj,...,"sigma",...) specify that the plot should be in sigma
            %   coordinates instead of z coordinates. Does not work for
            %   quiver plots yet.
            %
            %   see also: SigmaZetaMesh, plot3
            varargout = cell(1,nargout);
            if ~isscalar(obj)
                tiledlayout("flow");
                [varargout{:}] = obj.run_method('plot', varargin{:});
                return
            end
            inp = inputParser;
            inp.KeepUnmatched = true;
            inp.addOptional('ax', []);
            inp.addOptional('var',[]);
            inp.addParameter('AspectRatio',1,@(x) isscalar(x) && isfinite(x));
            inp.addParameter('FixAspectRatio',false,@(x) isscalar(x) && islogical(x));
            inp.addParameter('Sigma', false, @(x) isscalar(x) && islogical(x));
            inp.parse(varargin{:})
            aspect_ratio = inp.Results.AspectRatio;
            fix_ratio =inp.Results.FixAspectRatio;
            sigma = inp.Results.Sigma;
            if sigma
                fix_ratio = false;
            end
            plot_var=nan(obj.ncells,1);
            get_gca = true;
            for ca = 1:numel(varargin)
                carg = varargin{ca};
                if isa(carg,'double')
                    assert(size(plot_var,1)==obj.ncells,...
                        'SigmaZetaMesh:PlotVarWrongNRows',...
                        'Variable to plot should have same number of rows as cells in the mesh');
                    assert(ismember(size(plot_var,2),[1 3]),...
                        'SigmaZetaMesh:PlotVarWrongNCols',...
                        'Variable to plot should have either one or three columns')
                    plot_var = carg;
                elseif isa(carg,'matlab.graphics.axis.Axes')
                    assert(isscalar(carg),...
                        'SigmaZetaMesh:AxHandleNotScalar',...
                        'Only supports scalar axis handle')
                    ax = carg;
                    get_gca = false;
                elseif (isa(carg,'char') || isa(carg,'string')) &&...
             		    ismember(carg,{'AspectRatio',...
                        'FixAspectRatio', 'Sigma'})
                    varargin(ca)=[]; % skip next argument
                    if ca + 1 > numel(varargin)
                        break
                    end
                    continue
                end
            end
            if get_gca
                ax = nexttile;
            end
            hold_stat=get(ax,'NextPlot');
            if sigma
                zb = obj.nb_all * 0;
                zw = [1 1];
                zpatch = obj.sig_patch;
                zcenter = obj.sig_center;
            else
                zb = obj.zb_all;
                zw = obj.nw*0+obj.water_level;
                zpatch = obj.z_patch;
                zcenter = obj.z_center;
            end
            hbed = plot(ax,obj.nb_all,zb*aspect_ratio,'k','Linewidth',2);
            set(ax,'NextPlot','add')
            hwater = plot(ax,obj.nw,zw*aspect_ratio,'b','Linewidth',2);
            hmesh = patch(ax,obj.n_patch, zpatch*aspect_ratio, plot_var(:,1));
            if nargout>0
                varargout = {hbed, hwater, hmesh};
            end
            if size(plot_var,2)==3
                xl = get(gca,'xlim');
                yl = get(gca,'ylim');
                xpos = xl(1)+diff(xl)/10;
                ypos = yl(1)+diff(yl)/10;
                text(xpos,ypos,'0.2  m/s','HorizontalAlignment','left',...
                    VerticalAlignment='top',BackgroundColor='w')
                q=quiver([obj.n_middle(obj.col_to_cell)'; xpos], ...
                    [obj.z_center*aspect_ratio; ypos], ...
                    [plot_var(:,2); 0.2], [plot_var(:,3); 0], 'k');

                hquiv = quiver(ax, obj.n_middle(obj.col_to_cell)',...
                    zcenter*aspect_ratio, plot_var(:,2), plot_var(:,3),...
                    'Color','k');
                if nargout > 0
                    varargout = [varargout, {hquiv}];
                end
                shading(ax,'flat')
            end
            if fix_ratio
                axis equal
            end

            ylab = cellfun(@str2num, get(gca,'YTickLabel'));
            ylab = ylab/aspect_ratio;
            set(gca,'YTickLabel', ylab)
            set(gca,'YTickLabelMode','manual','YTickMode','manual')

            set(gca,'NextPlot',hold_stat);
        end

        function varargout=plot3(obj,varargin)
            % Plot the mesh optionally colored with a variable in 3D
            %
            %   plot3(obj) plots the mesh with the bed and water surface
            %
            %   plot3(obj,var) plot the mesh and color the cells with the varibale var
            %
            %   see also: SigmaZetaMesh, plot
            varargout = cell(1,nargout);
            if ~isscalar(obj)
                hold_stat = get(gca,'NextPlot');
                hold on
                [varargout{:}] = obj.run_method('plot3', varargin{:});
                set(gca,'NextPlot', hold_stat)
                return
            end
            plot_var=nan(obj.ncells,1);
            get_gca = true;
            for ca = 1:numel(varargin)
                carg = varargin{ca};
                if isa(carg,'double')
                    assert(size(plot_var,1)==obj.ncells,...
                        'SigmaZetaMesh:PlotVarWrongNRows',...
                        'Variable to plot should have same number of rows as cells in the mesh');
                    assert(ismember(size(plot_var,2),[1 3]),...
                        'SigmaZetaMesh:PlotVarWrongNCols',...
                        'Variable to plot should have either one or three columns')
                    plot_var = carg;
                elseif isa(carg,'matlab.graphics.axis.Axes')
                    assert(isscalar(carg),...
                        'SigmaZetaMesh:AxHandleNotScalar',...
                        'Only supports scalar axis handle')
                    ax = carg;
                    get_gca = false;
                end
            end
            if get_gca
                ax = gca;
            end
            hold_stat=get(ax,'NextPlot');
            hbed=plot3(ax,obj.xb_all,obj.yb_all,obj.zb_all, 'k', 'Linewidth',2);
            set(ax,'NextPlot','add')
            hwater=plot3(ax,obj.xw,obj.yw,obj.water_level+obj.xw*0, 'b', 'Linewidth',2);
            hmesh=patch(ax,obj.x_patch,obj.y_patch, obj.z_patch,plot_var(:,1));
            if nargout > 0
                varargout = {hbed, hwater, hmesh};
            end
            if size(plot_var,2)==3
                hquiv = quiver3(ax, obj.x_middle(obj.col_to_cell)',...
                    obj.y_middle(obj.col_to_cell)',...
                    obj.z_center,...
                    plot_var(:,1),...
                    plot_var(:,2), ...
                    plot_var(:,3),...
                    'Color','k');
                varargout = [varargout, {hquiv}];
                shading(ax,'flat')
            end
            set(ax,'NextPlot',hold_stat);
            set(ax,'Clipping', 'off')
            pbaspect([5 5 1])
            da = daspect;
            hrat = max(da(1:2))/da(3);
            daspect([hrat hrat 1])
            view(30,30)
        end
        function varargout=plot_neighbors(obj,varargin)
            varargout = cell(1,nargout);
            if ~isscalar(obj)
                tiledlayout("flow");
                [varargout{:}] = obj.run_method('plot_neighbors', varargin{:});
                return
            end
            ax = nexttile;
            obj.plot(ax, 'sigma', true);
            hold on
            n = obj.n_middle(obj.col_to_cell);
            sig = reshape(obj.sig_center,1,[]);
            plot(n,sig,'b.')
            text(n,sig,cellfun(@num2str,num2cell(1:obj.ncells),'UniformOutput',false))
            nb = obj.neighbors;
            cols = ax.ColorOrder;
            cidx = ax.ColorOrderIndex;
            ncols = size(cols,1);
            lwidth = [1 1 1 1];
            for cq = 1:4
                isf = isfinite(nb(:,cq));
                dn = n(nb(isf,cq)) - n(isf);
                dsig = sig(nb(isf,cq)) - sig(isf);
                plot([n(isf); n(isf) + dn/3],...
                    [sig(isf); sig(isf) + dsig/3],...
                    'color', cols(mod(cidx+cq, ncols)+1,:),...
                    'linewidth', lwidth(cq))
            end
        end
    end
    methods(Access=protected)
        function val=get_ncells(obj)
            val=size(obj.z_bottom_left,1);
        end
    end
    methods(Access=protected, Static)
        function sig=fit_sig(n,n0,n1,sig0,sig1)
            sig=(sig1-sig0)./(n1-n0).*(n-n0)+sig0;
        end
    end

end
