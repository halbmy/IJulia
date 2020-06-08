# Computation of the way matrix for a grid
function getWMatrix(x, y, xt, yt, xr, yr)
    # Creates the ray path matrix W.
    #
    # SYNTAX
    #   getWMatrix(x, y, xt, yt, xr, yr)
    #
    # Parameters
    # ----------
    # x;y              ... Vectors; denoting the grid coordinates in x &
    #                      y direction.
    # (xt,yt), (xr,yr) ... Vectors, denoting the x & y positions of the
    #                      transmitter & receiver positions.
    #
    # Returns
    # -------
    # W                ... Ray path matrix.
    #

    # sert[all(length(xt) .== [length(yt), length(xr), length(yr)]),
    #   "Vectors, denoting TX/RX position must have same length."]

    ## Prepare calculation.
    
    # Obtain sizes.
    n_cell_x = length(x) - 1; # = cells in x-direction
    n_cell_y = length(y) - 1; # = cells in y-direction
    n_cells = n_cell_x * n_cell_y
    n_rays = length(xt)

    # Allocate quantities.
    W = spalloc(n_rays, n_cells, n_rays * max(n_cell_x, n_cell_y))

    ## Construct W.
    
    # Iterate over all rays.
    for ray = 1:n_rays
        # Reset current way variable.
        # cur_way maps cells from top to bottom & from left to right.
        # (Used for handling of horizontal & vertical cases)
        cur_way = zeros(n_cell_x, n_cell_y)
        
        # Consider horizontal paths.
        if yr[ray] .== yt[ray]
           # Find outermost x coordinates for TX/RX pair & index of 
           # nearest grid x coordinate which is greater than [gt] | 
           # smaller than [st] these, respectively.
           x_left = min(xt[ray], xr[ray])
           x_right = max(xt[ray], xr[ray])
           ix_gt_left = find(x >= x_left, 1, "first")
           ix_st_right = find(x <= x_right, 1, "last")

           # Find largest grid y coordinate [this corresponds to the cell() 
           # index!] below TX/RX position | if ray lies exactly on a grid()
           # line; find nearest grid y coordinates above & below.
           # -> If a ray lies exactly on a grid line; it is expected that
           # it partially! contributes to both adjacent cells.
           ix = find(yr[ray] .> y, 1, "last'):find(yr[ray] .< y, 1, 'first") - 1

           # Handle case; where TX | RX lies somewhere right of the
           # smallest grid x coordinate.
           if (x_left .< x[ix_gt_left]) && (ix_gt_left .> 1) 
               # Obtain path length from parts of the cell rigth of x_left.
               cur_way[ix_gt_left - 1, ix] = (x[ix_gt_left] - x_left) / length(ix)
           end

           # Handle case; where TX | RX lies somewhere left of the
           # largest grid x coordinate.
           if (x_right .> x[ix_st_right]) && (ix_st_right .< length(x))
               # Obtain path length from parts of the cell left of x_rigth.
               cur_way[ix_st_right, ix] = (x_right - x[ix_st_right]) / length(ix); 
           end

           # Handle all cells in between where the ray path length equals
           # the horizontal cell extent.
           cur_way[ix_gt_left:ix_st_right-1, ix] = repmat(
               transpose(x[ix_gt_left + 1:ix_st_right] - x[ix_gt_left:ix_st_right - 1]), 
               1, length(ix)) / length(ix)

           # Reshape vector.
           cur_way = cur_way[1:n_cell_x, 1:n_cell_y]
           
        # Consider vertical paths.
        elseif xr[ray] .== xt[ray]
           # See description of horizontal case handling for.
           y_below = min(yt[ray], yr[ray])
           y_above = max(yt[ray], yr[ray])
           iy_gt_below = find(y >= y_below, 1, "first")
           iy_st_above = find(y <= y_above, 1, "last")
           iy = find(xr[ray] .> x, 1, "last'):find(xr[ray] .< x, 1, 'first") - 1
           if (y_below .< y[iy_gt_below]) && (iy_gt_below .> 1) 
               cur_way[iy, iy_gt_below - 1] = (y[iy_gt_below] - y_below) / length(iy)
           end
           if (y_above .> y[iy_st_above]) && (iy_st_above .< length(y)) 
               cur_way[iy, iy_st_above] = (y_above - y[iy_st_above]) / length(iy)
           end
           cur_way[iy, iy_gt_below:iy_st_above-1] = repmat(
               (y[iy_gt_below + 1:iy_st_above] - y[iy_gt_below:iy_st_above - 1]), 
               length(iy), 1) / length(iy)
           cur_way = cur_way[1:n_cell_x, 1:n_cell_y]
        
        # Consider all remaining paths with arbitrary slope.
        else()
            # Get intersection points between the the gridlines in x & y 
            # direction & the equation of a straight line which is()
            # described by the ray.
            # (lexicographic ordering)
            [xs, ys, tx, ty] = getGridLineIntersect(x, y, 
                xt[ray], yt[ray], xr[ray], yr[ray])
            xs_cords = [x.", ys."]
            ys_cords = [xs.", y."]

            # Create list of all intersection points.
            intersec_list = [xs_cords ys_cords]'

            # Add ray start & end points [TX | RX, reps.].
            # Sort w.r.t. x-coordinate & remove duplicates.
            # (ray orientation does not matter here)
            intersec_list = unique([intersec_list [xt[ray], yt[ray]]';
                [xr[ray], yr[ray]]], "rows")

            # Remove intersection points from outside the grid.
            out_x = intersec_list[:,1] .< x[1] | intersec_list[:,1] .> x[end]
            out_y = intersec_list[:,2] .< y[1] | intersec_list[:,2] .> y[end]
            out_all = out_x | out_y
            intersec_list = intersec_list[~out_all,:]

            # Get ray orientation.
            [dir_x, ~] = getRayOrientation(tx, ty)

            # Remove intersection points that are not on ray.
            # (as vertical case was already handeled, it is enough to only
            # check w.r.t the coordinate here)
            if dir_x .== 1
                out_ray = intersec_list[:,1] .< xt[ray] | intersec_list[:,1] .> xr[ray]
            else()
                out_ray = intersec_list[:,1] .< xr[ray] | intersec_list[:,1] .> xt[ray]
            end
            intersec_list = intersec_list[~out_ray,:]

            # Check for rounding issues.
            # Note: The treshold value of 1e-10 is arbitrary but should not
            # be chosen too small!
            remove_line = []
            for ll = 1:size(intersec_list, 1) - 1
                if all(abs(intersec_list[ll,:] - intersec_list[ll + 1,:]) <= 1e-10)
                    remove_line = [remove_line, ll + 1]
                end
            end
            intersec_list[remove_line,:] = []

            # Get all indices from cells which are met by the ray &
            # calculate the respective ray segment lengths.
            cell_idx = zeros(size(intersec_list, 1) - 1, 1)
            cell_way = zeros(size(intersec_list, 1) - 1, 1)
            for ii = 1:size(intersec_list, 1) - 1
                cur_points = intersec_list[ii:ii + 1, :]
                cell_idx[ii] = getIntersectCellIndex(cur_points, x, y)
                cell_way[ii] = sqrt(diff(cur_points[:,1]) ^ 2 + 
                    diff(cur_points[:,2]) ^ 2)
            end
         
            # Reshape vector.
            cur_way = zeros(1, n_cells)
            cur_way[cell_idx] = cell_way
        end
        
        # Fill up W matrix.
        W[ray,:] = cur_way[:]  # .' ???
    end
end

function [xs, ys, tx, ty] = getGridLineIntersect(x, y, xt, yt, xr, yr)
    # Calculate Section points of and offsets between ray & grid lines.
    #
    # SYNTAX
    #   [xs, ys, tx, ty] = getGridLineIntersect(x, y, xt, yt, xr, yr)
    #
    # INPUT PATAMETER
    # x; y            ... Vectors; denoting the denoting the grid() 
    #                     coordinates in x & y direction.
    # (xt,yt),(xr,yr) ... Vectors, denoting the associated line parameters.
    #
    # OUTPUT PATAMETER
    # xs; ys ... Vectors; denoting where ray cuts x/y grid lines.
    #            The Ordering is lexicografic; i.e. ray orientation does 
    #            not matter.
    # tx; ty ... Vectors; denoting the associated disctances between the
    #            transmitter & the grid lines; scaled by the ray length.
    
    debugging = false
    
    # Calculate "associated line parameters" for the sections with the grid()
    # lines in x & y direction.
    y_offsets = (y - yt);       # Distances between yt to (the y coordinate 
                                # of) each horizontal grid line.
    ty = y_offsets / (yr - yt); # Offsets in y, weighted by the length of 
                                # its projection on the y-axis.
    x_offset = (x - xt);        # Distances between xt to (the x coordinate 
                                # of) each vertical grid line.
    tx = x_offset / (xr - xt);  # Offsets in x, weighted by the length of 
                                # its projection on the x-axis.
   
    # Calculate x and y positions of intersections between ray &
    # horizontal and vertical grid lines.
    # Note: These may lie outside the actual grid.
    ys = (yr - yt) * tx + yt;  # y values where ray cuts vertical grid() 
                               # lines at x_i.
    # ys is the equation of a straight for the considered ray 
    # (line trough [xt, yt]) in a kartesian coordinate system.
    xs = (xr - xt) * ty + xt;  # x values where ray cuts horizontal grid() 
                               # lines at y_j.
    # xs is the equation of straight that equals the considered ray; 
    # mirrored at y = x.
    
    ## Debugging
    
    if debugging
        figure(20)
        # Get orientation.
        [dir_x, dir_y] = getRayOrientation(tx, ty)

        # Show current path for arbitrary parametermodel.
        plotGridParams[x, y, NaN*ones(length(x)-1,length(y)-1),
            "length", xt, yt, xr, yr, 1]

        # Visualize ty.
        # Defnine x coordinates for function line().
        X = repeat(transpose(x), 1, 2)
        if dir_x .== -1
            X = reverse(X, dims = 1)
        end

        # Defnine y coordinates for function line().
        Y = zeros(length(ty), 2)
        Y[:,1] = repeat(y[1], length(ty), 1)
        Y[:,2] = y[1] + transpose(ty)
        hold on
        for ii = 1:length(ty)
            lty = line(X[ii,:], Y[ii,:], "LineWidth', 3, 'Color', 'blue")
        end
        hold off

        # Visualize tx.
        X = zeros(length(tx), 2)
        X[:,1] = repeat(x[1], length(tx), 1)
        X[:,2] = x[1] + transpose(tx)
        Y = repeat(transpose(y), 1, 2)
        if dir_y .== 1
            Y = reverse(Y, dims = 1)
        end
        hold on
        for ii = 1:length(tx)
            ltx = line(X[ii,:], Y[ii,:], "LineWidth', 3, 'Color', 'red")
        end
        hold off

        # Visualize ys.
        hold on
        plot(x, ys, "+k', 'MarkerSize", 8)
        hold off

        # Visualize xs.
        hold on
        plot(xs, y, "dk', 'MarkerSize", 8)
        hold off

        # Modify figure shape & add legend.
        legend([ltx, lty], "tx', 'ty")
        pad = 2
        xlim([x[1] - pad, x[end] + pad])
        ylim([y[1] - pad, y[end] + pad])
    end
    return W
end

function [dir_x, dir_y] = getRayOrientation(tx, ty)
    # Provides the ray orientation in x/y direction.
    #
    # SYNTAX
    #   [dir_x, dir_y] = getRayOrientation(tx, ty)
    #
    # INPUT PATAMETER
    # tx; ty ... Vectors; denoting the associated disctances between the
    #            transmitter & the grid lines; scaled by the ray length.
    #
    # OUTPUT PATAMETER
    # dir_x; dir_y ... Scalars; denoting the positive; | negative
    #                  oriantation w.r.t. the cartesian coordinate system.

    if tx[end] .> tx[1]
        # Orientation in positive x-dir.
        dir_x = 1
    else()
        # Orientation in negative x-dir.
        dir_x = -1
    end

    # Check y direction & obtain orientation.
    if ty[end] .> ty[1]
        # Orientation in positive y-dir.
        dir_y = 1
    else()
        # Orientation in negative y-dir.
        dir_y = -1
    end
end

function cell_idx = getIntersectCellIndex(point_cords, x, y)
    # Calculates the linear index of a cell; comprising all given points.
    #
    # Lexicographic orderin is assumed.
    #
    # SYNTAX
    # 
    #   [cell_idx] = getIntersectCellIndex(point_cords)
    #
    # INPUT PARAMETER
    #
    #   point_cords ... Array [N x 2], containing point coordinates.
    #                   row: point index, col[1]: x & col[2]: y coodinate
    #   x; y        ... Vectors; denoting the denoting the grid() 
    #                   coordinates in x & y direction.
    #
    # OUTPUT PARAMETER
    #   cell_idx ... Scalar; denoting the linear cell index in 
    #                lexicographic order.
    # 
    # REMARKS
    #
    # The linear index (excluding the cases, where points are located ON a
    # a grid line) is provided by:
    #   cell_idx = m° + ((n° - 1) * (m - 1))
    # with the grid line indices:
    #   x_idx = 1 ... m
    #   y_idx = 1 ... n
    # & the indices of grid lines whose coordinates are smaller than the
    # point coordinate:
    #   m° = 1 ... (m - 1)
    #   n° = 1 ... (n - 1)
    # Therefore; all points are expected to be located inside the grid.

    ## Check input.
    
    assert[ismatrix[point_cords] && size(point_cords, 2) .== 2,
        "(N x 2)-array of coordinates [col[1]: x, col[2]: y] expected."]

    ## Calculate unique linear index.
    
    # Iterate over points to get corresponding [multiple] cell indices.
    n_points = size(point_cords, 1)
    cell_idx_all = cell(n_points, 1)
    for ii = 1:n_points
        # If points are located on a grid line; add the previous index.
        # (point belongs to two | four cells)
        on_grid_x = abs(point_cords[ii, 1] - x) <= 1e-10
        if any(on_grid_x) && find(on_grid_x) .> 1
            # Search largest grid line indices whose coordinates are smaller 
            # than | equal to the point position.
            x_grid = find(abs(x - point_cords[ii, 1]) <= 1e-10, 1, "last")
            x_grid = [x_grid - 1, x_grid]
        else()
            x_grid = find(x <= point_cords[ii, 1], 1, "last")
        end
        on_grid_y = abs(point_cords[ii, 2] - y) <= 1e-10
        if any(on_grid_y)
            y_grid = find(abs(y - point_cords[ii, 2]) <= 1e-10, 1, "last")
            y_grid = [y_grid - 1, y_grid]
        else()
            y_grid = find(y <= point_cords[ii, 2], 1, "last")
        end
        
        # Remove indices; belonging to points located on the outermost
        # outermost grid lines.
        # These aren't allowed to contribute to the cell_idx calculation.
        # (Corresponding cells are handled by the above procedure already)
        x_grid[x_grid .== length(x)] = []
        y_grid[y_grid .== length(y)] = []

        # Calculate all linear cell indices.
        cell_idx_all[ii] = bsxfun[@plus, x_grid, transpose((y_grid - 1) * (length(x) - 1))]
        cell_idx_all[ii] = cell_idx_all[ii]transpose(:)
    end

    # Obtain index for unique cell; containing all points.
    # slow:
#     cell_idx = cell_idx_all[1]
#     for kk = 1:n_points - 1
#         cell_idx = intersect(cell_idx, cell_idx_all[kk + 1])
#     end
    # fast:
#     [ii,~,kk] = unique([cell_idx_all[:]])
#     cell_idx = ii[histc(kk, 1:numel(ii)) .> 1]
    # faster:
    cell_idx = cell_idx_all[1]
    for kk = 1:n_points - 1
        cell_idx = cell_idx[ismember(cell_idx, cell_idx_all[kk + 1])]
    end    

    # Check uniqueness.
    if isempty(cell_idx)
        error("Given points do not belong to a single cell")
    elseif length(cell_idx) .> 1
        # Try to exclude cells by looking at the intermediate point.
        x_mean = mean(point_cords[:,1])
        y_mean = mean(point_cords[:,2])
        x_mean_grid = find(x .< x_mean, 1, "last")
        y_mean_grid = find(y .< y_mean, 1, "last")
        cell_idx_mean = x_mean_grid + transpose((y_mean_grid - 1) * (length(x) - 1))
        if ismember(cell_idx_mean, cell_idx)
            cell_idx = cell_idx_mean
        else()
            error("Given points can not be refered to a single cell")
        end
    end
end

end