function [IJKneighbors, interpWeighs] = idw3dInterp(fiberSet, volume)

    maxSlice_volume = volume.size(3);
    volume.SliceRASidx = zeros(1, maxSlice_volume); % (1XN)
    for i = 1:maxSlice_volume
        curPoint = ...
            ChangePtsCorSis(volume.origin, volume.spacing, volume.direction, [0 0 i], 1);
        volume.SliceRASidx(1, i) = curPoint(3);
    end

    maxX_volume = volume.size(1);
    volume.xRASidx = zeros(1, maxX_volume); % (1XN)
    for i = 1:maxX_volume
        curPoint = ...
            ChangePtsCorSis(volume.origin, volume.spacing, volume.direction, [i 0 0], 1);
        volume.xRASidx(1, i) = curPoint(1);
    end

    maxY_volume = volume.size(2);
    volume.yRASidx = zeros(1, maxY_volume); % (1XN)
    for i = 1:maxY_volume
        curPoint = ...
            ChangePtsCorSis(volume.origin, volume.spacing, volume.direction, [0 i 0], 1);
        volume.yRASidx(1, i) = curPoint(2);
    end

    % Loop through all fiber nodes (this operation wil take some time)

    IJKneighbors = struct(repmat(struct, [1 length(fiberSet)]));
    interpWeighs = IJKneighbors;

    for i = 1:length(fiberSet)

        IJKneighbors(i).mat = zeros(length(fiberSet(i).matrix), 3, 8);
        interpWeighs(i).mat = zeros(length(fiberSet(i).matrix), 8);

        for j = 1:length(fiberSet(i).matrix)

            curCoord = fiberSet(i).matrix(j, :);

            % Get closest 3 slices (along-z)
            [~, qz] = sort(sum(abs(bsxfun(@minus, volume.SliceRASidx, curCoord(3))), 1));
            qz = sort(qz(1:3));
            % Get closest 3 slices (along-x)
            [~, qx] = sort(sum(abs(bsxfun(@minus, volume.xRASidx, curCoord(1))), 1));
            qx = sort(qx(1:3));
            % Get closest 3 slices (along-y)
            [~, qy] = sort(sum(abs(bsxfun(@minus, volume.yRASidx, curCoord(2))), 1));
            qy = sort(qy(1:3));
            slice_points_ijk = zeros(3, 3, 9);
            slice_points_RAS = zeros(3, 3, 9);

            for k = 1:3
                for l = 1:3
                    for m = 1:3
                        slice_points_RAS(k, :, m + (l - 1) * 3) = [volume.xRASidx(qx(l)), ...
                                                                   volume.yRASidx(qy(m)), ...
                                                                   volume.SliceRASidx(qz(k))];
                        slice_points_ijk(k, :, m + (l - 1) * 3) = [qx(l) qy(m) qz(k)];
                    end
                end
            end

            % Search for 8 nodes.
            normHolder = zeros(27, 3);
            iter = 1;
            while iter <= 27
                for k = 1:3 % In each slice
                    for l = 1:9 % For each point
                        normHolder(iter, 1) = norm(slice_points_RAS(k, :, l) - curCoord);
                        normHolder(iter, 2) = k;
                        normHolder(iter, 3) = l;
                        iter = iter + 1;
                    end
                end
            end

            normHolder_sorted = sortrows(normHolder);
            normHolder_sorted = normHolder_sorted(1:8, :); % Here we have 8 of them

            weighs = power(normHolder_sorted(:, 1), -1);
            totweigh = sum(weighs);
            weighs = weighs ./ totweigh;

            for m = 1:8
                IJKneighbors(i).mat(j, :, m) = slice_points_ijk(normHolder_sorted(m, 2), ...
                                                                :, ...
                                                                normHolder_sorted(m, 3));
                interpWeighs(i).mat(j, m) = weighs(m);
            end

        end
    end

end
