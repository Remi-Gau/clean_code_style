% Points vector  : 3X1
% Spacing vector : 3X1
% Origin vector  : 3X1
% param is set to 1 for ijk to world
% param is set to 2 for world to ijk
function [outPoint] = ChangePtsCorSis(Origin, Spacing, Direction, inputPoints, param)

    S = zeros(3, 3);
    S(1, 1) = Spacing(1, 1);
    S(2, 2) = Spacing(1, 2);
    S(3, 3) = Spacing(1, 3);
    Origin = Origin';
    inputPoints = inputPoints';
    if param == 1  % ijk to world coord transform

        outPoint = (Direction * S * inputPoints) + Origin;

    end

    if param == 2  % World to ijk

        outPoint = (inv(Direction * S)) * (inputPoints - Origin);
        outPoint = ceil(outPoint);
    end

end
