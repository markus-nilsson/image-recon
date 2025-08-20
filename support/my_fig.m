function h = my_fig(tag)
%  Return a figure with a given tag, creating it if needed.
%
%   h = my_fig(tag) returns a handle to a figure that has the
%   specified Tag. If such a figure does not exist, it creates one.

% Try to find an existing figure with this tag
fig = findobj('Type','figure','Tag',tag);

if isempty(fig) || ~isvalid(fig)
    % Create a new figure with the given tag
    h = figure('Tag',tag,'Name',tag,'NumberTitle','off');
else
    % If multiple matches exist, pick the first
    h = fig(1);
end

set(h, 'Color', 'white');
