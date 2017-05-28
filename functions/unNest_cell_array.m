
function unnested_cell_array = unNest_cell_array(nested_cell_array)
% Custom function designed to "unnest" a cell array
%
% A "nested" cell array (for the sole purposes of this function) is a cell
% array where the contents of each cell are also cells (see below)
%
% nested_cell_array{1} = { 'text_kyle_wants' }
% nested_cell_array{2} = { 'more_text_kyle_wants' }
% nested_cell_array{n} = { 'even_more_text_kyle_wants' }
%
%   BECOMES
%
% unnested_cell_array{1} = 'text_kyle_wants'
% unnested_cell_array{2} = 'more_text_kyle_wants'
% unnested_cell_array{n} = 'even_more_text_kyle_wants'
%
% usageL unnested_cell_array = unNest_cell_array(nested_cell_array);
%
% Author: Kyle Kurkela, kyleakurkela@gmail.com
% Date: April, 2017

unnested_cell_array = nested_cell_array;

for i = 1:length(nested_cell_array)
    unnested_cell_array{i} = nested_cell_array{i}{:};
end