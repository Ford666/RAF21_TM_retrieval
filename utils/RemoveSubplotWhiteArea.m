function [] = RemoveSubplotWhiteArea(sub_row, sub_col, current_row, current_col)
% 设置OuterPosition
sub_axes_x = current_col*1/sub_col - 1/sub_col;
sub_axes_y = 1-current_row*1/sub_row; % y是从上往下的
sub_axes_w = 1/sub_col;
sub_axes_h = 1/sub_row;
set(gca, 'OuterPosition', [sub_axes_x, sub_axes_y, sub_axes_w, sub_axes_h]); % 重设OuterPosition

% TightInset的位置
inset_vectior = get(gca, 'TightInset');
inset_x = inset_vectior(1);
inset_y = inset_vectior(2);
inset_w = inset_vectior(3);
inset_h = inset_vectior(4);

% OuterPosition的位置
outer_vector = get(gca, 'OuterPosition');
pos_new_x = outer_vector(1) + inset_x; % 将Position的原点移到到TightInset的原点
pos_new_y = outer_vector(2) + inset_y;
pos_new_w = outer_vector(3) - inset_w - inset_x; % 重设Position的宽
pos_new_h = outer_vector(4) - inset_h - inset_y; % 重设Position的高

% 重设Position
set(gca, 'Position', [pos_new_x, pos_new_y, pos_new_w, pos_new_h]);