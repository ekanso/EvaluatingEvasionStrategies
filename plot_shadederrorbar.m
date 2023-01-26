function [patch,h] = plot_shadederrorbar(data,error,x_axis,lColor,aColor,lWidth)
alpha = 0.4;
x_vector = [x_axis, fliplr(x_axis)];
patch = fill(x_vector, [data + error,fliplr(data - error)], aColor);
patch.EdgeColor = 'none';
patch.FaceAlpha = alpha;
hold on;
h = plot(x_axis, data, 'color', lColor, 'LineWidth', lWidth);

end