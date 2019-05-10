
function [new_ref] = animation(r1T, r2T, r3T, r4T, x_ref, y_ref, lineObj, mode)
if mode == 'r'
    %set target position and pause
    set(lineObj.h1T,'xdata',[x_ref x_ref+r1T(1)],'ydata',[y_ref y_ref+r1T(2)])
    set(lineObj.h2T,'xdata',[x_ref+r1T(1) x_ref+r2T(1)],'ydata',[y_ref+r1T(2) y_ref+r2T(2)])
    set(lineObj.h3T,'xdata',[x_ref+r2T(1) x_ref+r3T(1)],'ydata',[y_ref+r2T(2) y_ref+r3T(2)])
    set(lineObj.h4T,'xdata',[x_ref+r3T(1) x_ref+r4T(1)],'ydata',[y_ref+r3T(2) y_ref+r4T(2)])
    drawnow;

    new_ref = x_ref+r4T(1);
elseif mode == 'f'
    %set target position and pause
    set(lineObj.h1T,'xdata',[x_ref+r1T(1) x_ref+r2T(1)],'ydata',[y_ref+r1T(2) y_ref+r2T(2)])
    set(lineObj.h2T,'xdata',[x_ref+r2T(1) x_ref+r3T(1)],'ydata',[y_ref+r2T(2) y_ref+r3T(2)])
    set(lineObj.h3T,'xdata',[x_ref+r4T(1) x_ref+r3T(1)],'ydata',[y_ref+r4T(2) y_ref+r3T(2)])
    set(lineObj.h4T,'xdata',[x_ref x_ref+r4T(1)],'ydata',[y_ref y_ref+r4T(2)])
    drawnow;

    new_ref = x_ref+r1T(1);
end
end
