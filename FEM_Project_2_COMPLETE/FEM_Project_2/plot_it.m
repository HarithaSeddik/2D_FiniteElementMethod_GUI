%%This is an appdesigner function that plots the mesh with colour and
%%checks necessary variables for conditions

function plot_it(app, X_1, Y_1, T, flux_x, flux_y, ...
    EdgeValue)
if(EdgeValue == 'On')
    switch app.plotType
        case 'Temperature'
            patch(app.UIAxes_2, X_1,Y_1,T,'FaceColor','interp');
            colormap(app.UIAxes_2, 'jet');
            colorbar(app.UIAxes_2);
            xlabel(app.UIAxes_2, 'X Coordinates');
            ylabel(app.UIAxes_2, 'Y Coordinates');
            title(app.UIAxes_2, ['Temperature ']);
            
        case 'Flux X'
            patch(app.UIAxes_2, X_1,Y_1,flux_x,'FaceColor','interp');
            colormap(app.UIAxes_2, 'jet');
            colorbar(app.UIAxes_2);
            xlabel(app.UIAxes_2, 'X Coordinates');
            ylabel(app.UIAxes_2, 'Y Coordinates');
            title(app.UIAxes_2, ['Flux in X Direction']);
            
        case 'Flux Y'
            patch(app.UIAxes_2, X_1,Y_1,flux_y,'FaceColor','interp');
            colormap(app.UIAxes_2, 'jet');
            colorbar(app.UIAxes_2);
            xlabel(app.UIAxes_2, 'X Coordinates');
            ylabel(app.UIAxes_2, 'Y Coordinates');
            title(app.UIAxes_2, ['Flux in Y Direction']);
    end
elseif(EdgeValue == 'Off')
    switch app.plotType
        case 'Temperature'
            patch(app.UIAxes_2, X_1,Y_1,T,'FaceColor','interp', 'EdgeColor', 'none');
            colormap(app.UIAxes_2, 'jet');
            colorbar(app.UIAxes_2);
            xlabel(app.UIAxes_2, 'X Coordinates');
            ylabel(app.UIAxes_2, 'Y Coordinates');
            title(app.UIAxes_2, ['Temperature ']);
            
        case 'Flux X'
            patch(app.UIAxes_2, X_1,Y_1,flux_x,'FaceColor','interp', 'EdgeColor', 'none');
            colormap(app.UIAxes_2, 'jet');
            colorbar(app.UIAxes_2);
            xlabel(app.UIAxes_2, 'X Coordinates');
            ylabel(app.UIAxes_2, 'Y Coordinates');
            title(app.UIAxes_2, ['Flux in X Direction']);
            
        case 'Flux Y'
            patch(app.UIAxes_2, X_1,Y_1,flux_y,'FaceColor','interp', 'EdgeColor', 'none');
            colormap(app.UIAxes_2, 'jet');
            colorbar(app.UIAxes_2);
            xlabel(app.UIAxes_2, 'X Coordinates');
            ylabel(app.UIAxes_2, 'Y Coordinates');
            title(app.UIAxes_2, ['Flux in Y Direction']);
    end
end

end