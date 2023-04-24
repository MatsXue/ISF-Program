%% 2D Principle Stiffness 
k11=[340 610 520 485;
     380 530 425 360;
     430 420 350 350;
     260 475 510 435
     ]';

%% Cubic Interpolation
x_1=[50 400 750 980];
y_1=[22.4 372.5 722.5 1009];
x_map=-10:25:1040;
y_map=-35:25:1015;

k1_spine=zeros(43,4);
k1_map=zeros(43,43);

% Get 4 y-direction spines on the map
for i=1:4
    k1_spine(1:43,i)=spline(x_1,k11(:,i),x_map);
end

% Get rows with 4 columns on the map
for i=1:43
    k1_map(i,1:43)=spline(y_1, k1_spine(i,:), y_map);
end

%% Display the Maps
k1_max=max(max(k1_map));
k1_min=min(min(k1_map));
for row = 1:43
    for column = 1:43
        plot_x1=-22.5+25*row;
        plot_y1=-47.5+25*column;
        plot_x2=2.5+25*row;
        plot_y2=-22.5+25*column;
        
        plot_x=[plot_x1 plot_x2 plot_x2 plot_x1 plot_x1];
        plot_y=[plot_y1 plot_y1 plot_y2 plot_y2 plot_x1];
        
        if k1_map(row,column)<k1_min+0.2*(k1_max-k1_min)
            p=fill(plot_x,plot_y,[1 1 1]);
        elseif k1_map(row,column)<k1_min+0.4*(k1_max-k1_min)
            p=fill(plot_x,plot_y,[1 0.75 0.75]);
        elseif k1_map(row,column)<k1_min+0.6*(k1_max-k1_min)
            p=fill(plot_x,plot_y,[1 0.5 0.5]);
        elseif k1_map(row,column)<k1_min+0.8*(k1_max-k1_min)
            p=fill(plot_x,plot_y,[1 0.25 0.25]);
        else
            p=fill(plot_x,plot_y,[1 0 0]);
        end
        p.EdgeColor=[1 1 1];
        xlim tight;
        ylim tight;
        hold on
    end
end
hold off
%title('Principle stiffness values K_{YY}')
xlabel('workpiece X axis')
ylabel('workpiece Y axis')
view(90,90)
set(gca,'XAxisLocation','bottom','YAxisLocation','right');
