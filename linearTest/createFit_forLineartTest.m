function [fitresult, gof] = createFit_forLineartTest(axis_x, Azimuth_mean_Seq_sort)
%CREATEFIT(AXIS_X,AZIMUTH_MEAN_SEQ_SORT)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : axis_x
%      Y Output: Azimuth_mean_Seq_sort
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 23-Mar-2018 09:40:35 自动生成


%% Fit: 'Degree of Linear'.
[xData, yData] = prepareCurveData( axis_x, Azimuth_mean_Seq_sort );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'name', 'Degree of Linear' );
h = plot( fitresult, xData, yData );
legend( h, 'Measured Azimuth', 'Linear Fitting', 'Location', 'NorthEast' );
% Label axes
xlabel Num.
ylabel Azimuth(Degree)
grid on
saveas(gcf,strcat( ImagesOutputPath, 'fig_LinearAzimuth_Distribution.emf'));
saveas(gcf,strcat( ImagesOutputPath, 'fig_LinearAzimuth_Distribution.png'));


