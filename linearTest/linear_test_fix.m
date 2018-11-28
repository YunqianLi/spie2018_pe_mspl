clear 
close all
clc

%% Read from file
ImagesInputPath = '.\croped\';
ImagesOutputPath = '.\';

ROI_Range1 = [ 870, 590, 930,610 ];
ROI_Range2 = [ 870, 590, 930,610 ];
ROI_Range3 = [ 870, 590, 930,610 ];
ROI_Range4 = [ 870, 590, 930,610 ];
ROI_Range5 = [ 870, 590, 930,610 ];
ROI_Range6 = [ 870, 590, 930,610 ];


GroupNumForAnalyze = 12;


% ImagesPath = strcat( '\', ImagesPath )  ;
ImagesPath1 = '';
ImagesPath1 = strcat( ImagesInputPath, num2str( 60 + GroupNumForAnalyze - 1 ) )  ;
ImagesPath1 = strcat( ImagesPath1,'\' )  ;

ImageName1 = 'dev_reg (1).jpg'; 
ImageName2 = 'dev_reg (2).jpg'; 
ImageName3 = 'dev_reg (3).jpg'; 
ImageName4 = 'dev_reg (4).jpg'; 
ImageName5 = 'dev_reg (5).jpg'; 
ImageName6 = 'dev_reg (6).jpg'; 



I1 = imread( strcat( ImagesPath1, ImageName1) );
I2 = imread( strcat( ImagesPath1, ImageName2) );
I3 = imread( strcat( ImagesPath1, ImageName3) );
I4 = imread( strcat( ImagesPath1, ImageName4) );
I5 = imread( strcat( ImagesPath1, ImageName5) );
I6 = imread( strcat( ImagesPath1, ImageName6) );



%% RGB to YCbCr
% RGB to YCbCr
I_rgb2rcbcr1 = rgb2ycbcr( I1 );
I_rgb2rcbcr2 = rgb2ycbcr( I2 );
I_rgb2rcbcr3 = rgb2ycbcr( I3 );
I_rgb2rcbcr4 = rgb2ycbcr( I4 );
I_rgb2rcbcr5 = rgb2ycbcr( I5 );
I_rgb2rcbcr6 = rgb2ycbcr( I6  );

% Extract the Y channels
I_Y1 = I_rgb2rcbcr1( :, :, 1 );
I_Y2 = I_rgb2rcbcr2( :, :, 1 );
I_Y3 = I_rgb2rcbcr3( :, :, 1 );
I_Y4 = I_rgb2rcbcr4( :, :, 1 );
I_Y5 = I_rgb2rcbcr5( :, :, 1 );
I_Y6 = I_rgb2rcbcr6( :, :, 1 );

I_Y1_ROI = I_Y1( ROI_Range1(2):ROI_Range1(4), ROI_Range1(1):ROI_Range1(3) );
I_Y2_ROI = I_Y2( ROI_Range2(2):ROI_Range2(4), ROI_Range2(1):ROI_Range2(3) );
I_Y3_ROI = I_Y3( ROI_Range3(2):ROI_Range3(4), ROI_Range3(1):ROI_Range3(3) );
I_Y4_ROI = I_Y4( ROI_Range4(2):ROI_Range4(4), ROI_Range4(1):ROI_Range4(3) );
I_Y5_ROI = I_Y5( ROI_Range5(2):ROI_Range5(4), ROI_Range5(1):ROI_Range5(3) );
I_Y6_ROI = I_Y6( ROI_Range6(2):ROI_Range6(4), ROI_Range6(1):ROI_Range6(3) );

ROI_images_size = size( I_Y1_ROI );

% 

% figure
% subplot(2,3,1)
% imshow( uint8(I_Y1)); title('I_Y1');
% hold on; rectangle('Position',[ROI_I1(1:2) , ROI_I1(3:4)  - ROI_I1(1:2) ],'EdgeColor','r'); hold off
% 
% subplot(2,3,2)
% imshow( uint8(I_Y2)); title('I_Y2');
% hold on; rectangle('Position',[ROI_I2(1:2) , ROI_I2(3:4)  - ROI_I2(1:2) ],'EdgeColor','r'); hold off
% 
% subplot(2,3,3)
% imshow( uint8(I_Y3)); title('I_Y3');
% hold on; rectangle('Position',[ROI_I3(1:2) , ROI_I3(3:4)  - ROI_I3(1:2) ],'EdgeColor','r'); hold off
% 
% subplot(2,3,4)
% imshow( uint8(I_Y4)); title('I_Y4');
% hold on; rectangle('Position',[ROI_I4(1:2) , ROI_I4(3:4)  - ROI_I4(1:2) ],'EdgeColor','r'); hold off
% 
% subplot(2,3,5)
% imshow( uint8(I_Y5)); title('I_Y5');
% hold on; rectangle('Position',[ROI_I5(1:2) , ROI_I5(3:4)  - ROI_I5(1:2) ],'EdgeColor','r'); hold off
% 
% subplot(2,3,6)
% imshow( uint8(I_Y6)); title('I_Y6');
% hold on; rectangle('Position',[ROI_I6(1:2) , ROI_I6(3:4)  - ROI_I6(1:2) ],'EdgeColor','r'); hold off
% 

%% Calculute the Stocks Vector
% cfv = [1;2];
% fv = circpol2pol(cfv);
% G = stokes(fv);

% Stocks vector
Ix = double( I_Y1_ROI( ) );    % the standard Cartesian basis (x,y)
Iy = double( I_Y2_ROI( ) );
Ia = double( I_Y5_ROI( ) );    % Cartesian basis rotated by 45��, (a,b)
Ib = double( I_Y4_ROI( ) );
Il = double( I_Y3_ROI( ) );    % The circular basis
Ir = double( I_Y6_ROI( ) );

% stocks_I = ( (E2x + E2y) + (E2a + E2b) + (E2l + E2r) ) /3;
stocks_s0 = (Ix + Iy);
stocks_s1 = Ix - Iy;
stocks_s2 = Ia - Ib;
stocks_s3 = Il - Ir;

%% �����ƫ���ƫ���
% ԭ���������͸���ƫ��֮��Ĳ�

% ����ƫ��ȣ�����ƫ��ⲿ�ֵ�ǿ��
stocks_s_polaried = sqrt( (stocks_s1 .* stocks_s1) + (stocks_s2 .* stocks_s2) + (stocks_s3 .* stocks_s3) );

% ����ƫ���
param_DoP  = stocks_s_polaried ./ (stocks_s0 + 0001) ;

% ����ƫ��ȣ��õ���ƫ����ͼ��ֵ
S_UnPolarized = ( stocks_s0 - stocks_s_polaried );



% �������ƫ���
Ix_polaried = Ix - S_UnPolarized/2;   
Iy_polaried = Iy - S_UnPolarized/2 ;
Ia_polaried = Ia - S_UnPolarized/2 ;    % Cartesian basis rotated by 45��, (a,b)
Ib_polaried = Ib - S_UnPolarized/2 ;
Il_polaried = Il - S_UnPolarized/2 ;    % The circular basis
Ir_polaried = Ir - S_UnPolarized/2 ;


%% ����ƫ����Ϣ��Polarization info
% ��ƫ�񲿷�x�ᡢy�����
param_linear_x = ( sqrt(Ix_polaried ) ) .* ( Ia_polaried - Ib_polaried)./(0.001+abs(Ia_polaried - Ib_polaried) ) ;
param_linear_y =  sqrt(Iy_polaried  )  ;
param_Azimuth =  acos( param_linear_x ./ sqrt(param_linear_x .* param_linear_x + param_linear_y .* param_linear_y + 0.001) )  ;


% ��λ��theta
param_theta = atan( stocks_s3 ./ stocks_s2 );

% ��Բƫ�����������
param_a = sqrt( Ix_polaried );

% ��Բƫ�����������
param_b = sqrt( Iy_polaried );

% Elliptical orientation angle
param_2psi = atan( stocks_s2./ stocks_s1 );

% Degree of Linear Polarization
param_DoLP = sqrt( stocks_s1.^2 + stocks_s2.^2 ) ./ stocks_s0;

% Degree of Circular Polarization
param_DoCP = sqrt( stocks_s3.^2 ) ./ stocks_s0;

% Eccentricity angle of ellipse
param_2chi = atan( stocks_s3./ sqrt( stocks_s1.^2 + stocks_s2.^2  ) );


%% �������ڶ������µĸ���
% Լ��� ���ڿռ���ĳ���꣬�ڽ�����ͼ��ƽ�����Ͻ������ڱ��������嵽������ľ���ʱ��������Ϊ���ڲ�ͬ�����������·ֽ�ʱ��ǿ�Ⱦ��⡣
% ������  0��-90��������� �� ������  45��-135�����������ȣ�
% ������������ left circle�� right circle ����
%   �� E2x + E2y ==   E2a + E2b  == E2l + E2rb == I0;

% ������  0��-90���������ԭʼ������ 
sum_Ixy = Ix + Iy ;

% ������  45��-135���������ԭʼ������
sum_Iab = Ia + Ib;

% ������ left circle�� right circle �������ԭʼ������
sum_Ilr = Il + Ir ;

% ����������ʧֵ��������ʧԽ�����ֵԽ��
% Energy loss distribution, ELD
EnergyLoss = (abs( sum_Ixy - sum_Iab ) + abs( sum_Iab - sum_Ilr) + abs( sum_Ilr - sum_Ixy ) );

% ����������ʧ��������ԭʼ������ǿ��
stocks_I_rect = ( sum_Ixy + sum_Iab + sum_Ilr - EnergyLoss ) /3;

% ��������ʧ��ϵ��
ImbalanceRate = (EnergyLoss)./ ( sum_Ixy + sum_Iab + sum_Ilr);

% ��������ƽ��ϵ��
correct_factor = abs( 1 - ImbalanceRate );
% rect_factor = ( max(max(EnergyLoss_factor)) - EnergyLoss_factor );

% ���½���˹�п�˹����
stocks_Q_rect = ( 1 - abs(sum_Ixy - stocks_I_rect) ./ stocks_I_rect ) .* stocks_s1 ;
stocks_U_rect = ( 1 - abs(sum_Iab - stocks_I_rect) ./ stocks_I_rect ) .* stocks_s2 ;
stocks_V_rect = ( 1 - abs(sum_Ilr - stocks_I_rect) ./ stocks_I_rect ) .* stocks_s3 ;

% Rectify the Degree of Polarization
param_DoP_rect = correct_factor .* param_DoP;

% Rectify the Degree of Linear Polarization
param_DoLP_rect = correct_factor .* param_DoLP;

% Rectify the Degree of Circular Polarization
param_DoCP_rect =  correct_factor .* param_DoCP;

% ��Բƫ�����������
param_a_rect = correct_factor .* param_a;

% ��Բƫ�����������
param_b_rect = correct_factor .* param_b;

% ����ƫ��ȣ�����ƫ��ⲿ�ֵ�ǿ��
stocks_I_polaried_rect =  correct_factor .* stocks_s_polaried;




%%  Display 


% ������ƫ��ʸ����
% Calculate the vector base on degree of polarization and polarization
% angle 

% Dop_x = param_DoP .* param_a ./ sqrt( param_linear_x );
% Dop_y = param_DoP .* param_b ./ sqrt( stocks_s_polaried );
% 
% Dop_x = param_DoLP .* cos( param_Azimuth );
% Dop_y = param_DoLP .* sin( param_Azimuth );
% 
% Dop_x_rect = param_DoLP_rect .* param_linear_x ./ sqrt( param_linear_x .* param_linear_x + param_linear_y .* param_linear_y + 0.0001) ;
% Dop_y_rect = param_DoLP_rect .* param_linear_y ./ sqrt( param_linear_x .* param_linear_x + param_linear_y .* param_linear_y + 0.0001);
% 
% 
% 
% gridSize = 10;
% dis_x_line = linspace(0,1,(gridSize+1));
% dis_y_line = linspace(0,1,(gridSize+1));
% dis_x_line = uint16 ( ROI_images_size(2) * dis_x_line(1:end-1) + 1);
% dis_y_line = uint16 ( ROI_images_size(1) * dis_y_line(1:end-1) + 1);
% 
% [dis_x, dis_y] = meshgrid( dis_x_line, dis_y_line);
% 
% figure
% imshow( I3 ); title( '��ƫ�񳡣�����ǰ����ɫ������������ɫ��'); hold on; 
% quiver( ROI_Range1(1) + dis_x,ROI_Range1(2) + dis_y, Dop_x( dis_y_line, dis_x_line  ), -Dop_y( dis_y_line, dis_x_line), 'r' );
% quiver( ROI_Range1(1) + dis_x,ROI_Range1(2) + dis_y, Dop_x_rect( dis_y_line, dis_x_line  ), -Dop_y_rect( dis_y_line, dis_x_line), 'g' );
% hold off; 
% 

%% ������
% ��������Region of interest����������������ƫ��Ƕȣ���ֵ�����Լ���׼��


% ƫ���һ����
ROI_DoLP = reshape( param_DoLP_rect,[],1 );
ROI_DoLP_mean = mean( ROI_DoLP );
ROI_DoLP_std = std( ROI_DoLP );



% ��ƫ��һ����
ROI_Azimuth = reshape( param_Azimuth,[],1 );

ROI_Azimuth = ROI_Azimuth * 180 / pi;
ROI_Azimuth_mean = mean( ROI_Azimuth );
ROI_Azimuth_std = std( ROI_Azimuth );




%  figure
% plot( ROI_Azimuth_mean * ones(1,length(ROI_Azimuth)) , '.' );
%     hold on;  plot( ROI_Azimuth, '.' ); xlabel( 'Number' ); ylabel( 'Azimuth Angle (unit: degree)' );
%     str = sprintf('Average angle:%f\nStd:%f',ROI_Azimuth_mean, ROI_Azimuth_std );
%     text( length(ROI_Azimuth)/2, ROI_Azimuth_mean, str,'FontSize',14,'Color','blue', 'EdgeColor','blue', 'BackgroundColor',[0.4 0.6 0.7]);
%     hold off; 
% 
%     
%% Save and output
DoLP_Seq( :,GroupNumForAnalyze ) = ROI_DoLP;
Azimuth_Seq( :,GroupNumForAnalyze ) = ROI_Azimuth;


DoLP_mean_Seq( GroupNumForAnalyze ) = ROI_DoLP_mean;
DoLP_std_Seq( GroupNumForAnalyze  ) = ROI_DoLP_std;
Azimuth_mean_Seq( GroupNumForAnalyze ) = ROI_Azimuth_mean;
Azimuth_std_Seq( GroupNumForAnalyze ) = ROI_Azimuth_std;





%% ƫ���һ���Է���

figure

ROI_start = [ ROI_Range1(1) ROI_Range1(2) ];
ROI_end   = [ ROI_Range1(3) ROI_Range1(4)];


Dop_x_rect = param_DoLP_rect .* param_linear_x ./ sqrt( param_linear_x .* param_linear_x + param_linear_y .* param_linear_y + 0.0001) ;
Dop_y_rect = param_DoLP_rect .* param_linear_y ./ sqrt( param_linear_x .* param_linear_x + param_linear_y .* param_linear_y + 0.0001);

gridSize = 20;
dis_x_line = linspace(0,1,(gridSize+1));
dis_y_line = linspace(0,1,(gridSize+1));
dis_x_line = uint16 ( ROI_images_size(2) * dis_x_line(1:end-1) + 1);
dis_y_line = uint16 ( ROI_images_size(1) * dis_y_line(1:end-1) + 1);

[dis_x, dis_y] = meshgrid( dis_x_line, dis_y_line);

imshow( I3 ); hold on; 
rectangle('Position',[ROI_start, ROI_end - ROI_start],'EdgeColor','r', 'LineWidth', 4); 
quiver(  (ROI_Range1(1) -1 + dis_x), ( ROI_Range1(2) -1 + dis_y), Dop_x_rect( dis_y_line, dis_x_line  ), -Dop_y_rect( dis_y_line, dis_x_line), 'g' );
saveas(gcf,strcat( ImagesOutputPath, 'Image_LinearPolarAzimuthVector_high.emf'));
saveas(gcf,strcat( ImagesOutputPath, 'Image_LinearPolarAzimuthVector_high.fig'));



% ƫ���һ����
ROI_DoLP = reshape( param_DoLP_rect,[],1 );
ROI_DoLP_mean = mean( ROI_DoLP );
ROI_DoLP_std = std( ROI_DoLP );

figure
imshow( uint8(I3) ); title('Regions of Interest'); 
    hold on; rectangle('Position',[ROI_start, ROI_end - ROI_start],'EdgeColor','r'); 
    hold off
 saveas(gcf,strcat( ImagesOutputPath, 'RegionsOf_Interest.emf'));
 
figure
h = plot( (1:length(ROI_DoLP))', ROI_DoLP ,(1:length(ROI_DoLP))',  ROI_DoLP_mean * ones(length (ROI_DoLP),1) ); 
xlabel( 'Number' ); ylabel( 'DoLP(%)' ); 
hold on; 
str = sprintf('Average:%f\nStd:%f',ROI_DoLP_mean, ROI_DoLP_std );
text( length(ROI_DoLP)/10, 1.3*max(ROI_DoLP), str,'FontSize',14,'Color','blue', 'EdgeColor','blue', 'BackgroundColor',[0.4 0.6 0.7]);
legend( h, 'Mean', 'Measured Data', 'Location', 'NorthEast' );
hold off; 
saveas(gcf,strcat( ImagesOutputPath, 'fig_ROI_DOLP_Distribution.emf'));



% ��ƫ��һ����
ROI_Azimuth = reshape( param_Azimuth,[],1 );

ROI_Azimuth = ROI_Azimuth * 180 / pi;
ROI_Azimuth_mean = mean( ROI_Azimuth );
ROI_Azimuth_std = std( ROI_Azimuth );


figure
h = plot( (1:length(ROI_Azimuth))', ROI_Azimuth ,(1:length(ROI_Azimuth))',  ROI_Azimuth_mean * ones(length (ROI_Azimuth),1) ); 
    hold on;
    xlabel( 'Num.' ); ylabel( 'Azimuth Angle (degrees)' );
    str = sprintf('Average angle:%f\nStd:%f',ROI_Azimuth_mean, ROI_Azimuth_std );
%     text( h , 'Location', 'NorthEast', str,'FontSize',14,'Color','blue', 'EdgeColor','blue', 'BackgroundColor',[0.4 0.6 0.7]);
    legend( h , str, 'Location', 'NorthEast');
    hold off; 
    saveas(gcf,strcat( ImagesOutputPath, 'fig_ROI_LinearAzimuth_Distribution.emf'));




