clear 
close all
clc

%% Read from file
WorkPath = '.\';

ImagesOutputPath = strcat(WorkPath,'outputs');
if ~exist(ImagesOutputPath,'dir') 
    mkdir(ImagesOutputPath);
    str = sprintf('�����ļ���:%s',ImagesOutputPath);
    disp(str)
end
ImagesOutputPath = strcat(ImagesOutputPath,'\');


% ԭʼ�ɼ�ͼ
subpath = '.\inputs';
CapImageName1 = [subpath '\' 'dev00_crop.png' ]; 
CapImageName2 = [subpath '\' 'dev01_crop.png' ]; 
CapImageName3 = [subpath '\' 'dev02_crop.png' ]; 
CapImageName4 = [subpath '\' 'dev03_crop.png' ]; 
CapImageName5 = [subpath '\' 'dev04_crop.png' ]; 
CapImageName6 = [subpath '\' 'dev05_crop.png' ];

I1_Cap = imread( strcat( WorkPath, CapImageName1) );
I2_Cap = imread( strcat( WorkPath, CapImageName2) );
I3_Cap = imread( strcat( WorkPath, CapImageName3) );
I4_Cap = imread( strcat( WorkPath, CapImageName4) );
I5_Cap = imread( strcat( WorkPath, CapImageName5) );
I6_Cap = imread( strcat( WorkPath, CapImageName6) );

% ������ͼ
ImageDepthName = [subpath '\' 'disp.png' ]; 
ImageName1 = [subpath '\' 'dev00_reg.png' ]; 
ImageName2 = [subpath '\' 'dev01_reg.png' ]; 
ImageName3 = [subpath '\' 'dev02_reg.png' ]; 
ImageName4 = [subpath '\' 'dev03_reg.png' ]; 
ImageName5 = [subpath '\' 'dev04_reg.png' ]; 
ImageName6 = [subpath '\' 'dev05_reg.png' ]; 
I1 = imread( strcat( WorkPath, ImageName1) );
I2 = imread( strcat( WorkPath, ImageName2) );
I3 = imread( strcat( WorkPath, ImageName3) );
I4 = imread( strcat( WorkPath, ImageName4) );
I5 = imread( strcat( WorkPath, ImageName5) );
I6 = imread( strcat( WorkPath, ImageName6) );
input_images_size = size( I1 );

%% ����RIO
ROI_start = [ 370 460 ];
ROI_end   = [ 399 489 ];

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


%% Calculute the Stocks Vector
% cfv = [1;2];
% fv = circpol2pol(cfv);
% G = stokes(fv);

% Stocks vector
Ix = double( I_Y1( ) );    % the standard Cartesian basis (x,y)
Iy = double( I_Y2( ) );
Ia = double( I_Y5( ) );    % Cartesian basis rotated by 45��, (a,b)
Ib = double( I_Y4( ) );
Il = double( I_Y3( ) );    % The circular basis
Ir = double( I_Y6( ) );

% stocks_I = ( (E2x + E2y) + (E2a + E2b) + (E2l + E2r) ) /3;
stocks_s0 = Ix + Iy;
stocks_s1 = Ix - Iy;
stocks_s2 = Ia - Ib;
stocks_s3 = Il - Ir;

%% �����ƫ���ƫ���
% ԭ���������͸���ƫ��֮��Ĳ�

% ����ƫ��ȣ�����ƫ��ⲿ�ֵ�ǿ��
stocks_s_polaried = sqrt( (stocks_s1 .* stocks_s1) + (stocks_s2 .* stocks_s2) + (stocks_s3 .* stocks_s3) );

% ����ƫ���
param_DoP  = stocks_s_polaried ./ (stocks_s0 + 0.001) ;

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
%param_theta1 = atan( stocks_s3 ./ ( stocks_s2) );
param_theta = acos( sqrt((stocks_s2 .* stocks_s2 ) ./ (0.01 + 4*Ix_polaried.*Iy_polaried)) );
%  figure;imshow( (param_theta1 - min(min(param_theta1)))/((max(max(param_theta1)) - min(min(param_theta1))) ) );

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




%%  Display and Output
% Input Images
figure; imshow( I1_Cap ); title('Camera1 Image'); saveas(gcf,strcat( ImagesOutputPath, 'InputImage_cam1.emf'));
figure; imshow( I2_Cap ); title('Camera2 Image'); saveas(gcf,strcat( ImagesOutputPath, 'InputImage_cam2.emf'));
figure; imshow( I3_Cap ); title('Camera3 Image'); saveas(gcf,strcat( ImagesOutputPath, 'InputImage_cam3.emf'));
figure; imshow( I4_Cap ); title('Camera4 Image'); saveas(gcf,strcat( ImagesOutputPath, 'InputImage_cam4.emf'));
figure; imshow( I5_Cap ); title('Camera5 Image'); saveas(gcf,strcat( ImagesOutputPath, 'InputImage_cam5.emf'));
figure; imshow( I6_Cap ); title('Camera6 Image'); saveas(gcf,strcat( ImagesOutputPath, 'InputImage_cam6.emf'));


% �������ͼ
figure;
imshow(ImageDepthName); colormap default;title('Scence Depth Image'); 
saveas(gcf,strcat( ImagesOutputPath, 'Image_Depth.emf'));


% Input Images
figure; imshow( I1 ); title('Registed Image1'); saveas(gcf,strcat( ImagesOutputPath, 'RegistedImage_cam1.emf'));
figure; imshow( I2 ); title('Registed Image2'); saveas(gcf,strcat( ImagesOutputPath, 'RegistedImage_cam2.emf'));
figure; imshow( I3 ); title('Registed Image3'); saveas(gcf,strcat( ImagesOutputPath, 'RegistedImage_cam3.emf'));
figure; imshow( I4 ); title('Registed Image4'); saveas(gcf,strcat( ImagesOutputPath, 'RegistedImage_cam4.emf'));
figure; imshow( I5 ); title('Registed Image5'); saveas(gcf,strcat( ImagesOutputPath, 'RegistedImage_cam5.emf'));
figure; imshow( I6 ); title('Registed Image6'); saveas(gcf,strcat( ImagesOutputPath, 'RegistedImage_cam6.emf'));


% Stock ����ͼ
figure;
imshow( uint8(stocks_s0)); title( 'Image of  Stocks S0');
saveas(gcf,strcat( ImagesOutputPath, 'Image_ResultantIntensity.emf'));
imwrite( uint8(stocks_s0), strcat( ImagesOutputPath, 'Image_StocksS0.png') ); 


figure;
stocks_s1_output = 255*((stocks_s1 - min(min(stocks_s1))) / (max(max(stocks_s1)) - min(min(stocks_s1))) );
imshow( uint8(stocks_s1_output)); title( 'Image of  Stocks S1');
saveas(gcf,strcat( ImagesOutputPath, 'Image_StocksS1.emf'));
imwrite( uint8(stocks_s1_output), strcat( ImagesOutputPath, 'Image_StocksS1.png') ); 

figure;
stocks_s2_output = 255*((stocks_s2 - min(min(stocks_s2))) / (max(max(stocks_s2)) - min(min(stocks_s2))) );
imshow( uint8(stocks_s2_output)); title( 'Image of  Stocks S2');
saveas(gcf,strcat( ImagesOutputPath, 'Image_StocksS2.emf'));
imwrite( uint8(stocks_s2_output), strcat( ImagesOutputPath, 'Image_StocksS2.png') ); 

figure;
stocks_s3_output = 255*((stocks_s3 - min(min(stocks_s3))) / (max(max(stocks_s3)) - min(min(stocks_s3))) );
imshow( uint8(stocks_s3_output)); title( 'Image of  Stocks S3');
saveas(gcf,strcat( ImagesOutputPath, 'Image_StocksS3.emf'));
imwrite( uint8(stocks_s3_output), strcat( ImagesOutputPath, 'Image_StocksS3.png') ); 



% ƫ��⡢��ƫ������Ч���Ա�

% ��ǿ��ͼ
YCbCr(:,:,1) = uint8( stocks_s0) ;
YCbCr(:,:,2) = uint8( I_rgb2rcbcr3(:,:,2) ) ;
YCbCr(:,:,3) = uint8( I_rgb2rcbcr3(:,:,3) ) ;
rgb = ycbcr2rgb( YCbCr ) ;
clear YCbCr;

figure;
imshow( uint8(stocks_s0)); title( 'Image of Resultant Intensity (gray)');
saveas(gcf,strcat( ImagesOutputPath, 'Image_ResultantIntensity.emf'));
imwrite( uint8(stocks_s0), strcat( ImagesOutputPath, 'Image_ResultantIntensity.png') ); 

figure;
imshow( uint8(rgb)); title( 'Image of Resultant Intensity(rgb)');
saveas(gcf,strcat( ImagesOutputPath, 'Image_ResultantIntensityImage_Rgb.emf'));
imwrite( uint8(rgb), strcat( ImagesOutputPath, 'Image_ResultantIntensityImage_Rgb.png') ); 
clear rgb;

% ��ƫ���ͼ
YCbCr(:,:,1) = uint8( S_UnPolarized) ;
YCbCr(:,:,2) = uint8( I_rgb2rcbcr3(:,:,2) ) ;
YCbCr(:,:,3) = uint8( I_rgb2rcbcr3(:,:,3) ) ;
rgb = ycbcr2rgb( YCbCr ) ;
clear YCbCr;

figure;
imshow( uint8(S_UnPolarized)); title( 'Image of Unpolarized Light');
saveas(gcf,strcat( ImagesOutputPath, 'Image_UnpolarizedLight.emf'));
imwrite( uint8(S_UnPolarized), strcat( ImagesOutputPath, 'Image_UnpolarizedLight.png') ); 

figure;
imshow( uint8(rgb)); title( 'Image of Unpolarized Light(rgb)');
saveas(gcf,strcat( ImagesOutputPath, 'Image_UnpolarizedLight_rgb.emf'));
imwrite( uint8(rgb), strcat( ImagesOutputPath, 'Image_UnpolarizedLight_rgb.png') ); 
clear rgb;

% ƫ���ͼ
YCbCr(:,:,1) = uint8( stocks_s_polaried) ;
YCbCr(:,:,2) = uint8( I_rgb2rcbcr3(:,:,2) ) ;
YCbCr(:,:,3) = uint8( I_rgb2rcbcr3(:,:,3) ) ;
rgb = ycbcr2rgb( YCbCr ) ;
clear YCbCr;

figure;
imshow( uint8(stocks_s_polaried)); title( 'Image of Polarized Light');
saveas(gcf,strcat( ImagesOutputPath, 'Image_PolarizedLight.emf'));
imwrite( uint8(stocks_s_polaried), strcat( ImagesOutputPath, 'Image_PolarizedLight.png') ); 

figure;
imshow( uint8(rgb)); title( 'Image of Polarized Light(rgb)');
saveas(gcf,strcat( ImagesOutputPath, 'Image_PolarizedLight_rgb.emf'));
imwrite( uint8(rgb), strcat( ImagesOutputPath, 'Image_PolarizedLight_rgb.png') ); 

% ������
YCbCr(:,:,1) = uint8( stocks_I_polaried_rect) ;
YCbCr(:,:,2) = uint8( I_rgb2rcbcr3(:,:,2) ) ;
YCbCr(:,:,3) = uint8( I_rgb2rcbcr3(:,:,3) ) ;
rgb = ycbcr2rgb( YCbCr ) ;
clear YCbCr;

figure;
imshow( uint8(stocks_I_polaried_rect)); title( 'Image of Polarized Light(rect)');
saveas(gcf,strcat( ImagesOutputPath, 'Image_PolarizedLight_rect.emf'));
imwrite( uint8(stocks_I_polaried_rect), strcat( ImagesOutputPath, 'Image_PolarizedLight_rect.png') ); 

figure;
imshow( uint8(rgb)); title( 'Image of Polarized Light (rect,rgb)');
saveas(gcf,strcat( ImagesOutputPath, 'Image_PolarizedLight_rect_rgb.emf'));
imwrite( uint8(rgb), strcat( ImagesOutputPath, 'Image_PolarizedLight_rect_rgb.png') ); 




% ǿ��ƽ��Լ�������Ա�ͼ
figure; 
imshow( uint8(255*ImbalanceRate) ); title('Imbalance Rate');
saveas(gcf,strcat( ImagesOutputPath, 'Image_ImbalanceRate.emf'));
imwrite( uint8(255*ImbalanceRate), strcat( ImagesOutputPath, 'Image_ImbalanceRate.png') ); 

figure; 
imshow( uint8(255*(1-ImbalanceRate)) ); title('Intensity Balance Rate');
saveas(gcf,strcat( ImagesOutputPath, 'Image_IntensityBalanceRate.emf'));
imwrite( uint8(255*(1-ImbalanceRate)), strcat( ImagesOutputPath, 'Image_IntensityBalanceRate.png') ); 


figure;
imshow( uint8(stocks_s_polaried) ); title('Image before Stocks S0 Imbalance Corrected');
saveas(gcf,strcat( ImagesOutputPath, 'Image_StocksS0_beforeImbalanceCorrected.emf'));
imwrite( uint8(stocks_s_polaried), strcat( ImagesOutputPath, 'Image_StocksS0_beforeImbalanceCorrected.png') ); 

figure; 
imshow( uint8(stocks_I_polaried_rect) ); title('Image of Stocks S0 After Imbalance Corrected');
saveas(gcf,strcat( ImagesOutputPath, 'Image_StocksS0_ImbalanceCorrected.emf'));
imwrite( uint8(stocks_I_polaried_rect), strcat( ImagesOutputPath, 'Image_StocksS0_ImbalanceCorrected.png') ); 



% ƫ�������ʾ
% Polarization parameters display
figure; 
imshow(  uint8(255*(( param_DoP_rect - min(min(param_DoP_rect))))/( max(max(param_DoP_rect)) - min(min(param_DoP_rect))) ) ); title('ƫ��ȣ�Degree of Polarization��');
saveas(gcf,strcat( ImagesOutputPath, 'Image_DOP.emf'));
imwrite( uint8(255*(( param_DoP_rect - min(min(param_DoP_rect))))/( max(max(param_DoP_rect)) - min(min(param_DoP_rect))) ) ,...
    strcat( ImagesOutputPath, 'Image_DOP.png') ); 

figure; 
imshow( uint8(255*(( param_2psi - min(min(param_2psi))))/( max(max(param_2psi)) - min(min(param_2psi))) ) ); title('Images of Ellipse Orientation (2psi)');
saveas(gcf,strcat( ImagesOutputPath, 'Image_EllipPola_2psi.emf'));
imwrite( uint8(255*(( param_2psi - min(min(param_2psi))))/( max(max(param_2psi)) - min(min(param_2psi))) ) ,...
    strcat( ImagesOutputPath, 'Image_EllipPola_2psi.png') ); 

figure; 
imshow( uint8(255*(( param_2chi - min(min(param_2chi))))/( max(max(param_2chi)) - min(min(param_2chi))) ) ); title('Images of Ellipse Ellipticity (2chi)');
saveas(gcf,strcat( ImagesOutputPath, 'Image_EllipPola_2chi.emf'));
imwrite( uint8(255*(( param_2chi - min(min(param_2chi))))/( max(max(param_2chi)) - min(min(param_2chi))) ) ,...
    strcat( ImagesOutputPath, 'Image_EllipPola_2chi.png') ); 

figure
imshow( uint8( 255 * ( (  1 - abs(correct_factor.*param_Azimuth + pi/2)/pi) ))) ; title('Image of Linear Plorization Azimuth');
saveas(gcf,strcat( ImagesOutputPath, 'Image_LinearPolarAzimuth.emf'));
imwrite( uint8( 255 * ( (  1 - abs(correct_factor.*param_Azimuth + pi/2)/pi) )) ,...
    strcat( ImagesOutputPath, 'Image_LinearPolarAzimuth.png') ); 

figure
imshow( uint8( 255 * param_DoP_rect .* ( (param_theta - min(min(param_theta)))/((max(max(param_theta)) - min(min(param_theta))) ) ) )) ; 
title('Image of Elipse Polarization Phase');
saveas(gcf,strcat( ImagesOutputPath, 'Image_ElipsePolarizationPhase.emf'));
imwrite( uint8( 255 * ( (param_theta - min(min(param_theta)))/((max(max(param_theta)) - min(min(param_theta))) ) )) ,...
    strcat( ImagesOutputPath, 'Image_ElipsePolarizationPhase.png') ); 

%  figure;imshow(  );
%  

% ������ƫ��ʸ����
% Calculate the vector base on degree of polarization and polarization
% angle 

% Dop_x = param_DoP .* param_a ./ sqrt( param_linear_x );
% Dop_y = param_DoP .* param_b ./ sqrt( stocks_s_polaried );

Dop_x = param_DoLP .* cos( param_Azimuth );
Dop_y = param_DoLP .* sin( param_Azimuth );

Dop_x_rect = param_DoLP_rect .* param_linear_x ./ sqrt( param_linear_x .* param_linear_x + param_linear_y .* param_linear_y + 0.0001) ;
Dop_y_rect = param_DoLP_rect .* param_linear_y ./ sqrt( param_linear_x .* param_linear_x + param_linear_y .* param_linear_y + 0.0001);

h_quiver1 = figure;
gridSize = 20;
dis_x_line = linspace(0,1,(gridSize+1));
dis_y_line = linspace(0,1,(gridSize+1));
dis_x_line = uint16 ( input_images_size(2) * dis_x_line(1:end-1) + 1);
dis_y_line = uint16 ( input_images_size(1) * dis_y_line(1:end-1) + 1);

[dis_x, dis_y] = meshgrid( dis_x_line, dis_y_line);

imshow( I3 ); title('Image_LinearPolarAzimuthVector(red),corrected(green)'); hold on; 
quiver( dis_x, dis_y, Dop_x( dis_y_line, dis_x_line  ), -Dop_y( dis_y_line, dis_x_line), 'r' );
quiver( dis_x, dis_y, Dop_x_rect( dis_y_line, dis_x_line  ), -Dop_y_rect( dis_y_line, dis_x_line), 'g' );
saveas(gcf,strcat( ImagesOutputPath, 'Image_LinearPolarAzimuthVector_compare.emf'));


% ��ͷ�ֲ��Ŵ�Ŵ�
h_quiver1 = figure;
gridSize = 200;
dis_x_line = linspace(0,1,(gridSize+1));
dis_y_line = linspace(0,1,(gridSize+1));
dis_x_line = uint16 ( input_images_size(2) * dis_x_line(1:end-1) + 1);
dis_y_line = uint16 ( input_images_size(1) * dis_y_line(1:end-1) + 1);
[dis_x, dis_y] = meshgrid( dis_x_line, dis_y_line);

figure
ROI_start1 = [ 370 460 ];
ROI_end1   = [ 399 489 ];

imshow( I3 ); hold on; 
rectangle('Position',[ROI_start1, ROI_end1 - ROI_start1],'EdgeColor','r', 'LineWidth', 4); 
quiver( dis_x, dis_y, Dop_x_rect( dis_y_line, dis_x_line  ), -Dop_y_rect( dis_y_line, dis_x_line), 'g' );
saveas(gcf,strcat( ImagesOutputPath, 'Image_LinearPolarAzimuthVector_high.emf'));
saveas(gcf,strcat( ImagesOutputPath, 'Image_LinearPolarAzimuthVector_high.fig'));




% ƫ��ͼ��
% Polarization Enhanced Image
Strength = 3;
PolarizationImage = ( double(I_Y6) .* ( 1 + Strength*param_DoP_rect ));
% PolarizationImages = PolarizationImages / max( max(PolarizationImages);
PolarizationImage_YCbCr(:,:,1) = uint8( PolarizationImage) ;
PolarizationImage_YCbCr(:,:,2) = uint8( I_rgb2rcbcr3(:,:,2) ) ;
PolarizationImage_YCbCr(:,:,3) = uint8( I_rgb2rcbcr3(:,:,3) ) ;
PolarizationImage_rgb = ycbcr2rgb( PolarizationImage_YCbCr ) ;
clear PolarizationImage_YCbCr;

figure; 
imshow( uint8( PolarizationImage_rgb ) ); title('Polarization Enhanced Image');
saveas(gcf,strcat( ImagesOutputPath, 'Image_PolarizationEnhanced.emf'));
imwrite( uint8(PolarizationImage_rgb), strcat( ImagesOutputPath, 'Image_PolarizationEnhanced.png') ); 

figure; imshow( uint8(255*param_DoLP_rect - 255*param_DoCP_rect) ); title('LinearPolarizedImage');
saveas(gcf,strcat( ImagesOutputPath, 'Image_LinearPolarized.emf'));
imwrite( uint8(255*param_DoLP_rect - 255*param_DoCP_rect), strcat( ImagesOutputPath, 'Image_LinearPolarized.png') ); 

figure;
imshow( uint8(255*param_DoCP_rect) ); title('Circle Polarized Image');
saveas(gcf,strcat( ImagesOutputPath, 'Image_CirclePolarizedImage.emf'));
imwrite( uint8(255*param_DoCP_rect), strcat( ImagesOutputPath, 'CirclePolarizedImage.png') ); 



%% ������
% ��������Region of interest����������������ƫ��Ƕȣ���ֵ�����Լ���׼��
% ROI_start = [ 360 460 ];
% ROI_end   = [ 400 500 ];

% ƫ���һ����
ROI_DoLP = reshape( param_DoLP_rect( (ROI_start(2):1:ROI_end(2) ), (ROI_start(1):1:ROI_end(1)) ),[],1 );
ROI_DoLP_mean = mean( ROI_DoLP );
ROI_DoLP_std = std( ROI_DoLP );

figure
imshow( uint8(I3) ); title('Regions of Interest'); 
    hold on; rectangle('Position',[ROI_start, ROI_end - ROI_start],'EdgeColor','r'); 
    hold off
 saveas(gcf,strcat( ImagesOutputPath, 'RegionsOf_Interest.emf'));
 
figure
h = plot( ROI_DoLP_mean * ones(1,length(ROI_DoLP)) , '.' );
    hold on;  plot( ROI_DoLP, '.' ); xlabel( 'Number' ); ylabel( 'DoLP(%)' ); axis([0 length(ROI_DoLP) 0 1 ]);
    str = sprintf('Average:%f\nStd:%f',ROI_DoLP_mean, ROI_DoLP_std );
%     text( length(ROI_DoLP)/10, 1.3*max(ROI_DoLP), str,'FontSize',14,'Color','blue', 'EdgeColor','blue', 'BackgroundColor',[0.4 0.6 0.7]);
    legend( h , str, 'Location', 'NorthEast');
    hold off; 
    saveas(gcf,strcat( ImagesOutputPath, 'fig_ROI_DOLP_Distribution.emf'));



% ��ƫ��һ����
ROI_Azimuth = reshape( param_Azimuth( (ROI_start(2):1:ROI_end(2) ), (ROI_start(1):1:ROI_end(1)) ),[],1 );

ROI_Azimuth = ROI_Azimuth * 180 / pi;
ROI_Azimuth_mean = mean( ROI_Azimuth );
ROI_Azimuth_std = std( ROI_Azimuth );


figure
h = plot( ROI_Azimuth_mean * ones(1,length(ROI_Azimuth)) , '.' );
    hold on;  plot( ROI_Azimuth, '.' ); xlabel( 'Num.' ); ylabel( 'Azimuth Angle (unit: degree)' );
    str = sprintf('Average angle:%f\nStd:%f',ROI_Azimuth_mean, ROI_Azimuth_std );
%     text( h , 'Location', 'NorthEast', str,'FontSize',14,'Color','blue', 'EdgeColor','blue', 'BackgroundColor',[0.4 0.6 0.7]);
    legend( h , str, 'Location', 'NorthEast');
    hold off; 
    saveas(gcf,strcat( ImagesOutputPath, 'fig_ROI_LinearAzimuth_Distribution.emf'));







