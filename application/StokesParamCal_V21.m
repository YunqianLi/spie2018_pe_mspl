clear 
close all
clc

%% Read from file
ImagesPath = 'F:\偏振采集系统\仿真\偏振\StokesParam\素材\论文实验项素材\3论文-应用1\';
ImagesOutputPath = ImagesPath ;

% 判断脚本和图片素材是否在一个文件夹，不在一个文件夹则不执行，避免覆盖
d1 = dir(ImagesPath);  
d2 = dir('.\');  
if ( ~isequal(d1,d2) )
    error('为了避免混淆，请将本脚本放在素材所在的文件夹执行！！！');
end


ImageName1 = 'dev00_reg.png'; 
ImageName2 = 'dev01_reg.png'; 
ImageName3 = 'dev02_reg.png'; 
ImageName4 = 'dev03_reg.png'; 
ImageName5 = 'dev04_reg.png'; 
ImageName6 = 'dev05_reg.png'; 

I1 = imread( strcat( ImagesPath, ImageName1) );
I2 = imread( strcat( ImagesPath, ImageName2) );
I3 = imread( strcat( ImagesPath, ImageName3) );
I4 = imread( strcat( ImagesPath, ImageName4) );
I5 = imread( strcat( ImagesPath, ImageName5) );
I6 = imread( strcat( ImagesPath, ImageName6) );

input_images_size = size( I1 );

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
Ia = double( I_Y5( ) );    % Cartesian basis rotated by 45°, (a,b)
Ib = double( I_Y4( ) );
Il = double( I_Y3( ) );    % The circular basis
Ir = double( I_Y6( ) );

% stocks_I = ( (E2x + E2y) + (E2a + E2b) + (E2l + E2r) ) /3;
stocks_s0 = Ix + Iy;
stocks_s1 = Ix - Iy;
stocks_s2 = Ia - Ib;
stocks_s3 = Il - Ir;

%% 分离非偏振和偏振光
% 原理：总能量和各个偏振之间的差

% 根据偏振度，计算偏振光部分的强度
stocks_s_polaried = sqrt( (stocks_s1 .* stocks_s1) + (stocks_s2 .* stocks_s2) + (stocks_s3 .* stocks_s3) );

% 计算偏振度
param_DoP  = stocks_s_polaried ./ (stocks_s0 + 0.001) ;

% 根据偏振度，得到非偏振光的图像值
S_UnPolarized = ( stocks_s0 - stocks_s_polaried );



% 分离出纯偏振光
Ix_polaried = Ix - S_UnPolarized/2;   
Iy_polaried = Iy - S_UnPolarized/2 ;
Ia_polaried = Ia - S_UnPolarized/2 ;    % Cartesian basis rotated by 45°, (a,b)
Ib_polaried = Ib - S_UnPolarized/2 ;
Il_polaried = Il - S_UnPolarized/2 ;    % The circular basis
Ir_polaried = Ir - S_UnPolarized/2 ;


%% 计算偏振信息，Polarization info
% 线偏振部分x轴、y轴分量
param_linear_x = ( sqrt(Ix_polaried ) ) .* ( Ia_polaried - Ib_polaried)./(0.001+abs(Ia_polaried - Ib_polaried) ) ;
param_linear_y =  sqrt(Iy_polaried  )  ;
param_Azimuth =  acos( param_linear_x ./ sqrt(param_linear_x .* param_linear_x + param_linear_y .* param_linear_y + 0.001) )  ;


% 相位差theta
param_theta = atan( stocks_s3 ./ stocks_s2 );

% 椭圆偏振参数，长轴
param_a = sqrt( Ix_polaried );

% 椭圆偏振参数，短轴
param_b = sqrt( Iy_polaried );

% Elliptical orientation angle
param_2psi = atan( stocks_s2./ stocks_s1 );

% Degree of Linear Polarization
param_DoLP = sqrt( stocks_s1.^2 + stocks_s2.^2 ) ./ stocks_s0;

% Degree of Circular Polarization
param_DoCP = sqrt( stocks_s3.^2 ) ./ stocks_s0;

% Eccentricity angle of ellipse
param_2chi = atan( stocks_s3./ sqrt( stocks_s1.^2 + stocks_s2.^2  ) );


%% 矫正由于对齐误差导致的干扰
% 约束项： 对于空间中某坐标，在接收器图像平面相距较近，大于被拍摄物体到成像面的距离时，我们认为光在不同的正交基底下分解时，强度均衡。
% 正交组  0°-90°的总能量 与 正交组  45°-135°的总能量相等，
% 并且与正交组 left circle和 right circle 等于
%   即 E2x + E2y ==   E2a + E2b  == E2l + E2rb == I0;

% 正交组  0°-90°所捕获的原始总能量 
sum_Ixy = Ix + Iy ;

% 正交组  45°-135°所捕获的原始总能量
sum_Iab = Ia + Ib;

% 正交组 left circle和 right circle 所捕获的原始总能量
sum_Ilr = Il + Ir ;

% 计算能量损失值，能量损失越大参数值越大
% Energy loss distribution, ELD
EnergyLoss = (abs( sum_Ixy - sum_Iab ) + abs( sum_Iab - sum_Ilr) + abs( sum_Ilr - sum_Ixy ) );

% 根据能量损失量，估计原始的能量强度
stocks_I_rect = ( sum_Ixy + sum_Iab + sum_Ilr - EnergyLoss ) /3;

% 计算能量失衡系数
ImbalanceRate = (EnergyLoss)./ ( sum_Ixy + sum_Iab + sum_Ilr);

% 计算能量平衡系数
correct_factor = abs( 1 - ImbalanceRate );
% rect_factor = ( max(max(EnergyLoss_factor)) - EnergyLoss_factor );

% 重新矫正斯托克斯参数
stocks_Q_rect = ( 1 - abs(sum_Ixy - stocks_I_rect) ./ stocks_I_rect ) .* stocks_s1 ;
stocks_U_rect = ( 1 - abs(sum_Iab - stocks_I_rect) ./ stocks_I_rect ) .* stocks_s2 ;
stocks_V_rect = ( 1 - abs(sum_Ilr - stocks_I_rect) ./ stocks_I_rect ) .* stocks_s3 ;

% Rectify the Degree of Polarization
param_DoP_rect = correct_factor .* param_DoP;

% Rectify the Degree of Linear Polarization
param_DoLP_rect = correct_factor .* param_DoLP;

% Rectify the Degree of Circular Polarization
param_DoCP_rect =  correct_factor .* param_DoCP;

% 椭圆偏振参数，长轴
param_a_rect = correct_factor .* param_a;

% 椭圆偏振参数，短轴
param_b_rect = correct_factor .* param_b;

% 根据偏振度，计算偏振光部分的强度
stocks_I_polaried_rect =  correct_factor .* stocks_s_polaried;




%%  Display and Output
% Input Images
figure; imshow( I1 ); title('Camera1 Image'); saveas(gcf,strcat( ImagesOutputPath, 'InputImage_cam1.emf'));
figure; imshow( I2 ); title('Camera2 Image'); saveas(gcf,strcat( ImagesOutputPath, 'InputImage_cam2.emf'));
figure; imshow( I3 ); title('Camera3 Image'); saveas(gcf,strcat( ImagesOutputPath, 'InputImage_cam3.emf'));
figure; imshow( I4 ); title('Camera4 Image'); saveas(gcf,strcat( ImagesOutputPath, 'InputImage_cam4.emf'));
figure; imshow( I5 ); title('Camera5 Image'); saveas(gcf,strcat( ImagesOutputPath, 'InputImage_cam5.emf'));
figure; imshow( I6 ); title('Camera6 Image'); saveas(gcf,strcat( ImagesOutputPath, 'InputImage_cam6.emf'));



% 偏振光、非偏振光分离效果对比

% 总强度图
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

% 非偏振光图
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

% 偏振光图
YCbCr(:,:,1) = uint8( stocks_I_polaried_rect) ;
YCbCr(:,:,2) = uint8( I_rgb2rcbcr3(:,:,2) ) ;
YCbCr(:,:,3) = uint8( I_rgb2rcbcr3(:,:,3) ) ;
rgb = ycbcr2rgb( YCbCr ) ;
clear YCbCr;

figure;
imshow( uint8(stocks_I_polaried_rect)); title( 'Image of Polarized Light');
saveas(gcf,strcat( ImagesOutputPath, 'Image_PolarizedLight.emf'));
imwrite( uint8(stocks_I_polaried_rect), strcat( ImagesOutputPath, 'Image_PolarizedLight.png') ); 

figure;
imshow( uint8(rgb)); title( 'Image of Polarized Light(rgb)');
saveas(gcf,strcat( ImagesOutputPath, 'Image_PolarizedLight_rgb.emf'));
imwrite( uint8(rgb), strcat( ImagesOutputPath, 'Image_PolarizedLight_rgb.png') ); 



% 强度平衡约束矫正对比图
figure; 
imshow( uint8(255*ImbalanceRate) ); title('Imbalance Rate');
saveas(gcf,strcat( ImagesOutputPath, 'Image_ImbalanceRate.emf'));
imwrite( uint8(255*ImbalanceRate), strcat( ImagesOutputPath, 'Image_ImbalanceRate.png') ); 

figure;
imshow( uint8(stocks_s_polaried) ); title('Image before Stocks S0 Imbalance Corrected');
saveas(gcf,strcat( ImagesOutputPath, 'Image_StocksS0_beforeImbalanceCorrected.emf'));
imwrite( uint8(stocks_s_polaried), strcat( ImagesOutputPath, 'Image_StocksS0_beforeImbalanceCorrected.png') ); 

figure; 
imshow( uint8(stocks_I_polaried_rect) ); title('Image of Stocks S0 After Imbalance Corrected');
saveas(gcf,strcat( ImagesOutputPath, 'Image_StocksS0_ImbalanceCorrected.emf'));
imwrite( uint8(stocks_I_polaried_rect), strcat( ImagesOutputPath, 'Image_StocksS0_ImbalanceCorrected.png') ); 



% 偏振参数显示
% Polarization parameters display
figure; 
imshow(  uint8(255*(( param_DoP_rect - min(min(param_DoP_rect))))/( max(max(param_DoP_rect)) - min(min(param_DoP_rect))) ) ); title('偏振度（Degree of Polarization）');
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
saveas(gcf,strcat( ImagesOutputPath, 'Image_EllipPola_2psi.emf'));
imwrite( uint8(255*(( param_2chi - min(min(param_2chi))))/( max(max(param_2chi)) - min(min(param_2chi))) ) ,...
    strcat( ImagesOutputPath, 'Image_EllipPola_2chi.png') ); 

figure
imshow( uint8( 255 * ( (  1 - abs(correct_factor.*param_Azimuth + pi/2)/pi) ))) ; title('Image of Linear Plorization Azimuth');
saveas(gcf,strcat( ImagesOutputPath, 'Image_LinearPolarAzimuth.emf'));
imwrite( uint8( 255 * ( (  1 - abs(correct_factor.*param_Azimuth + pi/2)/pi) )) ,...
    strcat( ImagesOutputPath, 'Image_LinearPolarAzimuth.png') ); 


% 绘制线偏振矢量场
% Calculate the vector base on degree of polarization and polarization
% angle 

% Dop_x = param_DoP .* param_a ./ sqrt( param_linear_x );
% Dop_y = param_DoP .* param_b ./ sqrt( stocks_s_polaried );

Dop_x = param_DoLP .* cos( param_Azimuth );
Dop_y = param_DoLP .* sin( param_Azimuth );

Dop_x_rect = param_DoLP_rect .* param_linear_x ./ sqrt( param_linear_x .* param_linear_x + param_linear_y .* param_linear_y + 0.0001) ;
Dop_y_rect = param_DoLP_rect .* param_linear_y ./ sqrt( param_linear_x .* param_linear_x + param_linear_y .* param_linear_y + 0.0001);

h_quiver1 = figure;
gridSize = 80;
dis_x_line = linspace(0,1,(gridSize+1));
dis_y_line = linspace(0,1,(gridSize+1));
dis_x_line = uint16 ( input_images_size(2) * dis_x_line(1:end-1) + 1);
dis_y_line = uint16 ( input_images_size(1) * dis_y_line(1:end-1) + 1);

[dis_x, dis_y] = meshgrid( dis_x_line, dis_y_line);

imshow( I3 ); 
%title('Image_LinearPolarAzimuthVector(red),corrected(green)'); 
hold on; 
quiver( dis_x, dis_y, Dop_x( dis_y_line, dis_x_line  ), -Dop_y( dis_y_line, dis_x_line), 'r' );
quiver( dis_x, dis_y, Dop_x_rect( dis_y_line, dis_x_line  ), -Dop_y_rect( dis_y_line, dis_x_line), 'g' );
saveas(gcf,strcat( ImagesOutputPath, 'Image_LinearPolarAzimuthVector_compare.emf'));

figure
imshow( I3 ); 
%title('Image_LinearPolarAzimuthVector(corrected)'); 
hold on; 
quiver( dis_x, dis_y, Dop_x_rect( dis_y_line, dis_x_line  ), -Dop_y_rect( dis_y_line, dis_x_line), 'g' );
saveas(gcf,strcat( ImagesOutputPath, 'Image_LinearPolarAzimuthVector_compare.emf'));



% 偏振图像
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



%% 分析项
% 给定区域（Region of interest），计算该区域的线偏振角度（均值），以及标准差
ROI_start = [ 360 460 ];
ROI_end   = [ 400 500 ];

% 偏振度一致性
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



% 线偏振一致性
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


% 应用 SPIE识别
param_Azimuth_1 =  acos(  sqrt(Ix_polaried ) ./ sqrt(Ix_polaried  + Iy_polaried + 0.001) )  ;
YCbCr(:,:,1) = uint8( 128*(param_DoP_rect) ) ;
YCbCr(:,:,2) = uint8( 128 - 40*( ( param_Azimuth_1 - min(min(param_Azimuth_1))) /( max(max(param_Azimuth_1)) - min(min(param_Azimuth_1))) ) ) ;
YCbCr(:,:,3) = uint8( 128 + 40*( ( param_Azimuth_1 - min(min(param_Azimuth_1))) /( max(max(param_Azimuth_1)) - min(min(param_Azimuth_1))) ) ) ;
rgb = ycbcr2rgb( YCbCr ) ;
figure
imshow(rgb);
%title('Application: recognition from mix polarization');
saveas(gcf,strcat( ImagesOutputPath, 'ApplicationRcgnitnSPIE.emf'));
imwrite( uint8(rgb), strcat( ImagesOutputPath, 'ApplicationRcgnitnSPIE.png') ); 
clear rgb
clear YCbCr




