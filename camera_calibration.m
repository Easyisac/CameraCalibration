%%Camera calibration
% Zhang's procedure implementation for camera calibration to estimate camera parameters.
% To check correctness of the estimation it's been computed the
% reprojection error and some solid have been superimposed over the
% calibration plane.

%% Importing and processing images
%loading the images and then using the built.in function detect all the
%points for later use in estimating homographies.
%I stores the image and XYpyxel stores the coordinates in pixels of the
%calibration points.

iimage=[3 4 7 10 13 16]; %indices of the images to be processed
  
  for ii=1:length(iimage)
      imageFileName = fullfile('images',['image' num2str(iimage(ii)) '.tif']);      
     imageData(ii).I = imread(imageFileName);
     imageData(ii).XYpixel = detectCheckerboardPoints(imageData(ii).I);
  end


%% Establishing correspondences
%Knowing the square size we can estabilish a correspondence between the
%pixels and the real coordinates tahat will be stored in XYmm.
squaresize = 30;
for ii=1:length(iimage)
    XYpixel=imageData(ii).XYpixel;
    clear Xmm Ymm
    for jj=1:length(XYpixel)
        [row,col]=ind2sub([12,13],jj);
        Xmm=(col-1)*squaresize;
        Ymm=(row-1)*squaresize;
        imageData(ii).XYmm(jj,:)=[Xmm Ymm];
    end
    pause(1)
end

%% Estimating homographies
%Using the direct inear transform method it's estimated the homography for
%each image that it's stored in the field H
for ii= 1:length(iimage)
    XYpixel=imageData(ii).XYpixel;
    XYmm=imageData(ii).XYmm;
    A=[];
    b=[];
    for jj=1:length(XYpixel)

        Xpixel=XYpixel(jj,1);
        Ypixel=XYpixel(jj,2);
        Xmm=XYmm(jj,1);
        Ymm=XYmm(jj,2);

        m=[Xmm; Ymm; 1];
        O=[0;0;0];
        A=[A; m' O' -Xpixel*m'; O' m' -Ypixel*m'];
        b=[b;0;0];
    end
    [U,S,V]=svd(A);
    h=V(:,end);

    imageData(ii).H=reshape(h,[3 3])';
end

%% 1) calibrate using Zhang procedure (find the intrinsic parameters $K$ and, 
% for each image, the pair of $R,t$ (extrinsic) 

%A linear system Vb = 0 is formed V has shape 2nx6 where n is the number of
%homographies in input. from the 6 values of the vector b the symmetric
%matrix B can be computed to solve the system we can use the singular
%vector decomposition as is known that the solution is the right singular
%vector of V

V = [];
c = [1 1; 1 2; 2 2;];
for ii=1:length(iimage)
    H = imageData(ii).H;
    v = zeros(6,3);
    for jj = 1:length(c)
        i = c(jj, 1);
        j = c(jj, 2);
        v(:,jj) = [H(1,i)*H(1,j);...
            H(1,i)*H(2,j)+H(2,i)*H(1,j);...
            H(2,i)*H(2,j);...
            H(3,i)*H(1,j)+H(1,i)*H(3,j);...
            H(3,i)*H(2,j)+H(2,i)*H(3,j);...
            H(3,i)*H(3,j)];
    end
    V=[V; v(:,2)'; (v(:,1)-v(:,3))'];
end

% Since cholesky can be applied only on positive definite matrices if B is
% not positive definite we invert it. then after the factorization K is
% divided by his last value since it's a scale factor of K

[U,Z,S] = svd(V);
b6 = S(:,end);
B = [b6(1), b6(2), b6(4); b6(2) b6(3) b6(5); b6(4) b6(5) b6(6)];
if b6(1) < 0 || b6(3) < 0 || b6(6) < 0
    B = -B;
end


L = chol(B);
K = inv(L);
K = K./K(3,3);

%After we've found the intrinsic K we can find the intrinsic values R and t
%than we can compute the projection matrix P

for ii=1:length(iimage)
    lambda = 1/norm(inv(K)*imageData(ii).H(:,1));
    r1 = lambda*inv(K)*imageData(ii).H(:,1);
    r2 = lambda*inv(K)*imageData(ii).H(:,2);
    t = lambda*inv(K)*imageData(ii).H(:,3);
    r3 = cross(r1, r2);
    R = [r1 r2 r3];
    imageData(ii).P = K*[R t];
end

%% 2) choose one of the calibration images and compute the total reprojection error

% Reprojection error is computed using teh square distance between each
% couple of points generated by the built-in detection function and the
% computed projection matrix, we assume that detectCheckerboardPoints()
% returns the actual points.

points = [];
for jj=0:12
    for ii=0:11    
        p = [jj*squaresize; ii*squaresize; 0; 1];
        points = [points p];
    end
end
total = 0;
image = imageData(1);
for ii=1:length(image.XYpixel)
    err = ((image.P(1,:) * points(:,ii))/(image.P(3,:) * points(:,ii)) - image.XYpixel(ii,1))^2 + ...
        ((image.P(2,:) * points(:,ii))/(image.P(3,:) * points(:,ii)) - image.XYpixel(ii,2))^2;
    total = total + err;
end
homPixelPoints = image.P*points;
pixelPoints = [homPixelPoints(1,:)./homPixelPoints(3,:); homPixelPoints(2,:)./ homPixelPoints(3,:)];

figure
imshow(image.I,'InitialMagnification',300)
hold on

for jj= 1:size(image.XYpixel, 1)
    x = image.XYpixel(jj,1);
    y = image.XYpixel(jj,2);
    plot(x,y, 'or')
    x = pixelPoints(1,jj);
    y = pixelPoints(2,jj);
    plot(x,y, '+b')
end
pause(1)
total_reprojection_error = total

%% 3) superimpose a cilinder (or a solid of your choice) to each image

%To check the uality of the estimation we superimpose a cylinder on the
%image, the object will appear placed dikrectly on the calibration plane.
%The cylinder is represented by a red circle placed at height 0 a blue
%circle placed at height 120 and the vertical face is represented in green.
%The 3D object is first constructed in real coordinates that are
%transformed in homogeneous coordinates to be projected, after that, to plot
%the values the points will ber recomputed as inhomogenous coordinates for
%the 2D image.
for ii=1:length(iimage)
    radius = 60;
    center = [120 120];
    height = 120;
    np = 100;
    num = np*2+1;
    figure
    imshow(imageData(ii).I,'InitialMagnification',200)
    hold on
    th = 0:pi/np:2*pi;
    xpoints = radius * cos(th) + center(1);
    ypoints = radius * sin(th) + center(2);
    points = [xpoints xpoints; ypoints ypoints; zeros(1, length(xpoints)) height+zeros(1, length(xpoints)); ones(1, length(xpoints)*2)];
    proj_hom = imageData(ii).P*points;
    proj = [proj_hom(1,:)./proj_hom(3,:); proj_hom(2,:)./proj_hom(3,:)];
    f = fill(proj(1,1:num),proj(2,1:num), 'r');
    alpha(f, .5);
    f = fill(proj(1,num+1:end),proj(2,num+1:end), 'b');
    alpha(f, .5);
    for jj=1:num-1
        f = fill([proj(1,jj) proj(1, jj+num) proj(1,jj+num+1) proj(1, jj+1)],...
            [proj(2,jj) proj(2, jj+num) proj(2,jj+num+1) proj(2, jj+1)], 'g', 'LineStyle','none');
        alpha(f, .5);
    end
    pause(1)
end