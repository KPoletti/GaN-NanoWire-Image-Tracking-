
close all; clear all;clc
%{
NAME: Nanowire Tracker V4.2.1

PURPOSE: Track the motion of a nanowire or similarly
rectangular/cylindrical object for 4 degrees of freedom. Degrees of freedom are "x"-position, "y"-position, "theta"-position, and "phi"-position. 
Determine average dimensional characteristics of object. Calculate windowed and averaged power spectra for four degrees of
freedom. 

INPUTS:
Suitable video file to track in .avi format.
 
PARAMETERS:
none

OUPUTS:

PROCESS:


CREATED: Max Chen, and Keith Poletti University of Colorado Boulder, June 14,2019

%}
%% Do you want to plot?
% if yes set to 1 (will increase run time)
% if no  set to 0
DoPlot = 0;

%%
%filename='GaNPow1_mov3_20fps_0_54V.avi';
filename='GaNNW5_MOV1_30fps_0_6V.avi';

%Determine total number of frames in the video file
videoFReader = vision.VideoFileReader(filename);

%Determine total time duration of video and frame rate of video
video=VideoReader(filename);
duration=video.Duration;
framerate=video.FrameRate;
totalframes=round(duration*framerate);

secondsperframe=duration/totalframes;
%{
TRACK PARAMETERS: Track parameter adjustment and optimization is critical
to achieving an accruate track. Assuming adequate illumination and by
extension and acceptable number of "black" frames, optimal track parameters
should yield in rougly a 10% error rate. Adjust track parameters while
enabling the optional TRACK IMAGE DISPLAY code to acheive a qualitatively
acceptable track, while also acheiving the quantitaive 10% error goal. 
%}
%% Parameters
%only need to change th it must be a whole number. Increasing it will make
%the box around the nanowire larger. Decreasing decreases the box size
mperpixel=(8.944*10^(-7));

%%
%Set frame length of segment based on total number of frames and desired number of
%windows/segments
numberofwindows=(10);

%Round total frame count to nearest 100s
roundedtotalframe = round(totalframes,-2);
%Determine number of frames out of rounded total per 100th percent tracked
numframetrackedper1percent = roundedtotalframe./100;

segmentlength=round(totalframes/numberofwindows);

%Preallocate data matrix for entire data set



errorframes = {};
outnew=[];
error=0;
saveang=[];
frame=0;
while ~isDone(videoFReader)
try
     frame=frame+1;  
     apple=videoFReader();
    im=double(rgb2gray(apple));
    [nr,nc]=size(im);
    if frame==1 
       [~,xloc]=max(max(im)); 
    else
        xloc=outnew(end,2);
    end
    %creates bins and counts for a histogram of each pixels brightness
    [counts,edge]=imhist(im);
    [~,idmaxcount]=max(counts);
    thresh=otsuthresh(counts);
    BW=imbinarize(im,thresh);
%    im(BW==0)=0;
    nanowiredim =regionprops(BW,'Centroid','MajorAxisLength','MinorAxisLength','Orientation');
    if length(nanowiredim)>1
        centx=extractfield(nanowiredim,'Centroid');
        centx=centx(1:2:2*length(nanowiredim));
        [~,s]=min(abs(centx-xloc));
    else
        s=1;
    end
    pixellength=nanowiredim(s).MajorAxisLength;
    [nr,~]=size(im);
%     if pixellength>=.5*nr
%         error("No particle detected")
%     end
    angularpos=deg2rad(nanowiredim(s).Orientation);
    
    im=im*256;
    %finds the max brightness of each column
    bigbright=max(im);
    [~,xbest]=max(bigbright);
    %finds the max brightness of each row
    bigbrightr=max(im,[],2);
    [~,ybest]=max(bigbrightr);
    bbc=find(bigbright>=256*thresh);
    bbr=find(bigbrightr>=256*thresh);
    if length(bbc)<=1
        bbc=find(bigbright>=256*thresh-1);
    elseif length(bbr)<=1
        bbr=find(bigbrightr>=256*thresh-1);
    end
    totbright=sum(sum(im(bbr,bbc)));
    cobx=sum(bbc.*sum(im(bbr,bbc)))/totbright;
    coby=sum(bbr.*sum(im(bbr,bbc),2))/totbright;
    
    middle=nanowiredim(s).Centroid;
    actuallinex=[middle(1)+pixellength/2*cos(angularpos),middle(1),middle(1)-pixellength/2*cos(angularpos)];
    actualliney=[middle(2)-pixellength/2*sin(angularpos),middle(2),middle(2)+pixellength/2*sin(angularpos)];

    angularpos=-1*sign(nanowiredim(s).Orientation)*(pi/2-abs(deg2rad(nanowiredim(s).Orientation)));
    if pixellength>=0.5*76
        t = linspace(0,2*pi,50);
        a = nanowiredim(s).MajorAxisLength/2;
        b = nanowiredim(s).MinorAxisLength/2;
        Xc =nanowiredim(s).Centroid(1);
        Yc = nanowiredim(s).Centroid(2);
        phi = deg2rad(-nanowiredim(s).Orientation);
        x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
        y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
        imshow(uint8(im),'InitialMagnification',525)
        drawnow;
        hold on;
        % plot(coblx(:,1),coblx(:,2),'LineWidth',0.25,'Color','red')
        % plot(cobly(:,1),cobly(:,2),'LineWidth',0.25,'Color','magenta')
         plot(middle(1),middle(2),'r.','MarkerSize',10,'Color','blue')
        plot(actuallinex,actualliney,'LineWidth',0.25,'Color','cyan')
        plot(x,y,'g','Linewidth',0.5)

        pause(.1);
        hold off

        error('No nanowire found')
    end
    outnew=[outnew;frame,middle,pixellength,angularpos,cobx,coby,xbest,ybest];
    if rem(frame./numframetrackedper1percent,1)==0
        disp('Track Progress %:')
        disp(frame./numframetrackedper1percent)
        if frame./numframetrackedper1percent==1
        else
            timeper1percent = toc;
            trackpercentremain = 100-(frame./numframetrackedper1percent);
            tracktimeremain = datestr(((timeper1percent*trackpercentremain)./(24*3600)),'HH:MM:SS');
            disp('Time Remaining:')
            disp(tracktimeremain)
        end
        tic
    end
%
if DoPlot ==1 
    t = linspace(0,2*pi,50);
    a = nanowiredim(s).MajorAxisLength/2;
    b = nanowiredim(s).MinorAxisLength/2;
    Xc =nanowiredim(s).Centroid(1);
    Yc = nanowiredim(s).Centroid(2);
    phi = deg2rad(-nanowiredim(s).Orientation);
    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    imshow(apple,'InitialMagnification',525)
    drawnow;
    hold on;
    % plot(coblx(:,1),coblx(:,2),'LineWidth',0.25,'Color','red')
    % plot(cobly(:,1),cobly(:,2),'LineWidth',0.25,'Color','magenta')
     plot(middle(1),middle(2),'r.','MarkerSize',10,'Color','blue')
    plot(actuallinex,actualliney,'LineWidth',0.75,'Color','cyan')
    plot(x,y,'g','Linewidth',0.75)

    pause(secondsperframe);
    hold off
end


catch
%Catch any error and populate current frame row with track data from last
%previously successful frame track. Index error counter and report what
%frame had error.
error=error+1;

% if  logical(frsttherror)==1
%     display('Error: nothing above first threshold');
% elseif frstthbndrynpks==1
%     display('Error: No peaks within boundary and above first threshold');
% end

errorframes = [errorframes frame];

% totalpos_lst(frame,1) = frame;
%     totalpos_lst(frame,2) = totalpos_lst((frame-1),2);
%     totalpos_lst(frame,3) = totalpos_lst((frame-1),3);
%     totalpos_lst(frame,4) = totalpos_lst((frame-1),4);
%     totalpos_lst(frame,5) = totalpos_lst((frame-1),5); 
[outr,~]=size(outnew);
% if outr==0
%     outnew=[outnew;frame,middle,pixellength,angularpos,cobx,coby,cobx,coby];
% else
%     outnew=[outnew;frame,middle,pixellength,angularpos,cobx,coby,outnew(outr,8:9)];
% end

end

  
end
errorframes = transpose(errorframes);
%%
totalframe=outnew(:,1);
totalxpixel = outnew(:,2);
totalypixel = outnew(:,3);
totallengthpixel = outnew(:,4);
totalangularpos = outnew(:,5);
totalxpixel2 = outnew(:,6);
totalypixel2 = outnew(:,7);
totalxpixelb = outnew(:,8);
totalypixelb= outnew(:,9);
%%
%Convert frame number to time and pixel location to spatial location.
%Convert nanowire pixel length to meter length.
totaltime=secondsperframe*(totalframe-1);
totalxpos=totalxpixel*mperpixel;
totalypos=totalypixel*mperpixel;
totalxpos2=totalxpixel2*mperpixel;
totalypos2=totalypixel2*mperpixel;
totalxposb=totalxpixelb*mperpixel;
totalyposb=totalypixelb*mperpixel;



totalmeterlength = totallengthpixel*mperpixel;
MeanNWLength = mean(totalmeterlength);
%HighMeanNWLength=mean(totalmeterlength(totalmeterlength>=MeanNWLength));
HighMeanNWLength=mean(findpeaks(totalmeterlength));

[~,maxlengthloc]=findpeaks(totalmeterlength(HighMeanNWLength>=totalmeterlength));
azimuth= acos(totalmeterlength(HighMeanNWLength>=totalmeterlength)/HighMeanNWLength);
inverseazimuth=max(azimuth)-azimuth;
[~,flatloc]=findpeaks(inverseazimuth);
for i=1:length(maxlengthloc)-1
    k=maxlengthloc(i);
    j=maxlengthloc(i+1);
   azimuth(k:j)=(-1)^i*azimuth(k:j); 
end

totaltimepadded=[ones(length(totaltime(HighMeanNWLength>=totalmeterlength)),1),totaltime(HighMeanNWLength>=totalmeterlength)];
mazbaz=totaltimepadded\azimuth;
azilinearfit=totaltimepadded*mazbaz;
azimuthcorrected=azimuth-azilinearfit;

%%
%Find linear fits for the x and y position data. Find linear fits for the
%length and angular position data.
totaltimepadded=[ones(length(totaltime),1),totaltime];
mxbx=totaltimepadded\totalxpos;
myby=totaltimepadded\totalypos;
xlinearfit=totaltimepadded*mxbx;
ylinearfit=totaltimepadded*myby;

mxbx=totaltimepadded\totalxpos2;
myby=totaltimepadded\totalypos2;
xlinearfit2=totaltimepadded*mxbx;
ylinearfit2=totaltimepadded*myby;

mxbx=totaltimepadded\totalxposb;
myby=totaltimepadded\totalyposb;
xlinearfitb=totaltimepadded*mxbx;
ylinearfitb=totaltimepadded*myby;

mlbl=totaltimepadded\totalmeterlength;
maba=totaltimepadded\totalangularpos;
lengthlinearfit=totaltimepadded*mlbl;
angularlinearfit=totaltimepadded*maba;

%Create x and y position lists that are corrected for linear drift. 
%Create length and angular position lists that are corrected for linear drift. 
xposlinearcorrected=totalxpos-xlinearfit;
yposlinearcorrected=totalypos-ylinearfit;

xposlinearcorrected2=totalxpos2-xlinearfit2;
yposlinearcorrected2=totalypos2-ylinearfit2;

xposlinearcorrectedb=totalxposb-xlinearfitb;
yposlinearcorrectedb=totalyposb-ylinearfitb;

lengthlinearcorrected=totalmeterlength-lengthlinearfit;
angularlinearcorrected=totalangularpos-angularlinearfit;




%Take Welch averaged PSDs for the linear corrected x and y positions.
%Take Welch averaged PSDs for the linear corrected length and angular
%position variations. 
[pwelchxcorrected,f]=pwelch(xposlinearcorrected,segmentlength,round(segmentlength/2),segmentlength,framerate);
[pwelchycorrected,f]=pwelch(yposlinearcorrected,segmentlength,round(segmentlength/2),segmentlength,framerate);
[pwelchlengthcorrected,f]=pwelch(lengthlinearcorrected,segmentlength,round(segmentlength/2),segmentlength,framerate);
[pwelchanglecorrected,f]=pwelch(angularlinearcorrected,segmentlength,round(segmentlength/2),segmentlength,framerate);

[pwelchxcorrected2,f]=pwelch(xposlinearcorrected2,segmentlength,round(segmentlength/2),segmentlength,framerate);
[pwelchycorrected2,f]=pwelch(yposlinearcorrected2,segmentlength,round(segmentlength/2),segmentlength,framerate);

[pwelchxcorrectedb,f]=pwelch(xposlinearcorrectedb,segmentlength,round(segmentlength/2),segmentlength,framerate);
[pwelchycorrectedb,f]=pwelch(yposlinearcorrectedb,segmentlength,round(segmentlength/2),segmentlength,framerate);

figure('Name','Linear Corrected: Welch PSD')
hold on
subplot(2,2,1)
loglog(f,pwelchxcorrected,f,pwelchxcorrected2,f,pwelchxcorrectedb)
legend('Built in Function','Centroid first method','Brightest X Location')


% ylim([10^-16 10^-11])
% xlim([10^-2 framerate/2])
xlabel('Frequency (Hz)')
ylabel('Power (m^2/Hz)')
title('Linear Corrected Signal Averaged Power in X')

subplot(2,2,2)
loglog(f,pwelchycorrected,f,pwelchycorrected2,f,pwelchycorrectedb)
legend('Built in Function','Centroid first method','Brightest Y Location')
% ylim([10^-16 10^-11])
% xlim([5*10^-2 framerate/2])
xlabel('Frequency (Hz)')
ylabel('Power (m^2/Hz)')
title('Linear Corrected Signal Averaged Power in Y')


subplot(2,2,3)
loglog(f,pwelchlengthcorrected)
% ylim([10^-15 10^-9])
% xlim([10^-2 framerate/2])
xlabel('Frequency (Hz)')
ylabel('Power (m^2/Hz)')
title('Linear Corrected Signal Averaged Power in length')

subplot(2,2,4)
loglog(f,pwelchanglecorrected)
% ylim([10^-15 10^-9])
% xlim([10^-2 framerate/2])
xlabel('Frequency (Hz)')
ylabel('Power (Rad^2/Hz)')
title('Linear Corrected Signal Averaged Power in angular position')
figure('Name','x, y, length, and angular Position')

subplot(2,2,1)
%scatter(totaltime,lengthlinearcorrected)
plot(totaltime,totalmeterlength)
hold on
plot(totaltime, HighMeanNWLength*ones(length(totaltime),1))
ax=gca;
% xlim([0 duration])
% ylim([26.5*10^-6 149.5*10^-6])
xlabel('Time (s)')
ylabel('Length (m)')
title('Length variation Linear Corrected')
ax.YAxis.Exponent=-6;
hold on

subplot(2,2,2)
plot(totaltime,angularlinearcorrected)
ax=gca;
% xlim([0 duration])
% ylim([1.5*10^-6 127.5*10^-6])
xlabel('Time (s)')
ylabel('Angle (rad)')
title('Angular Position Linear Corrected')
% ax.YAxis.Exponent=-6;
%%
figure;
subplot(3,2,1)
plot(totaltime,xposlinearcorrected)
xlim([0 duration])
xlabel('Time (s)')
ylabel('X pos (m)')
title('X Position Linear Corrected')
hold on;
subplot(3,2,3)
plot(totaltime,xposlinearcorrected2)
xlim([0 duration])
xlabel('Time (s)')
ylabel('X pos (m)')
title('X Position of Brightness Centroid Linear Corrected')
subplot(3,2,5)
plot(totaltime,xposlinearcorrectedb)
ax=gca;
xlim([0 duration])
ylim([-15*10^-6 15*10^-6])
xlabel('Time (s)')
ylabel('X pos (m)')
title('X Position of Brightest Point Linear Corrected')
%legend('Built in Function','Centroid first method','Brightest Y Location')
ax.YAxis.Exponent=-6;


% subplot(2,2,4)
% hold on;
% plot(totaltime,yposlinearcorrected)
% plot(totaltime,yposlinearcorrected2)
% plot(totaltime,yposlinearcorrectedb)
% ax=gca;
% xlim([0 duration])
% ylim([-15*10^-6 15*10^-6])
% xlabel('Time (s)')
% ylabel('Y pos (m)')
% title('Y Position Linear Corrected')
% legend('Built in Function','Centroid first method','Brightest Y Location')
% ax.YAxis.Exponent=-6;
% hold on
subplot(3,2,2)
plot(totaltime,yposlinearcorrected)
xlim([0 duration])
xlabel('Time (s)')
ylabel('Y pos (m)')
title('Y Position Linear Corrected')
hold on;
subplot(3,2,4)
plot(totaltime,yposlinearcorrected2)
xlim([0 duration])
xlabel('Time (s)')
ylabel('Y pos (m)')
title('Y Position of Brightness Centroid Linear Corrected')
subplot(3,2,6)
plot(totaltime,yposlinearcorrectedb)
ax=gca;
xlim([0 duration])
ylim([-15*10^-6 15*10^-6])
xlabel('Time (s)')
ylabel('Y pos (m)')
title('Y Position of Brightest Point Linear Corrected')
%legend('Built in Function','Centroid first method','Brightest Y Location')
ax.YAxis.Exponent=-6;

%DATA EXPORT: This is optional code to create and populate an Excel
%spreadsheet with the track data. Position data is tabulated, then PSD data
%is tabulated. Average dimensional characteristics are tabulated. The new
%file name must be specified by "TotalDatafilename". Destination folder is
%same as that containing the track code executed.

%{
PosData = table(totaltime,xposlinearcorrected,yposlinearcorrected,totalmeterlength,angularlinearcorrected);
PSDData = table(f,pwelchxcorrected,pwelchycorrected,pwelchlengthcorrected,pwelchanglecorrected);
CharacteristicData = table(MeanNWLength);
TotalDatafilename = strrep(filename,'avi','xlsx');
writetable(PosData,TotalDatafilename,'Sheet',1,'Range','A1')
writetable(PSDData,TotalDatafilename,'Sheet',1,'Range','F1')
writetable(CharacteristicData,TotalDatafilename, 'Sheet',2,'Range','A1') 
%}
