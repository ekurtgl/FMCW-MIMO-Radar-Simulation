function create_body_model(body_index)
    % body_index = 3x25
    
    % 0.001 m point intervals
    
    % kinect body index reference
    % 1 - spine base
    % 2 - spine mid
    % 3 - neck
    % 4 - head
    % 5 - left shoulder
    % 6 - left elbow
    % 7 - left wrist
    % 8 - left hand
    % 9 - right shoulder
    %10 - right elbow
    %11 - right wrist
    %12 - right hand
    %13 - left hip
    %14 - left knee
    %15 - left ankle
    %16 - left foot
    %17 - right hip
    %18 - right knee
    %19 - right ankle
    %20 - right foot
    %21 - spine shoulder
    %22 - left hand tip
    %23 - left thumb
    %24 - right hand tip
    %25 - right thumb
    
    head = body_index(:,4);
    neck = body_index(:,3);
    spine_top = body_index(:,21);
    lshoulder = body_index(:,5);
    rshoulder = body_index(:,9);
    base = body_index(:,1);
    spine_mid = body_index(:,2);
    lhip = body_index(:,13);
    rhip = body_index(:,17);
    lelbow = body_index(:,6);
    relbow = body_index(:,10);
    lwrist = body_index(:,7);
    rwrist = body_index(:,11);
    lhand = body_index(:,8);
    rhand = body_index(:,12);
    lknee = body_index(:,14);
    rknee = body_index(:,18);
    lankle = body_index(:,15);
    rankle = body_index(:,19);
    ltoe = body_index(:,16);
    rtoe = body_index(:,20);
    
    
    headlen = sqrt(sum((head(:,1)-neck(:,1)).^2));
    shoulderlen = sqrt(sum((neck(:,1)-lshoulder(:,1)).^2)); 
    torsolen = sqrt(sum((neck(:,1)-base(:,1)).^2));
    hiplen = sqrt(sum((lhip(:,1)-base(:,1)).^2));
    upperleglen =  sqrt(sum((lhip(:,1)-lknee(:,1)).^2));
    lowerleglen = sqrt(sum((lankle(:,1)-lknee(:,1)).^2));
    footlen = sqrt(sum((lankle(:,1)-ltoe(:,1)).^2));
    upperarmlen = sqrt(sum((lelbow(:,1)-lshoulder(:,1)).^2));
    lowerarmlen = sqrt(sum((lelbow(:,1)-lhand(:,1)).^2));


%     for i = 1:size(koordinatlar,2)
%         if mod(i,160)==0
i=1;
            clf
            hold on
            [x,y,z] = cylinder2P(neck(:,i)',head(:,i)',0.04,10);
            surf(x,y,z) % neck
            [x,y,z] = ellipsoid2P(base(:,i)',neck(:,i)',0.15,0.15,torsolen/2,20);
            surf(x,y,z) % torso
            [x,y,z] = ellipsoid2P(lankle(:,i)',lknee(:,i)',0.06,0.06,lowerleglen/2,10);
            surf(x,y,z) % left lower leg
            [x,y,z] = ellipsoid2P(rankle(:,i)',rknee(:,i)',0.06,0.06,lowerleglen/2,10);
            surf(x,y,z) % right lower leg
            [x,y,z] = ellipsoid2P(ltoe(:,i)',lankle(:,i)',0.05,0.05,footlen/2,10);
            surf(x,y,z) % left foot
            [x,y,z] = ellipsoid2P(rtoe(:,i)',rankle(:,i)',0.05,0.05,footlen/2,10);
            surf(x,y,z) % right foot
            [x,y,z] = ellipsoid2P(lknee(:,i)',lhip(:,i)',0.07,0.07,upperleglen/2,10);
            surf(x,y,z) % left upper leg
            [x,y,z] = ellipsoid2P(rknee(:,i)',rhip(:,i)',0.07,0.07,upperleglen/2,10);
            surf(x,y,z) % right upper leg
            [x,y,z] = ellipsoid2P(base(:,i)',lhip(:,i)',0.07,0.07,hiplen/2,10);
            surf(x,y,z) % left hip
            [x,y,z] = ellipsoid2P(base(:,i)',rhip(:,i)',0.07,0.07,hiplen/2,10);
            surf(x,y,z) % right hip
            [x,y,z] = ellipsoid2P(lelbow(:,i)',lshoulder(:,i)',0.06,0.06,upperarmlen/2,10);
            surf(x,y,z) % left upper arm
            [x,y,z] = ellipsoid2P(relbow(:,i)',rshoulder(:,i)',0.06,0.06,upperarmlen/2,10);
            surf(x,y,z) % right upper arm
            [x,y,z] = ellipsoid2P(neck(:,i)',lshoulder(:,i)',0.06,0.06, shoulderlen/2,10);
            surf(x,y,z) % left shoulder
            [x,y,z] = ellipsoid2P(neck(:,i)',rshoulder(:,i)',0.06,0.06, shoulderlen/2,10);
            surf(x,y,z) % right shoulder
            [x,y,z] = ellipsoid2P(lhand(:,i)',lelbow(:,i)',0.05,0.05,lowerarmlen/2,10);
            surf(x,y,z) % left lower arm
%             sqrt(sum((lelbow(i,1)-lhand(i,1)).^2))
            

            [x,y,z] = ellipsoid2P(rhand(:,i)',relbow(:,i)',0.05,0.05,lowerarmlen/2,10);
            surf(x,y,z) % right lower arm

            [x,y,z] = sphere(10);
            surf(x*0.02+ltoe(1,i),y*0.02+ltoe(2,i),z*0.01+ltoe(3,i)) % left toe

            [x,y,z] = sphere(10);
            surf(x*0.02+rtoe(1,i),y*0.02+rtoe(2,i),z*0.01+rtoe(3,i)) % right toe

            [x,y,z] = sphere(10);
            surf(x*0.05+lankle(1,i),y*0.05+lankle(2,i),z*0.05+lankle(3,i)) % lankle

            [x,y,z] = sphere(10);
            surf(x*0.05+rankle(1,i),y*0.05+rankle(2,i),z*0.05+rankle(3,i)) % rankle

            [x,y,z] = sphere(10);
            surf(x*0.05+lknee(1,i),y*0.05+lknee(2,i),z*0.05+lknee(3,i)) % lknee

            [x,y,z] = sphere(10);
            surf(x*0.05+rknee(1,i),y*0.05+rknee(2,i),z*0.05+rknee(3,i)) % rknee

            [x,y,z] = sphere(10);
            surf(x*0.1+head(1,i),y*0.1+head(2,i),z*headlen/2+head(3,i)) % head

            [x,y,z] = sphere(10);
            surf(x*0.05+lhip(1,i),y*0.05+lhip(2,i),z*0.05+lhip(3,i)) % lhip

            [x,y,z] = sphere(10);
            surf(x*0.05+rhip(1,i),y*0.05+rhip(2,i),z*0.05+rhip(3,i)) % rhip

            [x,y,z] = sphere(10); % left shoulder
            surf(x*0.05+lshoulder(1,i),y*0.05+lshoulder(2,i),z*0.05+lshoulder(3,i))

            [x,y,z] = sphere(10); % right shoulder
            surf(x*0.05+rshoulder(1,i),y*0.05+rshoulder(2,i),z*0.05+rshoulder(3,i))

            [x,y,z] = sphere(10);
            surf(x*0.05+lelbow(1,i),y*0.05+lelbow(2,i),z*0.05+lelbow(3,i)) % lelbow

            [x,y,z] = sphere(10);
            surf(x*0.05+relbow(1,i),y*0.05+relbow(2,i),z*0.05+relbow(3,i)); % relbow

            [x,y,z] = sphere(10);
            surf(x*0.05+lhand(1,i),y*0.05+lhand(2,i),z*0.05+lhand(3,i)); % lhand

            [x,y,z] = sphere(10);
            surf(x*0.05+rhand(1,i),y*0.05+rhand(2,i),z*0.05+rhand(3,i)); % rhand
%             xlim([-1.2 1.2]);ylim([-1.2 1]);zlim([-2 3]);
            view([180,-90])
            light
            lighting none
            shading interp
            axis equal

            xlabel('x')
            ylabel('y')
            zlabel('z')
            grid off
            set(gcf,'Color',[1 1 1])
            view([180,-90])
%             pause(0.03);
            drawnow
            
%         end
%     end

end

