% Photodiode analysis
figure()
plot (time, data(:,3));
figure()
plot (time, data(:,2));

% beginning and end of frames
a = find (data(:,2)>4.99, 1,"first");
frames_on = time (find (data(:,2)>4.99, 1,"first"));
frames_off = time (find (data(:,2)>4.99, 1,"last"));
frames_t = frames_off - frames_on;

%gray screen and gratings on and off
first = max (data (1:a+2000,3));
t_first = time (find (data(1:a+2000,3)== first));
gray_on = t_first(end);

second = max (data (a+2001:a+12000,3));
t_second = time (find (data(1:a+11000,3)== second));
gratings_on= t_second(end);

third = max (data (a+12001:122000,3));
t_third = time (find (data(:,3)== third));
gratings_off = t_third (end);

gray_t = gratings_on - gray_on;
gratings_t = gratings_off - gratings_on;

gray_t_2 = gratings_on - frames_on;
gratings_t_2 = gratings_off - frames_on;

timing (1,1) = gray_t;
timing (2,1) = gratings_t;
timing (3,1) = gray_t_2;
timing (4,1) = gratings_t_2;
