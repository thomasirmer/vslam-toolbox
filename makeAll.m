path = strcat(pwd, '\CInterface\');
disp(['Going to ', path]);
cd(path);
make;
cd ..

path = strcat(pwd, '\EDLinesExtractor\');
disp(['Going to ', path]);
cd(path);
makeLineFeatures_d;
cd ..

path = strcat(pwd, '\FrameTransforms\');
disp(['Going to ', path]);
cd(path);
make;
cd ..

path = strcat(pwd, '\Math\');
disp(['Going to ', path]);
cd(path);
make;
cd ..

path = strcat(pwd, '\Points\');
disp(['Going to ', path]);
cd(path);
make;
cd ..