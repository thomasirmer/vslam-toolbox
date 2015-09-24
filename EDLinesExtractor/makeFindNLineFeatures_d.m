disp('Compiling FindNLineFeatures.cpp');
mexOpenCV 'FindNLineFeatures.cpp' './LineMatching/EDLineDetector.cpp' './LineMatching/LineDescriptor.cpp'
