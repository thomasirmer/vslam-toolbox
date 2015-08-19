disp('Compiling FindLineFeatures.cpp');
mexOpenCV -g 'FindNLineFeatures.cpp' 'EDLineDetector.cpp' 'LineDescriptor.cpp'
mexOpenCV -g 'MatchLines.cpp'