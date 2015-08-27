disp('Compiling FindNLineFeatures.cpp');
mexOpenCV -g 'FindNLineFeatures.cpp' 'EDLineDetector.cpp' 'LineDescriptor.cpp'
