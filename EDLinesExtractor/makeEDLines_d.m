disp('Compiling EDLinesExtractor.cpp');
mexOpenCV -g 'EDLinesExtractor.cpp' 'EDLineDetector.cpp' 'LineDescriptor.cpp' 'PairwiseLineMatching.cpp'