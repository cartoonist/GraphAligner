#ifndef Aligner_h
#define Aligner_h

#include <string>
#include <vector>
#include "AlignmentGraph.h"
#include "vg.pb.h"
#include "AlignmentSelection.h"

struct AlignerParams
{
	std::string graphFile;
	std::vector<std::string> fastqFiles;
	size_t numThreads;
	size_t initialBandwidth;
	size_t rampBandwidth;
	int dynamicRowStart;
	size_t maxCellsPerSlice;
	std::vector<std::string> seedFiles;
	std::string outputGAMFile;
	std::string outputJSONFile;
	std::string outputGAFFile;
	std::string outputCorrectedFile;
	std::string outputCorrectedClippedFile;
	bool verboseMode;
	bool tryAllSeeds;
	bool highMemory;
	size_t psiLength;
	size_t psiChunkSize;
	size_t psiDistance;
	size_t psiPathCount;
	size_t psiContext;
	size_t psiStep;
	size_t psiGoccThreshold;
	size_t mxmLength;
	size_t mumCount;
	size_t memCount;
	std::string seederCachePrefix;
	AlignmentSelection::SelectionMethod alignmentSelectionMethod;
	double selectionECutoff;
	bool forceGlobal;
	bool compressCorrected;
	bool compressClipped;
	bool preciseClipping;
	size_t minimizerLength;
	size_t minimizerWindowSize;
	double minimizerSeedDensity;
	size_t seedClusterMinSize;
	double minimizerDiscardMostNumerousFraction;
	double seedExtendDensity;
	bool nondeterministicOptimizations;
	bool optimalDijkstra;
};

void alignReads(AlignerParams params);
void replaceDigraphNodeIdsWithOriginalNodeIds(vg::Alignment& alignment, const AlignmentGraph& graph);

#endif
