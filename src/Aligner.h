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
	size_t readMinInsertSize;
	size_t readMaxInsertSize;
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
	bool interleaved;
};

/* tags */
struct SingleEnd {};
struct PairedEnd {};

/* queue entry wrapper class */
template<typename TSpec, typename TValue>
struct Entry;

/* queue entry wrapper class specialised for SingleEnd reads */
template<typename TValue>
struct Entry<SingleEnd, TValue> {
	using end_type = std::shared_ptr<TValue>;
	using type = end_type;

	static inline end_type&
	getNullEnd()
	{
		static end_type null_end = nullptr;
		return null_end;
	}
};

template<typename TValue>
inline std::shared_ptr<TValue>&
getOneEnd(std::shared_ptr<TValue>& entry)
{
	return entry;
}

template<typename TValue>
inline std::shared_ptr<TValue>&
getOtherEnd(std::shared_ptr<TValue>& entry)
{
	return Entry<SingleEnd, TValue>::getNullEnd();
}

template<typename TValue>
inline bool
isEmpty(std::shared_ptr<TValue>& entry)
{
	return entry == nullptr;
}

template<typename TValue>
inline void
clear(std::shared_ptr<TValue>& entry)
{
	entry = nullptr;
}

template<typename TValue>
constexpr inline bool
isSingle(std::shared_ptr<TValue>& entry)
{
	return true;
}

/* queue entry wrapper class specialised for PairedEnd reads */
template<typename TValue>
struct Entry<PairedEnd, TValue> {
	using end_type = std::shared_ptr<TValue>;
	using type = std::pair<end_type, end_type>;
};

template<typename TValue>
inline std::shared_ptr<TValue>&
getOneEnd(std::pair<std::shared_ptr<TValue>, std::shared_ptr<TValue>>& entry)
{
	return entry.first;
}

template<typename TValue>
inline std::shared_ptr<TValue>&
getOtherEnd(std::pair<std::shared_ptr<TValue>, std::shared_ptr<TValue>>& entry)
{
	return entry.second;
}

template<typename TValue>
inline bool
isEmpty(std::pair<std::shared_ptr<TValue>, std::shared_ptr<TValue>>& entry)
{
	return entry.first == nullptr;
}

template<typename TValue>
inline void
clear(std::pair<std::shared_ptr<TValue>, std::shared_ptr<TValue>>& entry)
{
	entry.first = nullptr;
	entry.second = nullptr;
}

template<typename TValue>
constexpr inline bool
isSingle(std::pair<std::shared_ptr<TValue>, std::shared_ptr<TValue>>& entry)
{
	return false;
}

template<typename TSpec = SingleEnd>
void alignReads(AlignerParams params, TSpec=TSpec{});

void replaceDigraphNodeIdsWithOriginalNodeIds(vg::Alignment& alignment, const AlignmentGraph& graph);

#endif
