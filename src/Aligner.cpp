#include <sstream>
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <utility>
#include <thread>
#include <csignal>
#include <concurrentqueue.h> //https://github.com/cameron314/concurrentqueue
#include <google/protobuf/util/json_util.h>
#include <psi/seed_finder.hpp>
#include <gum/io_utils.hpp>
#include "Aligner.h"
#include "CommonUtils.h"
#include "stream.hpp"
#include "fastqloader.h"
#include "BigraphToDigraph.h"
#include "ThreadReadAssertion.h"
#include "GraphAlignerWrapper.h"
#include "MummerSeeder.h"
#include "ReadCorrection.h"
#include "MinimizerSeeder.h"
#include "AlignmentSelection.h"

namespace psi {
#ifdef PSI_STATS
	using seeder_type = SeedFinder<WithStats>;
#else
	using seeder_type = SeedFinder<NoStats>;
#endif
	using graph_type = seeder_type::graph_type;
}  /* --- end of namespace psi --- */

struct Seeder
{
	enum Mode
	{
		None, File, Mum, Mem, Minimizer, Psi
	};
	Mode mode;
	size_t mumCount;
	size_t memCount;
	size_t mxmLength;
	size_t minimizerLength;
	size_t minimizerWindowSize;
	double minimizerSeedDensity;
	const MummerSeeder* mummerSeeder;
	const MinimizerSeeder* minimizerSeeder;
	const psi::seeder_type* psiSeeder;
	const std::unordered_map<std::string, std::vector<SeedHit>>* fileSeeds;
	Seeder(const AlignerParams& params, const std::unordered_map<std::string, std::vector<SeedHit>>* fileSeeds, const MummerSeeder* mummerSeeder, const MinimizerSeeder* minimizerSeeder, const psi::seeder_type* psiSeeder) :
		mumCount(params.mumCount),
		memCount(params.memCount),
		mxmLength(params.mxmLength),
		minimizerLength(params.minimizerLength),
		minimizerWindowSize(params.minimizerWindowSize),
		minimizerSeedDensity(params.minimizerSeedDensity),
		mummerSeeder(mummerSeeder),
		minimizerSeeder(minimizerSeeder),
		psiSeeder(psiSeeder),
		fileSeeds(fileSeeds)
	{
		mode = Mode::None;
		if (fileSeeds != nullptr)
		{
			assert(psiSeeder == nullptr);
			assert(minimizerSeeder == nullptr);
			assert(mummerSeeder == nullptr);
			assert(mumCount == 0);
			assert(memCount == 0);
			assert(minimizerSeedDensity == 0);
			mode = Mode::File;
		}
		if (psiSeeder != nullptr)
		{
			assert(minimizerSeeder == nullptr);
			assert(fileSeeds == nullptr);
			assert(mummerSeeder == nullptr);
			assert(mumCount == 0);
			assert(memCount == 0);
			assert(minimizerSeedDensity == 0);
			mode = Mode::Psi;
		}
		if (minimizerSeeder != nullptr)
		{
			assert(psiSeeder == nullptr);
			assert(mummerSeeder == nullptr);
			assert(mumCount == 0);
			assert(memCount == 0);
			assert(minimizerSeedDensity != 0);
			mode = Mode::Minimizer;
		}
		if (mummerSeeder != nullptr)
		{
			assert(minimizerSeeder == nullptr);
			assert(psiSeeder == nullptr);
			assert(fileSeeds == nullptr);
			assert(mumCount != 0 || memCount != 0);
			assert(minimizerSeedDensity == 0);
			if (mumCount != 0)
			{
				mode = Mode::Mum;
				assert(memCount == 0);
			}
			if (memCount != 0)
			{
				mode = Mode::Mem;
				assert(mumCount == 0);
			}
		}
	}
	/**
	 *  @brief  Verify the distance between two alignments of a paired-end read.
	 *
	 *  @params aln1 the first alignment
	 *  @params aln2 the second alignment
	 *  @return `true` if they meet the distance constraints; i.e. there is a path of the
	 *  length in the range [m, M], where m and M is minimum and maximum insert size
	 *  respectively. Otherwise, it returns `false`.
	 *
	 *  NOTE: This function assumes that the input reads are forward-reversed.
	 */
	bool verifyDistance(AlignmentResult::AlignmentItem& aln1, AlignmentResult::AlignmentItem& aln2) const
	{
		auto fwd = aln1.trace->trace.back().DPposition;
		auto bwd = aln2.trace->trace.back().DPposition;
		if (fwd.node % 2 == 1) std::swap(fwd, bwd);
		auto fwd_head_id = this->psiSeeder->get_graph_ptr()->id_by_coordinate(fwd.node / 2);
		auto bwd_head_id = this->psiSeeder->get_graph_ptr()->id_by_coordinate(bwd.node / 2);
		auto bwd_head_offset = this->psiSeeder->get_graph_ptr()->node_length(bwd_head_id) - bwd.nodeOffset - 1;
		return psiSeeder->verify_distance(fwd_head_id, fwd.nodeOffset, bwd_head_id, bwd_head_offset);
	}
	std::vector<std::vector<SeedHit>> getSeeds(const psi::seeder_type::readsrecord_type& chunk, psi::seeder_type::readsrecord_type& seeds, psi::seeder_type::traverser_type& traverser, const AlignerParams& params) const
	{
		assert(mode == Mode::Psi);
		std::vector<std::vector<SeedHit>> chunk_hits(params.psiChunkSize/2);
		psiSeeder->get_seeds(seeds, chunk, params.psiDistance);
		auto sindex = psiSeeder->index_reads(seeds);
		auto callback =
				[this, &chunk, &chunk_hits, &params](const auto& hit) {
					auto id = this->psiSeeder->get_graph_ptr()->coordinate_id(hit.node_id);
					bool reversed = hit.read_id % 2;
					auto read_offset = (reversed ? length(chunk.str[hit.read_id]) - hit.read_offset - 1 : hit.read_offset);
					auto node_offset = (reversed ? this->psiSeeder->get_graph_ptr()->node_length(hit.node_id) - hit.node_offset - 1 : hit.node_offset);
					chunk_hits[hit.read_id/2].push_back(SeedHit(id, node_offset, read_offset, params.psiLength, params.psiLength, reversed));
				};
		psiSeeder->seeds_all(seeds, sindex, traverser, callback);
		return chunk_hits;
	}
	std::vector<SeedHit> getSeeds(const std::string& seqName, const std::string& seq) const
	{
		switch(mode)
		{
			case Mode::File:
				assert(fileSeeds != nullptr);
				if (fileSeeds->count(seqName) == 0) return std::vector<SeedHit>{};
				return fileSeeds->at(seqName);
			case Mode::Mum:
				assert(mummerSeeder != nullptr);
				return mummerSeeder->getMumSeeds(seq, mumCount, mxmLength);
			case Mode::Mem:
				assert(mummerSeeder != nullptr);
				return mummerSeeder->getMemSeeds(seq, memCount, mxmLength);
			case Mode::Minimizer:
				assert(minimizerSeeder != nullptr);
				return minimizerSeeder->getSeeds(seq, minimizerSeedDensity);
			case Mode::Psi:
				assert(psiSeeder != nullptr);
				{
					std::vector<SeedHit> hits;
					auto callback_fwd =
							[this, &hits](const auto& hit) {
								auto id = this->psiSeeder->get_graph_ptr()->coordinate_id(hit.node_id);
								hits.push_back(SeedHit(id, hit.node_offset, hit.read_offset, hit.match_len, (hit.match_len/hit.gocc+1), false));
							};
					auto callback_rvs =
							[this, &hits, seqlen=seq.size()](const auto& hit) {
								auto id = this->psiSeeder->get_graph_ptr()->coordinate_id(hit.node_id);
								auto read_offset = seqlen - hit.read_offset - 1;
								auto node_offset = this->psiSeeder->get_graph_ptr()->node_length(hit.node_id) - hit.node_offset - 1;
								hits.push_back(SeedHit(id, node_offset, read_offset, hit.match_len, (hit.match_len/hit.gocc+1), true));
							};
					psiSeeder->seeds_on_paths(seq, callback_fwd);
					auto seq_rc = CommonUtils::ReverseComplement(seq);
					psiSeeder->seeds_on_paths(seq_rc, callback_rvs);
					return hits;
				}
			case Mode::None:
				assert(false);
		}
		return std::vector<SeedHit>{};
	}
};

struct AlignmentStats
{
	AlignmentStats() :
	reads(0),
	seeds(0),
	seedsFound(0),
	seedsExtended(0),
	readsWithASeed(0),
	alignments(0),
	fullLengthAlignments(0),
	readsWithAnAlignment(0),
	bpInReads(0),
	bpInReadsWithASeed(0),
	bpInAlignments(0),
	bpInFullAlignments(0),
	allAlignmentsCount(0),
	nofChunks(0),
	assertionBroke(false)
	{
	}
	std::atomic<size_t> reads;
	std::atomic<size_t> seeds;
	std::atomic<size_t> seedsFound;
	std::atomic<size_t> seedsExtended;
	std::atomic<size_t> readsWithASeed;
	std::atomic<size_t> alignments;
	std::atomic<size_t> fullLengthAlignments;
	std::atomic<size_t> readsWithAnAlignment;
	std::atomic<size_t> bpInReads;
	std::atomic<size_t> bpInReadsWithASeed;
	std::atomic<size_t> bpInAlignments;
	std::atomic<size_t> bpInFullAlignments;
	std::atomic<size_t> allAlignmentsCount;
	std::atomic<size_t> nofChunks;
	std::atomic<bool> assertionBroke;
};

bool is_file_exist(std::string fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}

void replaceDigraphNodeIdsWithOriginalNodeIds(vg::Alignment& alignment, const AlignmentGraph& graph)
{
	for (int i = 0; i < alignment.path().mapping_size(); i++)
	{
		int digraphNodeId = alignment.path().mapping(i).position().node_id();
		int originalNodeId = digraphNodeId / 2;
		alignment.mutable_path()->mutable_mapping(i)->mutable_position()->set_node_id(originalNodeId);
		std::string name = graph.OriginalNodeName(digraphNodeId);
		if (name.size() > 0)
		{
			alignment.mutable_path()->mutable_mapping(i)->mutable_position()->set_name(name);
		}
	}
}

void readFastqs(const std::vector<std::string>& filenames, moodycamel::ConcurrentQueue<std::pair<std::shared_ptr<FastQ>, std::shared_ptr<FastQ>>>& writequeue, std::atomic<bool>& readStreamingFinished)
{
	assertSetNoRead("Paired-end read streamer");
	std::pair<std::shared_ptr<FastQ>, std::shared_ptr<FastQ>> ends;
	ends.first = nullptr;
	for (auto filename : filenames)
	{
		assert(ends.first == nullptr);	// paired reads cannot be split across files
		FastQ::streamFastqFromFile(filename, false, [&writequeue, &ends](FastQ& read)
		{
			if (ends.first == nullptr) {
				ends.first = std::make_shared<FastQ>();
				std::swap(*ends.first, read);
				ends.first->seq_id.pop_back();  // remove delimiter
				ends.first->seq_id.pop_back();  // remove pair number
				return;
			}
			else {
				ends.second = std::make_shared<FastQ>();
				std::swap(*ends.second, read);
				ends.second->seq_id.pop_back();  // remove delimiter
				ends.second->seq_id.pop_back();  // remove pair number
			}
			size_t slept = 0;
			while (writequeue.size_approx() > 200)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				slept++;
				if (slept > 100) break;
			}
			writequeue.enqueue(ends);
			ends.first = nullptr;
		});
	}
	readStreamingFinished = true;
}

void readFastqs(const std::vector<std::string>& filenames, moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>>& writequeue, std::atomic<bool>& readStreamingFinished)
{
	assertSetNoRead("Read streamer");
	for (auto filename : filenames)
	{
		FastQ::streamFastqFromFile(filename, false, [&writequeue](FastQ& read)
		{
			std::shared_ptr<FastQ> ptr = std::make_shared<FastQ>();
			std::swap(*ptr, read);
			size_t slept = 0;
			while (writequeue.size_approx() > 200)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				slept++;
				if (slept > 100) break;
			}
			writequeue.enqueue(ptr);
		});
	}
	readStreamingFinished = true;
}

void consumeBytesAndWrite(const std::string& filename, moodycamel::ConcurrentQueue<std::string*>& writequeue, moodycamel::ConcurrentQueue<std::string*>& deallocqueue, std::atomic<bool>& allThreadsDone, std::atomic<bool>& allWriteDone, bool verboseMode, bool textMode)
{
	assertSetNoRead("Writer");
	auto openmode = std::ios::out;
	if (!textMode) openmode |= std::ios::binary;
	std::ofstream outfile { filename, openmode };

	bool wroteAny = false;

	std::string* alns[100] {};

	BufferedWriter coutoutput;
	if (verboseMode)
	{
		coutoutput = {std::cout};
	}

	while (true)
	{
		size_t gotAlns = writequeue.try_dequeue_bulk(alns, 100);
		if (gotAlns == 0)
		{
			if (!writequeue.try_dequeue(alns[0]))
			{
				if (allThreadsDone) break;
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				continue;
			}
			gotAlns = 1;
		}
		coutoutput << "write " << gotAlns << ", " << writequeue.size_approx() << " left" << BufferedWriter::Flush;
		for (size_t i = 0; i < gotAlns; i++)
		{
			outfile.write(alns[i]->data(), alns[i]->size());
		}
		deallocqueue.enqueue_bulk(alns, gotAlns);
		wroteAny = true;
	}

	if (!textMode && !wroteAny)
	{
		::google::protobuf::io::ZeroCopyOutputStream *raw_out =
					new ::google::protobuf::io::OstreamOutputStream(&outfile);
		::google::protobuf::io::GzipOutputStream *gzip_out =
					new ::google::protobuf::io::GzipOutputStream(raw_out);
		::google::protobuf::io::CodedOutputStream *coded_out =
					new ::google::protobuf::io::CodedOutputStream(gzip_out);
		coded_out->WriteVarint64(0);
		delete coded_out;
		delete gzip_out;
		delete raw_out;
	}

	allWriteDone = true;
}

void QueueInsertSlowly(moodycamel::ProducerToken& token, moodycamel::ConcurrentQueue<std::string*>& queue, std::string&& str)
{
	std::string* write = new std::string { std::move(str) };
	size_t waited = 0;
	while (!queue.try_enqueue(token, write) && !queue.try_enqueue(token, write))
	{
		if (queue.size_approx() < 100 && queue.enqueue(token, write)) break;
		std::this_thread::sleep_for(std::chrono::milliseconds(10));
		waited++;
		if (waited >= 10)
		{
			if (queue.enqueue(token, write)) break;
		}
	}
}

void writeGAMToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, moodycamel::ConcurrentQueue<std::string*>& alignmentsOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	::google::protobuf::io::ZeroCopyOutputStream *raw_out = new ::google::protobuf::io::OstreamOutputStream(&strstr);
	::google::protobuf::io::GzipOutputStream *gzip_out = new ::google::protobuf::io::GzipOutputStream(raw_out);
	::google::protobuf::io::CodedOutputStream *coded_out = new ::google::protobuf::io::CodedOutputStream(gzip_out);
	coded_out->WriteVarint64(alignments.alignments.size());
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].alignment != nullptr);
		std::string s;
		alignments.alignments[i].alignment->SerializeToString(&s);
		coded_out->WriteVarint32(s.size());
		coded_out->WriteRaw(s.data(), s.size());
	}
	delete coded_out;
	delete gzip_out;
	delete raw_out;
	QueueInsertSlowly(token, alignmentsOut, strstr.str());
}

void writeJSONToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, moodycamel::ConcurrentQueue<std::string*>& alignmentsOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	google::protobuf::util::JsonPrintOptions options;
	options.preserve_proto_field_names = true;
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].alignment != nullptr);
		std::string s;
		google::protobuf::util::MessageToJsonString(*alignments.alignments[i].alignment, &s, options);
		strstr << s;
		strstr << '\n';
	}
	QueueInsertSlowly(token, alignmentsOut, strstr.str());
}

void writeGAFToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, moodycamel::ConcurrentQueue<std::string*>& alignmentsOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].GAFline.size() > 0);
		strstr << alignments.alignments[i].GAFline;
		strstr << '\n';
	}
	QueueInsertSlowly(token, alignmentsOut, strstr.str());
}

void writeCorrectedToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, const std::string& readName, const std::string& original, size_t maxOverlap, moodycamel::ConcurrentQueue<std::string*>& correctedOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	zstr::ostream *compressed;
	if (params.compressCorrected)
	{
		compressed = new zstr::ostream(strstr);
	}
	std::vector<Correction> corrections;
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].corrected.size() > 0);
		corrections.emplace_back();
		corrections.back().startIndex = alignments.alignments[i].alignmentStart;
		corrections.back().endIndex = alignments.alignments[i].alignmentEnd;
		corrections.back().corrected = alignments.alignments[i].corrected;
	}
	std::string corrected = getCorrected(original, corrections, maxOverlap);
	if (params.compressCorrected)
	{
		(*compressed) << ">" << readName << std::endl;
		(*compressed) << corrected << std::endl;
		delete compressed;
	}
	else
	{
		strstr << ">" << readName << std::endl;
		strstr << corrected << std::endl;
	}
	QueueInsertSlowly(token, correctedOut, strstr.str());
}

void writeCorrectedClippedToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, moodycamel::ConcurrentQueue<std::string*>& correctedClippedOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	zstr::ostream *compressed;
	if (params.compressClipped)
	{
		compressed = new zstr::ostream(strstr);
	}
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].corrected.size() > 0);
		if (params.compressClipped)
		{
			(*compressed) << ">" << alignments.readName << "_" << i << "_" << alignments.alignments[i].alignmentStart << "_" << alignments.alignments[i].alignmentEnd << std::endl;
			(*compressed) << alignments.alignments[i].corrected << std::endl;
		}
		else
		{
			strstr << ">" << alignments.readName << "_" << i << "_" << alignments.alignments[i].alignmentStart << "_" << alignments.alignments[i].alignmentEnd << std::endl;
			strstr << alignments.alignments[i].corrected << std::endl;
		}
	}
	if (params.compressClipped)
	{
		delete compressed;
	}
	QueueInsertSlowly(token, correctedClippedOut, strstr.str());
}

template<typename TEntry>
void runComponentMappings(const AlignmentGraph& alignmentGraph, moodycamel::ConcurrentQueue<TEntry>& readFastqsQueue, std::atomic<bool>& readStreamingFinished, int threadnum, const Seeder& seeder, AlignerParams params, moodycamel::ConcurrentQueue<std::string*>& GAMOut, moodycamel::ConcurrentQueue<std::string*>& JSONOut, moodycamel::ConcurrentQueue<std::string*>& GAFOut, moodycamel::ConcurrentQueue<std::string*>& correctedOut, moodycamel::ConcurrentQueue<std::string*>& correctedClippedOut, moodycamel::ConcurrentQueue<std::string*>& deallocqueue, AlignmentStats& stats)
{
	moodycamel::ProducerToken GAMToken { GAMOut };
	moodycamel::ProducerToken JSONToken { JSONOut };
	moodycamel::ProducerToken GAFToken { GAFOut };
	moodycamel::ProducerToken correctedToken { correctedOut };
	moodycamel::ProducerToken clippedToken { correctedClippedOut };
	assertSetNoRead("Before any read");
	GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState reusableState { alignmentGraph, std::max(params.initialBandwidth, params.rampBandwidth), !params.highMemory };
	AlignmentSelection::SelectionOptions selectionOptions;
	selectionOptions.method = params.alignmentSelectionMethod;
	selectionOptions.graphSize = alignmentGraph.SizeInBP();
	selectionOptions.ECutoff = params.selectionECutoff;
	BufferedWriter cerroutput;
	BufferedWriter coutoutput;
	if (params.verboseMode)
	{
		cerroutput = {std::cerr};
		coutoutput = {std::cout};
	}
	psi::seeder_type::readsrecord_type chunk;
	psi::seeder_type::readsrecord_type seeds;
	psi::seeder_type::traverser_type traverser;
	TEntry entry;
	bool has_mate = true;
	std::shared_ptr<AlignmentResult> alignments = std::make_shared<AlignmentResult>();
	std::shared_ptr<AlignmentResult> alns_1st = nullptr;
	std::shared_ptr<AlignmentResult> alns_2nd = nullptr;
	std::vector<std::vector<SeedHit>> chunk_hits;
	auto thread_id = psi::get_thread_id();

	if (seeder.mode == Seeder::Mode::Psi)
	{
		chunk = seeder.psiSeeder->create_readrecord();
		seeds = seeder.psiSeeder->create_readrecord();
		traverser = seeder.psiSeeder->create_traverser();

		assert(params.psiChunkSize != 0);
		params.psiChunkSize = (2 - isSingle(entry)) * params.psiChunkSize * 2 /* reverse complements */;
	}

	std::size_t rid = 0;  /* rid=1..|chunk| -> processing; rid=0 -> stop; rid=|chunk|+1 -> start */
	while (true)
	{
		std::string* dealloc;
		while (deallocqueue.try_dequeue(dealloc))
		{
			delete dealloc;
		}
		std::shared_ptr<FastQ> fastq = nullptr;
		std::vector<SeedHit> hits;

		/* If the chunk size is met, start mapping the chunk */
		if (seeder.mode == Seeder::Mode::Psi && !rid && length(chunk) == params.psiChunkSize)
		{
			assert(params.psiChunkSize != 0);
			rid = length(chunk) + 1;
		}

		/* Toggle `has_mate` if paired end, otherwise `false`. */
		has_mate = !isSingle(entry) && !has_mate;

		if (rid)
		{
			if (rid == length(chunk) + 1) {
				auto chunkno = ++stats.nofChunks;
				std::cout << "PSI: Thread " << thread_id << " started processing chunk " << chunkno << " with " << length(chunk) << " reads" << std::endl;
				assert(chunk_hits.empty());
				auto timeStart = std::chrono::system_clock::now();
				chunk_hits = seeder.getSeeds(chunk, seeds, traverser, params);
				auto timeEnd = std::chrono::system_clock::now();
				size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
				std::cout << "PSI: Thread " << thread_id << " finished seed finding of chunk " << chunkno << " in " << time << "ms" << std::endl;
			}
			else if (rid == 1 && !has_mate)
			{
				clear(chunk);
				rid = 0;
				has_mate = true;
				continue;
			}
			fastq = getOtherEnd(entry);
			if (!has_mate)
			{
				FastQ read;
				rid -= 2;  /* NOTE: Also skips the reverse complement */
				read.seq_id = toCString(chunk.name[rid-1]);
				read.sequence = toCString(seqan::String<char, seqan::CStyle>(chunk.str[rid-1]));
				setOneEnd(entry, std::move(read));
				if (!isSingle(entry)) {
					rid -= 2;  /* NOTE: Also skips the reverse complement */
					read.seq_id = toCString(chunk.name[rid-1]);
					read.sequence = toCString(seqan::String<char, seqan::CStyle>(chunk.str[rid-1]));
					setOtherEnd(entry, std::move(read));
				}
				fastq = getOneEnd(entry);
			}
			else assert(fastq != nullptr);
			//getQualities(cur_read.quality, chunk.str[rid-1]);
			chunk_hits.back().swap(hits);
			chunk_hits.pop_back();
		}
		else
		{
			fastq = getOtherEnd(entry);
			if (!has_mate)
			{
				clear(entry);
				while (isEmpty(entry) && !readFastqsQueue.try_dequeue(entry))
				{
					bool tryBreaking = readStreamingFinished;
					if (!readFastqsQueue.try_dequeue(entry) && tryBreaking) break;
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
				}
				fastq = getOneEnd(entry);
			}
			else assert(fastq != nullptr);
		}

		if (fastq == nullptr)
		{
			has_mate = true;
			if (seeder.mode != Seeder::Mode::Psi || empty(chunk)) break;
			else rid = length(chunk) + 1;
			continue;
		}

		assertSetNoRead(fastq->seq_id);
		coutoutput << "Read " << fastq->seq_id << " size " << fastq->sequence.size() << "bp" << BufferedWriter::Flush;
		selectionOptions.readSize = fastq->sequence.size();

		size_t alntimems = 0;
		size_t clustertimems = 0;
		try
		{
			if (seeder.mode == Seeder::Mode::Psi)
			{
				if (!rid)
				{
					auto timeStart = std::chrono::system_clock::now();
					hits = seeder.getSeeds(fastq->seq_id, fastq->sequence);
					auto timeEnd = std::chrono::system_clock::now();
					size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
					coutoutput << "Read " << fastq->seq_id << " seeding took " << time << "ms" << BufferedWriter::Flush;
					if (hits.empty())
					{
						/* Append forward */
						appendValue(chunk.name, fastq->seq_id);
						appendValue(chunk.str, fastq->sequence);
						/* Append reverse complement */
						appendValue(chunk.name, fastq->seq_id);
						appendValue(chunk.str, CommonUtils::ReverseComplement(fastq->sequence));
						if (!isSingle(entry))
						{
							if (!has_mate)
							{
								fastq = getOtherEnd(entry);
								has_mate = true;
							}
							else
							{
								fastq = getOneEnd(entry);
								alns_1st = nullptr;  // discard alignments of the first end
							}
							/* Append forward */
							appendValue(chunk.name, fastq->seq_id);
							appendValue(chunk.str, fastq->sequence);
							/* Append reverse complement */
							appendValue(chunk.name, fastq->seq_id);
							appendValue(chunk.str, CommonUtils::ReverseComplement(fastq->sequence));
						}
						coutoutput << "Read " << fastq->seq_id << " added to the current chunk for traversal" << BufferedWriter::Flush;
						cerroutput << "Read " << fastq->seq_id << " added to the current chunk for traversal" << BufferedWriter::Flush;
					}
				}
				stats.seeds += hits.size();
				if (hits.size() == 0)
				{
					coutoutput << "Read " << fastq->seq_id << " has no seed hits" << BufferedWriter::Flush;
					cerroutput << "Read " << fastq->seq_id << " has no seed hits" << BufferedWriter::Flush;
					coutoutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
					cerroutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
					continue;
				}
				stats.seedsFound += hits.size();
				stats.readsWithASeed += 1;
				stats.bpInReadsWithASeed += fastq->sequence.size();
				auto clusterTimeStart = std::chrono::system_clock::now();
				{
					[[maybe_unused]] auto timer = seeder.psiSeeder->get_stats().timeit_ts("seed-cluster");
					OrderSeeds(alignmentGraph, hits);
				}
				auto clusterTimeEnd = std::chrono::system_clock::now();
				size_t clusterTime = std::chrono::duration_cast<std::chrono::milliseconds>(clusterTimeEnd - clusterTimeStart).count();
				coutoutput << "Read " << fastq->seq_id << " clustering took " << clusterTime << "ms" << BufferedWriter::Flush;
				auto alntimeStart = std::chrono::system_clock::now();
				{
					[[maybe_unused]] auto timer = seeder.psiSeeder->get_stats().timeit_ts("seed-extend");
					*alignments = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, params.maxCellsPerSlice, !params.verboseMode, !params.tryAllSeeds, hits, reusableState, !params.highMemory, params.forceGlobal, params.preciseClipping, params.seedClusterMinSize, params.seedExtendDensity, params.nondeterministicOptimizations);
				}
				auto alntimeEnd = std::chrono::system_clock::now();
				alntimems = std::chrono::duration_cast<std::chrono::milliseconds>(alntimeEnd - alntimeStart).count();
			}
			else if (seeder.mode != Seeder::Mode::None)
			{
				auto timeStart = std::chrono::system_clock::now();
				std::vector<SeedHit> hits = seeder.getSeeds(fastq->seq_id, fastq->sequence);
				auto timeEnd = std::chrono::system_clock::now();
				size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
				coutoutput << "Read " << fastq->seq_id << " seeding took " << time << "ms" << BufferedWriter::Flush;
				stats.seeds += hits.size();
				if (hits.size() == 0)
				{
					coutoutput << "Read " << fastq->seq_id << " has no seed hits" << BufferedWriter::Flush;
					cerroutput << "Read " << fastq->seq_id << " has no seed hits" << BufferedWriter::Flush;
					coutoutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
					cerroutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
					if (params.outputCorrectedFile != "") writeCorrectedToQueue(correctedToken, params, fastq->seq_id, fastq->sequence, alignmentGraph.getDBGoverlap(), correctedOut, *alignments);
					continue;
				}
				stats.seedsFound += hits.size();
				stats.readsWithASeed += 1;
				stats.bpInReadsWithASeed += fastq->sequence.size();
				auto clusterTimeStart = std::chrono::system_clock::now();
				OrderSeeds(alignmentGraph, hits);
				auto clusterTimeEnd = std::chrono::system_clock::now();
				size_t clusterTime = std::chrono::duration_cast<std::chrono::milliseconds>(clusterTimeEnd - clusterTimeStart).count();
				coutoutput << "Read " << fastq->seq_id << " clustering took " << clusterTime << "ms" << BufferedWriter::Flush;
				auto alntimeStart = std::chrono::system_clock::now();
				*alignments = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, params.maxCellsPerSlice, !params.verboseMode, !params.tryAllSeeds, hits, reusableState, !params.highMemory, params.forceGlobal, params.preciseClipping, params.seedClusterMinSize, params.seedExtendDensity, params.nondeterministicOptimizations);
				auto alntimeEnd = std::chrono::system_clock::now();
				alntimems = std::chrono::duration_cast<std::chrono::milliseconds>(alntimeEnd - alntimeStart).count();
			}
			else if (params.optimalDijkstra)
			{
				auto alntimeStart = std::chrono::system_clock::now();
				*alignments = AlignOneWayDijkstra(alignmentGraph, fastq->seq_id, fastq->sequence, !params.verboseMode, reusableState, params.forceGlobal, params.preciseClipping);
				auto alntimeEnd = std::chrono::system_clock::now();
				alntimems = std::chrono::duration_cast<std::chrono::milliseconds>(alntimeEnd - alntimeStart).count();
			}
			else
			{
				auto alntimeStart = std::chrono::system_clock::now();
				*alignments = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, !params.verboseMode, reusableState, !params.highMemory, params.forceGlobal, params.preciseClipping, params.nondeterministicOptimizations);
				auto alntimeEnd = std::chrono::system_clock::now();
				alntimems = std::chrono::duration_cast<std::chrono::milliseconds>(alntimeEnd - alntimeStart).count();
			}
		}
		catch (const ThreadReadAssertion::AssertionFailure& a)
		{
			coutoutput << "Read " << fastq->seq_id << " alignment failed (assertion!)" << BufferedWriter::Flush;
			cerroutput << "Read " << fastq->seq_id << " alignment failed (assertion!)" << BufferedWriter::Flush;
			reusableState.clear();
			stats.assertionBroke = true;
			continue;
		}

		stats.reads += 1;
		stats.bpInReads += fastq->sequence.size();
		stats.allAlignmentsCount += alignments->alignments.size();

		coutoutput << "Read " << fastq->seq_id << " alignment took " << alntimems << "ms" << BufferedWriter::Flush;
		if (alignments->alignments.size() > 0) alignments->alignments = AlignmentSelection::SelectAlignments(alignments->alignments, selectionOptions);

		//failed alignment, don't output
		if (alignments->alignments.size() == 0)
		{
			coutoutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			cerroutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			try
			{
				if (params.outputCorrectedFile != "") writeCorrectedToQueue(correctedToken, params, fastq->seq_id, fastq->sequence, alignmentGraph.getDBGoverlap(), correctedOut, *alignments);
			}
			catch (const ThreadReadAssertion::AssertionFailure& a)
			{
				reusableState.clear();
				stats.assertionBroke = true;
				continue;
			}
			if (isSingle(entry) || !has_mate || alns_1st == nullptr) continue;
		}
		else
		{
			stats.seedsExtended += alignments->seedsExtended;
			stats.readsWithAnAlignment += 1;

			size_t totalcells = 0;
			for (size_t i = 0; i < alignments->alignments.size(); i++)
			{
				totalcells += alignments->alignments[i].cellsProcessed;
			}
		
			std::sort(alignments->alignments.begin(), alignments->alignments.end(), [](const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right) { return left.alignmentStart < right.alignmentStart; });
		}

		if (!isSingle(entry) && seeder.mode == Seeder::Mode::Psi) // paired end
		{
			if (!has_mate)
			{
				assert(alns_1st == nullptr);
				alns_1st = alignments;
				alignments = std::make_shared<AlignmentResult>();
				continue;
			}
			if (alns_1st == nullptr)  // first mate had no alignment
			{
				alns_1st = std::make_shared<AlignmentResult>();
				alns_1st->readName = alns_2nd->readName;
			}
			/* else */
			assert(alns_2nd == nullptr);
			alns_2nd = alignments;
			alignments = std::make_shared<AlignmentResult>();
			alignments->seedsExtended = alns_1st->seedsExtended + alns_2nd->seedsExtended;
			alignments->readName = alns_1st->readName;
			if (alns_1st->readName != alns_2nd->readName) std::cerr << "The name of pairs are different (" << alns_1st->readName << " != " << alns_2nd->readName << ")" << std::endl;
			std::vector<std::pair<AlignmentResult::AlignmentItem, AlignmentResult::AlignmentItem>> paired_alignments;
			std::vector<bool> mark1(alns_1st->alignments.size(), 0);
			std::vector<bool> mark2(alns_2nd->alignments.size(), 0);
			auto mate_fastq = getOneEnd(entry);
			for (size_t i = 0; i < alns_1st->alignments.size(); ++i)
			{
				size_t alignmentSize = alns_1st->alignments[i].alignmentEnd - alns_1st->alignments[i].alignmentStart;
				if (alignmentSize == mate_fastq->sequence.size())
				{
					stats.fullLengthAlignments += 1;
					stats.bpInFullAlignments += alignmentSize;
				}
				stats.bpInAlignments += alignmentSize;
			}
			for (size_t j = 0; j < alns_2nd->alignments.size(); ++j)
			{
				size_t alignmentSize = alns_2nd->alignments[j].alignmentEnd - alns_2nd->alignments[j].alignmentStart;
				if (alignmentSize == fastq->sequence.size())
				{
					stats.fullLengthAlignments += 1;
					stats.bpInFullAlignments += alignmentSize;
				}
				stats.bpInAlignments += alignmentSize;
			}

			if (params.outputGAMFile != "" || params.outputJSONFile != "")
			{
				for (size_t i = 0; i < alns_1st->alignments.size(); ++i)
				{
					AddAlignment(mate_fastq->seq_id, mate_fastq->sequence, alns_1st->alignments[i]);
					replaceDigraphNodeIdsWithOriginalNodeIds(*alns_1st->alignments[i].alignment, alignmentGraph);
				}
				for (size_t j = 0; j < alns_2nd->alignments.size(); ++j)
				{
					AddAlignment(fastq->seq_id, fastq->sequence, alns_2nd->alignments[j]);
					replaceDigraphNodeIdsWithOriginalNodeIds(*alns_2nd->alignments[j].alignment, alignmentGraph);
				}
			}

			if (params.outputGAFFile != "")
			{
				for (size_t i = 0; i < alns_1st->alignments.size(); ++i)
				{
					AddGAFLine(alignmentGraph, mate_fastq->seq_id, mate_fastq->sequence, alns_1st->alignments[i]);
				}
				for (size_t j = 0; j < alns_2nd->alignments.size(); ++j)
				{
					AddGAFLine(alignmentGraph, fastq->seq_id, fastq->sequence, alns_2nd->alignments[j]);
				}
			}

			for (size_t i = 0; i < alns_1st->alignments.size(); ++i)
			{
				for (size_t j = 0; j < alns_2nd->alignments.size(); ++j)
				{
					if (seeder.verifyDistance(alns_1st->alignments[i], alns_2nd->alignments[j]))
					{
						mark1[i] = true;
						mark2[j] = true;
						auto p = std::make_pair(alns_1st->alignments[i], alns_2nd->alignments[j]);
						if (p.first.alignmentStart > p.second.alignmentStart) std::swap(p.first, p.second);
						if (params.outputGAMFile != "" || params.outputJSONFile != "")
						{
							p.first.alignment->mutable_fragment_next()->set_name(fastq->seq_id);
							p.second.alignment->mutable_fragment_prev()->set_name(mate_fastq->seq_id);
						}
						if (params.outputGAFFile != "")
						{
							p.first.GAFline += "\tfn:Z:" + fastq->seq_id;
							p.second.GAFline += "\tfp:Z:" + mate_fastq->seq_id;
						}
						paired_alignments.push_back(std::move(p));
					}
				}
			}

			std::sort(paired_alignments.begin(), paired_alignments.end(), [](const auto& l, const auto& r) { return l.first.alignmentStart < r.first.alignmentStart; });

			for (size_t i = 0; i < paired_alignments.size(); i++)
			{
				alignments->alignments.push_back(std::move(paired_alignments[i].first));
				alignments->alignments.push_back(std::move(paired_alignments[i].second));
			}
			for (size_t i = 0; i < alns_1st->alignments.size(); ++i)
			{
				if (mark1[i]) continue;
				alignments->alignments.push_back(alns_1st->alignments[i]);
			}
			for (size_t j = 0; j < alns_2nd->alignments.size(); ++j)
			{
				if (mark2[j]) continue;
				alignments->alignments.push_back(alns_2nd->alignments[j]);
			}

			std::sort(alignments->alignments.begin() + 2*paired_alignments.size(), alignments->alignments.end(), [](const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right) { return left.alignmentStart < right.alignmentStart; });

			alns_1st = nullptr;
			alns_2nd = nullptr;

			stats.alignments += alignments->alignments.size();

			try
			{
				if (params.outputGAMFile != "") writeGAMToQueue(GAMToken, params, GAMOut, *alignments);
				if (params.outputJSONFile != "") writeJSONToQueue(JSONToken, params, JSONOut, *alignments);
				if (params.outputGAFFile != "") writeGAFToQueue(GAFToken, params, GAFOut, *alignments);
				if (params.outputCorrectedFile != "") std::runtime_error("not implemented for paired-end reads");
				if (params.outputCorrectedClippedFile != "") writeCorrectedClippedToQueue(clippedToken, params, correctedClippedOut, *alignments);
			}
			catch (const ThreadReadAssertion::AssertionFailure& a)
			{
				reusableState.clear();
				stats.assertionBroke = true;
				continue;
			}
		}
		else
		{
			if (params.outputGAMFile != "" || params.outputJSONFile != "")
			{
				for (size_t i = 0; i < alignments->alignments.size(); i++)
				{
					AddAlignment(fastq->seq_id, fastq->sequence, alignments->alignments[i]);
					replaceDigraphNodeIdsWithOriginalNodeIds(*alignments->alignments[i].alignment, alignmentGraph);
				}
			}

			if (params.outputGAFFile != "")
			{
				for (size_t i = 0; i < alignments->alignments.size(); i++)
				{
					AddGAFLine(alignmentGraph, fastq->seq_id, fastq->sequence, alignments->alignments[i]);
				}
			}
		
			std::sort(alignments->alignments.begin(), alignments->alignments.end(), [](const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right) { return left.alignmentStart < right.alignmentStart; });

			std::string alignmentpositions;

			for (size_t i = 0; i < alignments->alignments.size(); i++)
			{
				stats.alignments += 1;
				size_t alignmentSize = alignments->alignments[i].alignmentEnd - alignments->alignments[i].alignmentStart;
				if (alignmentSize == fastq->sequence.size())
				{
					stats.fullLengthAlignments += 1;
					stats.bpInFullAlignments += alignmentSize;
				}
				stats.bpInAlignments += alignmentSize;
				if (params.outputCorrectedFile != "" || params.outputCorrectedClippedFile != "") AddCorrected(alignments->alignments[i]);
				alignmentpositions += std::to_string(alignments->alignments[i].alignmentStart) + "-" + std::to_string(alignments->alignments[i].alignmentEnd) + ", ";
			}

			alignmentpositions.pop_back();
			alignmentpositions.pop_back();
			coutoutput << "Read " << fastq->seq_id << " aligned by thread " << threadnum << " with positions: " << alignmentpositions << " (read " << fastq->sequence.size() << "bp)" << BufferedWriter::Flush;

			try
			{
				if (params.outputGAMFile != "") writeGAMToQueue(GAMToken, params, GAMOut, *alignments);
				if (params.outputJSONFile != "") writeJSONToQueue(JSONToken, params, JSONOut, *alignments);
				if (params.outputGAFFile != "") writeGAFToQueue(GAFToken, params, GAFOut, *alignments);
				if (params.outputCorrectedFile != "") writeCorrectedToQueue(correctedToken, params, fastq->seq_id, fastq->sequence, alignmentGraph.getDBGoverlap(), correctedOut, *alignments);
				if (params.outputCorrectedClippedFile != "") writeCorrectedClippedToQueue(clippedToken, params, correctedClippedOut, *alignments);
			}
			catch (const ThreadReadAssertion::AssertionFailure& a)
			{
				reusableState.clear();
				stats.assertionBroke = true;
				continue;
			}
		}
	}
	assertSetNoRead("After all reads");
	coutoutput << "Thread " << threadnum << " finished" << BufferedWriter::Flush;
}

AlignmentGraph getGraph(std::string graphFile, MummerSeeder** mxmSeeder, psi::graph_type& seedGraph, psi::seeder_type** psiseeder, const AlignerParams& params)
{
	bool loadMxmSeeder = params.mumCount > 0 || params.memCount > 0;
	bool loadPsi = params.psiLength != 0;
	if (is_file_exist(graphFile)){
		std::cout << "Load graph from " << graphFile << std::endl;
	}
	else{
		std::cerr << "No graph file exists" << std::endl;
		std::exit(0);
	}
	try
	{
		if (loadPsi)
		{
			[[maybe_unused]] auto timer = psi::seeder_type::stats_type::timer_type("psi-prepare");
			std::cout << "Initializing Pan-genome Seed Index (PSI)" << std::endl;
			gum::util::load(seedGraph, graphFile, true);
			*psiseeder = new psi::seeder_type(seedGraph, params.psiLength, params.psiGoccThreshold, params.psiMaxMem);
			std::signal(SIGUSR1, psi::seeder_type::stats_type::signal_handler);
#ifdef PSI_STATS
			std::cout << "PSI: Custom signal handler for SIGUSR1 has been set" << std::endl;
			std::cout << "PSI: You can get the PSI::SeedFinder status report by sending SIGUSR1" << std::endl;
#endif
			/* Load the genome-wide path index for the graph if available. */
			if ((*psiseeder)->load_path_index(params.seederCachePrefix, params.psiContext, params.psiStep, params.readMinInsertSize, params.readMaxInsertSize))
			{
				std::cout << "PSI: Loaded existing path index" << std::endl;
			}
			/* No genome-wide path index requested. */
			else if (params.psiPathCount == 0)
			{
				std::cout << "PSI: Skip path indexing since no path specified" << std::endl;
			}
			else
			{
				std::cout << "PSI: Create a path index with " << params.psiPathCount << " paths" << std::endl;
				(*psiseeder)->create_path_index(params.psiPathCount, /*patched=*/true, params.psiContext, params.psiStep, params.readMinInsertSize, params.readMaxInsertSize,
																				[](auto msg){
																					std::cout << msg << std::endl;
																				},
																				[](auto msg){
																					std::cerr << msg << std::endl;
																				});
				/* Serialize the indexed paths. */
				if (!params.seederCachePrefix.empty())
				{
					(*psiseeder)->serialize_path_index(params.seederCachePrefix, params.psiStep);
				}
			}
			std::cout << "PSI: Number of uncovered loci: " << (*psiseeder)->get_starting_loci().size() << std::endl;
		}
		if (graphFile.substr(graphFile.size()-3) == ".vg")
		{
			if (loadMxmSeeder)
			{
				auto graph = CommonUtils::LoadVGGraph(graphFile);
				if (loadMxmSeeder)
				{
					std::cout << "Build MUM/MEM seeder from the graph" << std::endl;
					*mxmSeeder = new MummerSeeder { graph, params.seederCachePrefix };
				}
				std::cout << "Build alignment graph" << std::endl;
				auto result = DirectedGraph::BuildFromVG(graph);
				return result;
			}
			else
			{
				return DirectedGraph::StreamVGGraphFromFile(graphFile);
			}
		}
		else if (graphFile.substr(graphFile.size() - 4) == ".gfa")
		{
			auto graph = GfaGraph::LoadFromFile(graphFile, true);
			if (loadMxmSeeder)
			{
				std::cout << "Build MUM/MEM seeder from the graph" << std::endl;
				*mxmSeeder = new MummerSeeder { graph, params.seederCachePrefix };
			}
			std::cout << "Build alignment graph" << std::endl;
			auto result = DirectedGraph::BuildFromGFA(graph);
			return result;
		}
		else
		{
			std::cerr << "Unknown graph type (" << graphFile << ")" << std::endl;
			std::exit(0);
		}
	}
	catch (const CommonUtils::InvalidGraphException& e)
	{
		std::cout << "Error in the graph: " << e.what() << std::endl;
		std::cerr << "Error in the graph: " << e.what() << std::endl;
		std::exit(1);
	}
}

template<typename TSpec>
void alignReads(AlignerParams params, TSpec)
{
	using entry_type = typename Entry<TSpec, FastQ>::type;

	assertSetNoRead("Preprocessing");
	auto preProcStart = std::chrono::system_clock::now();

	const std::unordered_map<std::string, std::vector<SeedHit>>* seedHitsToThreads = nullptr;
	std::unordered_map<std::string, std::vector<SeedHit>> seedHits;
	MummerSeeder* mummerseeder = nullptr;
	psi::seeder_type* psiseeder = nullptr;
	psi::graph_type seedGraph;
	auto alignmentGraph = getGraph(params.graphFile, &mummerseeder, seedGraph, &psiseeder, params);
	bool loadMinimizerSeeder = params.minimizerSeedDensity != 0;
	MinimizerSeeder* minimizerseeder = nullptr;
	if (loadMinimizerSeeder)
	{
		std::cout << "Build minimizer seeder from the graph" << std::endl;
		minimizerseeder = new MinimizerSeeder(alignmentGraph, params.minimizerLength, params.minimizerWindowSize, params.numThreads, 1.0 - params.minimizerDiscardMostNumerousFraction);
		if (!minimizerseeder->canSeed())
		{
			std::cout << "Warning: Minimizer seeder has no seed hits. Reads cannot be aligned. Try unchopping the graph with vg or a different seeding mode" << std::endl;
		}
	}

	if (params.seedFiles.size() > 0)
	{
		for (auto file : params.seedFiles)
		{
			if (is_file_exist(file)){
				std::cout << "Load seeds from " << file << std::endl;
				std::ifstream seedfile { file, std::ios::in | std::ios::binary };
				size_t numSeeds = 0;
				std::function<void(vg::Alignment&)> alignmentLambda = [&seedHits, &numSeeds](vg::Alignment& seedhit) {
					seedHits[seedhit.name()].emplace_back(seedhit.path().mapping(0).position().node_id(), seedhit.path().mapping(0).position().offset(), seedhit.query_position(), seedhit.path().mapping(0).edit(0).from_length(), seedhit.path().mapping(0).edit(0).from_length(), seedhit.path().mapping(0).position().is_reverse());
					numSeeds += 1;
				};
				stream::for_each(seedfile, alignmentLambda);
				std::cout << numSeeds << " seeds" << std::endl;
			}
			else {
				std::cerr << "No seeds file exists" << std::endl;
				std::exit(0);
			}
		}
		seedHitsToThreads = &seedHits;
	}

	Seeder seeder { params, seedHitsToThreads, mummerseeder, minimizerseeder, psiseeder };

	switch(seeder.mode)
	{
		case Seeder::Mode::File:
			std::cout << "Seeds from file" << std::endl;
			break;
		case Seeder::Mode::Psi:
			std::cout << "PSI seeds, seed length " << params.psiLength
								<< ", chunk size " << params.psiChunkSize
								<< ", distance " << params.psiDistance
								<< ", paths " << params.psiPathCount
								<< ", context " << params.psiContext
								<< ", step " << params.psiStep
								<< ", GOCC threshold " << params.psiGoccThreshold
								<< ", maximum MEMs on paths " << params.psiMaxMem
								<< ", distance index minimum insert size " << params.readMinInsertSize
								<< ", distance index maximum insert size " << params.readMaxInsertSize
								<< std::endl;
			break;
		case Seeder::Mode::Mum:
			std::cout << "MUM seeds, min length " << seeder.mxmLength;
			if (seeder.mumCount != std::numeric_limits<size_t>::max()) std::cout << ", max count " << seeder.mumCount;
			std::cout << std::endl;
			break;
		case Seeder::Mode::Mem:
			std::cout << "MEM seeds, min length " << seeder.mxmLength;
			if (seeder.memCount != std::numeric_limits<size_t>::max()) std::cout << ", max count " << seeder.memCount;
			std::cout << std::endl;
			break;
		case Seeder::Mode::Minimizer:
			std::cout << "Minimizer seeds, length " << seeder.minimizerLength << ", window size " << seeder.minimizerWindowSize << ", density " << seeder.minimizerSeedDensity << std::endl;
			break;
		case Seeder::Mode::None:
			if (params.optimalDijkstra)
			{
				std::cout << "Optimal alignment. VERY SLOW!" << std::endl;
			}
			else
			{
				std::cout << "No seeds, calculate the entire first row. VERY SLOW!" << std::endl;
			}
			break;
	}
	if (seeder.mode != Seeder::Mode::None) std::cout << "Seed cluster size " << params.seedClusterMinSize << std::endl;
	if (seeder.mode != Seeder::Mode::None && params.seedExtendDensity != -1) std::cout << "Extend up to best " << params.seedExtendDensity << " fraction of seeds" << std::endl;

	if (!params.optimalDijkstra) std::cout << "Initial bandwidth " << params.initialBandwidth;
	if (params.rampBandwidth > 0) std::cout << ", ramp bandwidth " << params.rampBandwidth;
	if (params.maxCellsPerSlice != std::numeric_limits<size_t>::max()) std::cout << ", tangle effort " << params.maxCellsPerSlice;
	std::cout << std::endl;

	if (params.selectionECutoff != -1) std::cout << "Discard alignments with an E-value > " << params.selectionECutoff << std::endl;

	if (params.outputGAMFile != "") std::cout << "write alignments to " << params.outputGAMFile << std::endl;
	if (params.outputJSONFile != "") std::cout << "write alignments to " << params.outputJSONFile << std::endl;
	if (params.outputGAFFile != "") std::cout << "write alignments to " << params.outputGAFFile << std::endl;
	if (params.outputCorrectedFile != "") std::cout << "write corrected reads to " << params.outputCorrectedFile << std::endl;
	if (params.outputCorrectedClippedFile != "") std::cout << "write corrected & clipped reads to " << params.outputCorrectedClippedFile << std::endl;

	std::vector<std::thread> threads;

	auto preProcEnd = std::chrono::system_clock::now();
	auto preprocms = std::chrono::duration_cast<std::chrono::milliseconds>(preProcEnd - preProcStart).count();
	std::cout << "Preprocessing took " << preprocms << "ms" << std::endl;

	assertSetNoRead("Running alignments");

	moodycamel::ConcurrentQueue<std::string*> outputGAM { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<std::string*> outputGAF { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<std::string*> outputJSON { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<std::string*> deallocAlns;
	moodycamel::ConcurrentQueue<std::string*> outputCorrected { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<std::string*> outputCorrectedClipped { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<entry_type> readFastqsQueue;
	std::atomic<bool> readStreamingFinished { false };
	std::atomic<bool> allThreadsDone { false };
	std::atomic<bool> GAMWriteDone { false };
	std::atomic<bool> GAFWriteDone { false };
	std::atomic<bool> JSONWriteDone { false };
	std::atomic<bool> correctedWriteDone { false };
	std::atomic<bool> correctedClippedWriteDone { false };

	std::cout << "Align" << std::endl;
	AlignmentStats stats;
	std::thread fastqThread { [files=params.fastqFiles, &readFastqsQueue, &readStreamingFinished]() { readFastqs(files, readFastqsQueue, readStreamingFinished); } };
	std::thread GAMwriterThread { [file=params.outputGAMFile, &outputGAM, &deallocAlns, &allThreadsDone, &GAMWriteDone, verboseMode=params.verboseMode]() { if (file != "") consumeBytesAndWrite(file, outputGAM, deallocAlns, allThreadsDone, GAMWriteDone, verboseMode, false); else GAMWriteDone = true; } };
	std::thread GAFwriterThread { [file=params.outputGAFFile, &outputGAF, &deallocAlns, &allThreadsDone, &GAFWriteDone, verboseMode=params.verboseMode]() { if (file != "") consumeBytesAndWrite(file, outputGAF, deallocAlns, allThreadsDone, GAFWriteDone, verboseMode, false); else GAFWriteDone = true; } };
	std::thread JSONwriterThread { [file=params.outputJSONFile, &outputJSON, &deallocAlns, &allThreadsDone, &JSONWriteDone, verboseMode=params.verboseMode]() { if (file != "") consumeBytesAndWrite(file, outputJSON, deallocAlns, allThreadsDone, JSONWriteDone, verboseMode, true); else JSONWriteDone = true; } };
	std::thread correctedWriterThread { [file=params.outputCorrectedFile, &outputCorrected, &deallocAlns, &allThreadsDone, &correctedWriteDone, verboseMode=params.verboseMode, uncompressed=!params.compressCorrected]() { if (file != "") consumeBytesAndWrite(file, outputCorrected, deallocAlns, allThreadsDone, correctedWriteDone, verboseMode, uncompressed); else correctedWriteDone = true; } };
	std::thread correctedClippedWriterThread { [file=params.outputCorrectedClippedFile, &outputCorrectedClipped, &deallocAlns, &allThreadsDone, &correctedClippedWriteDone, verboseMode=params.verboseMode, uncompressed=!params.compressClipped]() { if (file != "") consumeBytesAndWrite(file, outputCorrectedClipped, deallocAlns, allThreadsDone, correctedClippedWriteDone, verboseMode, uncompressed); else correctedClippedWriteDone = true; } };

	{
		[[maybe_unused]] auto timer = psi::seeder_type::stats_type::timer_type("alignment");
		for (size_t i = 0; i < params.numThreads; i++)
		{
			threads.emplace_back([&alignmentGraph, &readFastqsQueue, &readStreamingFinished, i, seeder, params, &outputGAM, &outputJSON, &outputGAF, &outputCorrected, &outputCorrectedClipped, &deallocAlns, &stats]() { runComponentMappings(alignmentGraph, readFastqsQueue, readStreamingFinished, i, seeder, params, outputGAM, outputJSON, outputGAF, outputCorrected, outputCorrectedClipped, deallocAlns, stats); });
		}

		for (size_t i = 0; i < params.numThreads; i++)
		{
			threads[i].join();
		}
	}

	assertSetNoRead("Postprocessing");

	allThreadsDone = true;

	GAMwriterThread.join();
	GAFwriterThread.join();
	JSONwriterThread.join();
	correctedWriterThread.join();
	correctedClippedWriterThread.join();
	fastqThread.join();

	if (mummerseeder != nullptr) delete mummerseeder;
	if (minimizerseeder != nullptr) delete minimizerseeder;
	if (psiseeder != nullptr) delete psiseeder;

	std::string* dealloc;
	while (deallocAlns.try_dequeue(dealloc))
	{
		delete dealloc;
	}

	std::cout << "Alignment finished" << std::endl;
	std::cout << "Input reads: " << stats.reads << " (" << stats.bpInReads << "bp)" << std::endl;
	std::cout << "Seeds found: " << stats.seedsFound << std::endl;
	std::cout << "Seeds extended: " << stats.seedsExtended << std::endl;
	std::cout << "Reads with a seed: " << stats.readsWithASeed << " (" << stats.bpInReadsWithASeed << "bp)" << std::endl;
	std::cout << "Reads with an alignment: " << stats.readsWithAnAlignment << std::endl;
	std::cout << "Alignments: " << stats.alignments << " (" << stats.bpInAlignments << "bp)";
	if (stats.allAlignmentsCount > stats.alignments) std::cout << " (" << (stats.allAlignmentsCount - stats.alignments) << " additional alignments discarded)";
	std::cout << std::endl;
	std::cout << "End-to-end alignments: " << stats.fullLengthAlignments << " (" << stats.bpInFullAlignments << "bp)" << std::endl;
	std::cout << "Chunks processed: " << stats.nofChunks << " (in " << params.numThreads << " threads)" << std::endl;
	if (stats.assertionBroke)
	{
		std::cout << "Alignment broke with some reads. Look at stderr output." << std::endl;
	}

	for (const auto& timer : psi::seeder_type::stats_type::timer_type::get_timers())
	{
		std::cout << "PSI timer '" << timer.first << "': " << timer.second.str() << std::endl;
	}
}

template void alignReads< SingleEnd >(AlignerParams params, SingleEnd);
template void alignReads< PairedEnd >(AlignerParams params, PairedEnd);
