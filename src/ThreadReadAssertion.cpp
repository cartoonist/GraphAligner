#include <iostream>
#include <sstream>
#include <string_view>
#include "ThreadReadAssertion.h"

namespace ThreadReadAssertion
{
	thread_local int currentnodeID;
	thread_local bool currentreverse;
	thread_local size_t currentseqPos;
	thread_local size_t currentmatchLen;
	thread_local size_t currentnodeOffset;
	thread_local std::string_view currentRead;
	void signal(int signal)
	{
		std::stringstream msg;
		msg << "Signal " << signal << ". Read: " << currentRead << ". Seed: " << assertGetSeedInfo();
		std::cerr << msg.str() << std::endl;
		std::abort();
	}

	void bt_sighandler(int sig, siginfo_t* info, void* secret)
	{
		void* trace[16];
		char** messages = NULL;
		int i, trace_size = 0;
		ucontext_t* uc = (ucontext_t*)secret;

		char* sigstr = strsignal(sig);
		void* regip = NULL;
#ifdef __x86_64__
		regip = (void*) uc->uc_mcontext.gregs[ REG_RIP ];
#else
		regip = (void*) uc->uc_mcontext.gregs[ REG_EIP ];
#endif
		if (sig == SIGSEGV)
		{
			printf("Got signal '%s' (%d), faulty address is %p, from %p\n", sigstr, sig,
							info->si_addr, regip);
		}
		else printf("Got signal '%s' (%d)\n", sigstr, sig);

		trace_size = backtrace(trace, 16);
		/* overwrite sigaction with caller's address */
		trace[1] = regip;

		messages = backtrace_symbols(trace, trace_size);
		/* skip first stack frame (points here) */
		for (i = 1; i < trace_size; ++i)
		{
			Dl_info dlinfo;
			if (dladdr(trace[i], &dlinfo))
			{
				printf("[bt] %-3d %s\n", i, messages[i]);
				char syscom[256];
				sprintf(syscom, "addr2line -fCpe %s %p | sed -e 's/^/[bt]     /'",
								 dlinfo.dli_fname,
								 (void*)((char*) trace[i] - (char*) dlinfo.dli_fbase));
				system(syscom);
			}
			else printf("[bt] %-3d %s\n", i, messages[i]);
		}
		std::cerr << "Read: " << currentRead << ". Seed: " << assertGetSeedInfo() << std::endl;
	}

	void setRead(const std::string& readName)
	{
		currentRead = std::string_view(readName.data(), readName.size());
	}
	void setSeed(int nodeID, bool reverse, size_t seqPos, size_t matchLen, size_t nodeOffset)
	{
		currentnodeID = nodeID;
		currentreverse = reverse;
		currentseqPos = seqPos;
		currentmatchLen = matchLen;
		currentnodeOffset = nodeOffset;
	}
	void assertFailed(const char* expression, const char* file, int line)
	{
		std::stringstream msg;
		msg << file << ":" << line << ": Assertion '" << expression << "' failed. Read: " << currentRead << ". Seed: " << assertGetSeedInfo();
		std::cerr << msg.str() << std::endl;
		throw AssertionFailure {};
	}	
	std::string assertGetSeedInfo()
	{
		return std::to_string(currentnodeID) + (currentreverse ? "-" : "+") + "," + std::to_string(currentseqPos) + "," + std::to_string(currentmatchLen) + "," + std::to_string(currentnodeOffset);
	}
}
