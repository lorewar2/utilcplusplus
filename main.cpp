#include "include/kmc_runner.h"
#include <iostream>

int main()
{
    try
    {
        std::vector<std::string> inputFiles { "/data1/phasstphase_test/potato/hera/SRR15115655_1.fastq.gz", "/data1/phasstphase_test/potato/hera/SRR15115655_2.fastq.gz"};
        KMC::Stage1Params stage1Params;

        stage1Params.SetInputFiles(inputFiles)
            .SetKmerLen(31)
            .SetNThreads(8)
            .SetMaxRamGB(10)
            .SetEstimateHistogramCfg(KMC::EstimateHistogramCfg::ESTIMATE_AND_COUNT_KMERS);

        KMC::Runner kmc;

        auto stage1Results = kmc.RunStage1(stage1Params);

        uint32_t cutoffMin = 5;

        uint64_t nUniqCountedKmersEst{};
        for(uint32_t i = cutoffMin ; i < stage1Results.estimatedHistogram.size() ; ++i)
            nUniqCountedKmersEst += stage1Results.estimatedHistogram[i];
        
        std::cout << "#uniq counted kmers estimate: " << nUniqCountedKmersEst << "\n";

        //now lests assume the total amount of memory that we allow KMC to use depends on the total number of unique k-mers
        //and per each unique k-mer we want at most bitsPerUniqieKmers
        double bitsPerUniqieKmers = 20;
        double ramBits = bitsPerUniqieKmers * nUniqCountedKmersEst;
        double ramBytes = ramBits / 8;
        uint32_t ramForStage2  = ramBytes / 1000 / 1000 / 1000;    

        if (ramForStage2 < 2)
            ramForStage2 = 2; //at least 2 GB needed for kmc

        std::cout << "ram for Stage 2 (GB): " << ramForStage2 << "\n";

        KMC::Stage2Params stage2Params;

        stage2Params.SetNThreads(8)
            .SetMaxRamGB(ramForStage2)
            .SetCutoffMin(cutoffMin)
            .SetOutputFileName("kmers").
            SetStrictMemoryMode(true);

        auto stage2Results = kmc.RunStage2(stage2Params);

        std::cout << "#total counted k-mers: " << stage2Results.nTotalKmers << "\n";
        std::cout << "#unique k-mers: " << stage2Results.nUniqueKmers << "\n";
        std::cout << "#unique counted k-mers: " << stage2Results.nUniqueKmers - stage2Results.nBelowCutoffMin - stage2Results.nAboveCutoffMax <<"\n";
        std::cout << "#sequences: " << stage1Results.nSeqences << "\n";


    }
    catch(const std::runtime_error& err)
    {
        std::cerr << err.what() << "\n";
    }
}