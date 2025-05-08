#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <set>
#include <vector>
#include <random>
#include <fstream>
#include "conv.h"

// Generate conv encoder
void conv::create_encoder(int verbose) {
    // std::cout << "K = " << K << std::endl;
}

// Read in offset tables
void conv::readFileToRows(const std::string& filename, RowList& rows) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error: Unable to open " << filename << "\n";
        return;
    }

    std::string line;
    int rowIndex = 0;

    while (std::getline(infile, line) && rowIndex < numPairs) {
        std::istringstream iss(line);
        int value;
        while (iss >> value) {
            if (value != 3) {
                rows[rowIndex].push_back(value);
            }
        }
        ++rowIndex;
    }
}

// Main function to read messages and codewords
void conv::readFiles(const std::string& folderPath, RowList& msg, RowList& cwd, int lowDist, int highDist) {
    for (int fileNum = lowDist; fileNum <= highDist; fileNum += 4) {
        std::string filename_m = "m_" + folderPath + "/" + std::to_string(fileNum) + ".txt";
        std::string filename_c = "c_" + folderPath + "/" + std::to_string(fileNum) + ".txt";

        readFileToRows(filename_m, msg);
        readFileToRows(filename_c, cwd);

        std::cout << "read in offset " << fileNum << "\n";
    }
}


// SLVD decoding
MessageInformation conv::decode(fltvec &llr_in, intvec punctured_indices, int verbose = 0) {
    // TBCC
    if (CODE_TYPE == 'T'){
        return decoder.SSD_SLVD_TB(llr_in, punctured_indices);
    }
    else if (CODE_TYPE == 'Z'){
        return decoder.SSD_SLVD_ZT(llr_in);
    }
    // returns the correct result we get from squared distance metric
    // return decoder.decode(llr_in, punctured_indices);
}
    


// Encode info bitvec into codeword bitvec
void conv::encode(intvec &info, intvec &cw) {
    if (CODE_TYPE == 'T'){
        cw = trellis.encode(info);
    }
    else if (CODE_TYPE == 'Z'){
        // Append ZT bits
        for (int i=0; i<V; i++){
            info.push_back(0);
        }
        // ZT encoder
        cw = trellis.encode_zt(info);
    }
    
}

