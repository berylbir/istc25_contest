#ifndef conv_H
#define conv_H

#include "FeedForwardTrellis.h"
#include "consts.h"
#include "tbcc_types.h"
#include "minHeap.h"
#include "tbcc_namespace.h"
#include "lowRateListDecoder.h"
#include <string>
#include <vector>

using intvec = std::vector<int>;
using fltvec = std::vector<float>;
using dblvec = std::vector<double>;

// Class for generating, encoding, and decoding binary conv codes
class conv
{
  private:
    // basic parameters
    int k_, n_, v_;
    std::vector<int> numerators_;

    /* - Trellis setup - */
    FeedForwardTrellis trellis;

    /* - Decoder setup - */
    LowRateListDecoder decoder;

  public:

    // Sparse binary matrix given by list row nad col positions of non-zero elements
    intvec row, col;
    
    // Parity generator matrix
    std::vector<std::vector<int>> parity_generator;

    // Constructor
    conv() : k_(k), n_(n), v_(V), numerators_(numerators), 
              trellis(k, n, V, numerators), decoder(trellis, 1e7, M+1, CRC, 'T'){}


    // Generate encoder
    void create_encoder(int verbose = 0);

    // Read in offset lists
    void readFileToRows(const std::string& filename, RowList& rows);
    void readFiles(const std::string& folderPath, RowList& msg, RowList& cwd, int lowDist, int highDist);

    // SLVD decoding
    MessageInformation decode(fltvec &llr_in, intvec punctured_indices, int verbose);
 
    // Encode K info bits into N codeword bits
    void encode(intvec &info, intvec &cw);
};

// MORE TEST FUNCTIONS
void test_no_error(conv &code);
// void test_single_error(conv &code);
// void test_gaussian_noise(conv &code, float esno);

#endif // conv_H

