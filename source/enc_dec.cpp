#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include "enc_dec.h"
#include "conv.h"

conv code;

RowList messageList(numPairs);
RowList codewordList(numPairs);

// Setup for [n,k] code
int enc_dec::init(int K, int N, bool opt_avg_latency) {
    // Contestants should replace this code
    //   This code should initialize the encoder-decoder
    std::cout << "k = " << k << ", n = " << n << ", v = " << V << std::endl;
    std::cout << "message length: " << K << ", blocklength: " << N << std::endl;
    // Setup encoder
    code.create_encoder();

    // Setup offset tables for decoding
    std::string filename = "k" + std::to_string(k) + 
                       "n" + std::to_string(n) + 
                       "v" + std::to_string(V) + 
                       "m" + std::to_string(M) + 
                       "K" + std::to_string(K);
    code.readFiles(filename, messageList, codewordList, lowDist, highDist);

    // Decoding iterations
    return 0;
}

// Encode K info bits into N codeword bits
void enc_dec::encode(bitvec &info, bitvec &cw) {

    // add crc
	crc::crc_calculation(info, M+1, CRC);
    // Encode
    bitvec cw_unpunctured(N + PUNCTURING_INDICES.size());
    code.encode(info, cw_unpunctured);
    
    // std::cout << "printing unpunctured encoded word: ";
    // utils::print_int_vector(cw_unpunctured);
    // std::cout << std::endl;

    // puncture
    bitvec cw_punctured = {};
    for (int i = 0; i < cw_unpunctured.size(); i++) {
        if (std::find(PUNCTURING_INDICES.begin(), PUNCTURING_INDICES.end(), i) == PUNCTURING_INDICES.end()) {
            cw_punctured.push_back(cw_unpunctured[i]);
        }
    }
    cw = cw_punctured;
}

// Decode N llrs into N codeword bits and K info bits, return -1 if detected error
int enc_dec::decode(fltvec &llr, bitvec &cw_est, bitvec &info_est) {
    
    fltvec unpunctured_symbols(N + PUNCTURING_INDICES.size());
    auto llr_ptr = llr.begin();
    for (int i = 0; i < unpunctured_symbols.size(); i++) {
        if (find(PUNCTURING_INDICES.begin(), PUNCTURING_INDICES.end(), i) == PUNCTURING_INDICES.end()) {
            // this position is not punctured
            unpunctured_symbols[i] = *llr_ptr;
            llr_ptr++;
        } else {
            unpunctured_symbols[i] = 0.0f;
        }
    }    
    /*
    // projecting onto the codeword sphere
    float received_word_energy = utils::compute_vector_energy(unpunctured_symbols);
    float energy_normalize_factor = std::sqrt((float)N / received_word_energy);  // normalizing received message
    std::vector<float> projected_received_word(unpunctured_symbols.size(), 0.0);
    for (size_t i = 0; i < unpunctured_symbols.size(); i++) {
      projected_received_word[i] = unpunctured_symbols[i] * energy_normalize_factor;
    }
    // std::cout << "unpuncture the bits" << std::endl;
    // utils::print_double_vector(unpunctured_symbols);
    // std::cout << std::endl;

    MessageInformation mi_result = code.decode(projected_received_word, 4, PUNCTURING_INDICES, 0);
    */

    MessageInformation mi_result = code.decode(unpunctured_symbols, PUNCTURING_INDICES, 0);
    
    info_est = mi_result.message;
    int result = 1;
    return result;
}

