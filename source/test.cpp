#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <algorithm>
#include "conv.h"

// Function declarations for tests
// void test_no_error(conv &code, int verbose = 0);
// void test_single_error(conv &code, float llr_mag, int verbose = 0);
// int test_gaussian_noise(conv &code, float esno, int verbose = 0);
// void test_tbcc_encoder(conv &code, int verbose = 0);

// void test_no_error(conv &code, int verbose) {
//     intvec info(N, 0); // Initialize info bits to zero
//     intvec cw(N);
//     fltvec llr(N);
//     intvec cw_est(N, 0);
//     fltvec llr_est(N);

//     if (verbose) {
//         std::cout << "Running Test No Error..." << std::endl;
//         std::cout << "Initial info bits: ";
//         for (const auto &bit : info) std::cout << bit << " ";
//         std::cout << std::endl;
//         std::cout << "Encoding info bits..." << std::endl;
//     }

//     // Encode
//     code.encode(info, cw);
//     for (size_t i = 0; i < cw.size(); ++i) {
//        llr[i] = (cw[i] == 0 ? 1.0f : -1.0f);
//     }

//     if (verbose) {
//         std::cout << "Encoded codeword: ";
//         for (const auto &bit : cw) std::cout << bit << " ";
//         std::cout << std::endl;
//     }

//     // Decode
//     MessageInformation result = code.decode(llr, 0.87, PUNCTURING_INDICES, 1);


//     if (result == 1) {
//        std::cout << "Test No Error: Passed" << std::endl;
//     } else {
//         std::cout << "Test No Error: Failed" << std::endl;
//     }
// }

// void test_single_error(conv &code, float llr_mag, int verbose) {
//     intvec info(code.n_cols, 0); // Initialize info bits to zero
//     intvec cw(code.n_cols);
//     fltvec llr(code.n_cols);
//     fltvec llr_out(code.n_cols);
//     intvec cw_est(code.n_cols);

//     code.encode(info, cw);

//     if (verbose) {
//         std::cout << "Running Test Single Error..." << std::endl;
//         std::cout << "Initial info bits: ";
//         for (const auto &bit : info) std::cout << bit << " ";
//         std::cout << std::endl;

//         std::cout << "Encoded codeword: ";
//         for (const auto &bit : cw) std::cout << bit << " ";
//         std::cout << std::endl;

//         std::fill(llr.begin(), llr.end(), llr_mag);
//         llr[0] =  -llr[0];
//         std::cout << "Decoding..." << std::endl;
//     }

//     int result = code.decode(llr, 20, llr_out, verbose);

//     if (verbose) {
//         std::cout << "LLR output from decoder: ";
//         for (const auto &llr_value : llr_out) std::cout << llr_value << " ";
//         std::cout << std::endl;
//     }

//     if (result == 1) {
//          std::cout << "Test Single Error: Passed" << std::endl;
//     } else {
//          std::cout << "Test Single Error: Failed" << std::endl;
//     }
// }

// int test_gaussian_noise(conv &code, float esno, int verbose) {
//     // Setup
//     intvec info(N, 0); // Initialize info bits to zero
//     intvec cw;
//     fltvec llr(N);
//     fltvec llr_out(N, 0.0f);
//     intvec cw_est(N);

//     // Setup binary RNG and encoding
//     std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
//     std::uniform_int_distribution<int> rand_bit(0, 1);
//     code.encode(info, cw);
//     std::cout << "cw size" << cw.size() << std::endl;
 
//     // Generate LLRs
//     std::normal_distribution<float> distribution(4*esno, std::sqrt(8*esno));
//     for (size_t i = 0; i < cw.size(); ++i) {
//         llr[i] = (cw[i] == 0 ? 1.0f : -1.0f) * distribution(generator);
//     }

//     if (verbose) {
//         std::cout << "Running Test Gaussian Noise..." << std::endl;
//         std::cout << "Initial info bits: ";
//         for (const auto &bit : info) std::cout << bit << " ";
//         std::cout << std::endl;

//         std::cout << "Encoded codeword: ";
//         for (const auto &bit : cw) std::cout << bit << " ";
//         std::cout << std::endl;

//         std::cout << "LLRs with Gaussian noise: ";
//         for (const auto &value : llr) std::cout << value << " ";
//         std::cout << std::endl;
//     }

//     int result = code.decode(llr, 0.87, PUNCTURING_INDICES, verbose);

//     if (verbose) {
//         std::cout << "LLR output from decoder: ";
//         for (const auto &value : llr_out) {
//             if (value <= 0.0f) {
//                 std::cout << "XXX" << value << " XXX ";
//             } else
//                 std::cout << value << " ";
//         }
//         std::cout << std::endl;
 
//         if (result == 1) {
//             std::cout << "Test Gaussian Noise: Passed" << std::endl;
//         } else {
//             std::cout << "Test Gaussian Noise: Failed" << std::endl;
//         }
//     }

//     return result;
//     // return 0;
// }


// // Test conv::create_encoder and conv::encode
// void test_tbcc_encode(conv &code, int verbose) {

//     // Call conv::create_encoder to create the parity generator
//     code.create_encoder(verbose);


//     // Setup binary RNG
//     std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
//     std::uniform_int_distribution<int> distribution(0, 1);

//     // Generate random information bit string and call conv::encode
//     intvec info(k);
//     std::generate(info.begin(), info.end(), [&]() { return distribution(generator); });
//     intvec cw(n);
//     code.encode(info, cw);
// }


int main(int argc, char* argv[])
{
    return 0;
    conv code;
    
    // Test encoder
    // test_tbcc_encode(code, 1);

    // Generate short conv code

    // Run test functions
    // test_tbcc_encode(code, 0);
    // test_no_error(code, 0); 
    // test_single_error(code, 3.0f, 0);
    // test_gaussian_noise(code, 0.8, 1); // Example ESNO value

    // // Test encoder
    // test_tbcc_encode(code, 0);

    // // Test single error
    // test_single_error(code, 3.0f, 1);

    // // Test Gaussian noise
    // test_gaussian_noise(code, 0.72, 0);

    // // Test gaussian noise
    // int count = 0;
    // for (int i=0; i<100; ++i) {
    //     if (!test_gaussian_noise(code, 0.72, 0))
    //         count++;
    // }
    // std::cout << count << " errors out of 100 trials.\n"; 
}

