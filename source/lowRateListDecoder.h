
#ifndef LOWRATELISTDECODER_H
#define LOWRATELISTDECODER_H

#include <climits>
#include <iostream>

#include "FeedForwardTrellis.h"
#include "minHeap.h"
#include "tbcc_types.h"
#include "tbcc_namespace.h"
#include "consts.h"

class LowRateListDecoder{
public:
	LowRateListDecoder(FeedForwardTrellis FT, int listSize, int crcDegree, int crc, char code_type);
	MessageInformation lowRateDecoding_MaxListsize(std::vector<float> receivedMessage, std::vector<int> punctured_indices);
	MessageInformation SSD_SLVD_TB(std::vector<float> receivedMessage, std::vector<int> punctured_indices);
	MessageInformation SSD_SLVD_ZT(std::vector<float> receivedMessage);
	MessageInformation decode(std::vector<float> receivedMessage, std::vector<int> punctured_indices);


private:
	int numForwardPaths;
	int listSize;
	int crcDegree;
	int crc;
	int n;
	char code_type;

	std::vector<std::vector<int>> lowrate_nextStates;
	std::vector<std::vector<int>> lowrate_outputs;
	std::vector<std::vector<int>> neighboring_cwds; // ${listSize} x 516 matrix
	std::vector<std::vector<int>> neighboring_msgs;  // ${listSize} x 43 matrix
	std::vector<std::vector<int>> path_ie_state;
	int lowrate_numStates;
	int lowrate_symbolLength;
	int lowrate_pathLength;

	struct cell {
		int optimalFatherState = -1;
		int suboptimalFatherState = -1;
		float pathMetric = INT_MAX;
		float suboptimalPathMetric = INT_MAX;
		bool init = false;
	};

  std::vector<int> pathToMessage(std::vector<int>); 
  std::vector<int> pathToMessage_ZT(std::vector<int> path);
  std::vector<int> pathToCodeword(std::vector<int>); 
  void readNeighborList(std::string path);
  std::vector<std::vector<cell>> constructLowRateTrellis(std::vector<float> receivedMessage);
  std::vector<std::vector<cell>> constructLowRateTrellis_ZT(std::vector<float> receivedMessage);
  std::vector<std::vector<cell>> constructLowRateTrellis_Punctured(std::vector<float> receivedMessage, std::vector<int> punctured_indices);
};


#endif
