#include "lowRateListDecoder.h"

LowRateListDecoder::LowRateListDecoder(FeedForwardTrellis feedforwardTrellis, int listSize, int crcDegree, int crc, char code_type) {
  this->lowrate_nextStates    = feedforwardTrellis.getNextStates();
	this->lowrate_outputs       = feedforwardTrellis.getOutputs();
	this->lowrate_numStates     = feedforwardTrellis.getNumStates();
	this->lowrate_symbolLength  = feedforwardTrellis.getN();
	this->numForwardPaths       = lowrate_nextStates[0].size();
	this->listSize              = listSize;
	this->crcDegree             = crcDegree;
	this->crc                   = crc;
	this->code_type				= code_type;
	
	// int V = feedforwardTrellis.getV();
}

MessageInformation LowRateListDecoder::decode(std::vector<float> receivedMessage, std::vector<int> punctured_indices) {
	/** construct either ZT or TB trellis accoding to parameter passed into the constructor
	 * 
	 */
	if (this->code_type == 'Z') {
		return SSD_SLVD_ZT(receivedMessage);

	} else if (this->code_type == 'T') {
		return SSD_SLVD_TB(receivedMessage, punctured_indices);
	}
	throw std::invalid_argument("INVALID DECODING CHOICE!");
}

MessageInformation LowRateListDecoder::lowRateDecoding_MaxListsize(std::vector<float> receivedMessage, std::vector<int> punctured_indices){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis_Punctured(receivedMessage, punctured_indices);

	// start search
	MessageInformation output;
	//RBTree detourTree;
	MinHeap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	int TBPathsSearched = 0;
  
	while(numPathsSearched < this->listSize){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		float forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			float suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			float prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			float suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			float currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			float prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		
		previousPaths.push_back(path);

		std::vector<int> message = pathToMessage(path);
		std::vector<int> codeword = pathToCodeword(path);
		
		// one trellis decoding requires both a tb and crc check
		if(path[0] == path[lowrate_pathLength - 1] && crc::crc_check(message, crcDegree, crc) && numPathsSearched <= this->listSize){
			output.message = message;
			output.path = path;
		 	output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			output.TBListSize = TBPathsSearched + 1;
		 	return output;
		}

		numPathsSearched++;
		if(path[0] == path[lowrate_pathLength - 1])
			TBPathsSearched++;
	} // while(numPathsSearched < this->listSize)

	output.listSizeExceeded = true;
	return output;
}

// SSD-SLVD for tail-biting punctured convolutional codes 
MessageInformation LowRateListDecoder::SSD_SLVD_TB(std::vector<float> receivedMessage, std::vector<int> punctured_indices){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis_Punctured(receivedMessage, punctured_indices);

	// start search
	MessageInformation output;
	//RBTree detourTree;
	MinHeap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	int TBPathsSearched = 0;
	float metric_to_beat = INT_MAX;
	int cur_neighbor_msg_bit; 
    int cur_neighbor_cwd_bit; 
	std::vector<int> linearity_message(K, 0);
	int real_N = N + PUNCTURING_INDICES.size();
	std::vector<int> linearity_cwd(real_N, 0);
  
	while(numPathsSearched < this->listSize){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		float forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			float suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			float prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			float suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			float currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			float prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		
		previousPaths.push_back(path);

		std::vector<int> message = pathToMessage(path);
		std::vector<int> codeword = pathToCodeword(path);
		int ESD = path[0]^path[lowrate_pathLength - 1];
		int ED = crc::crc_remainder(message, crcDegree, crc);

		// return if SLVD finds ELF-TB codeword
		if(!ESD && !ED){	
			output.message = message;
			output.path = path;
		 	output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			output.TBListSize = TBPathsSearched + 1;
			output.SSD = false;
		 	return output;
		}
		// apply offset table
		// puncture codeword
		for (int b = 0; b < punctured_indices.size(); b++) {
            codeword[punctured_indices[b]] = 0;
        }
		int neighbor_idx = ESD*(1<<M) + ED;
        int numNeighbor = messageList[neighbor_idx].size()/K;
		for (int nei=0; nei<numNeighbor; nei ++){ 
			// get ELF-TB neighbor codeword and calculate metric
            float linearity_metric = 0.0;
			for (int ii = 0; ii < real_N; ii++) {
				cur_neighbor_cwd_bit = codewordList[neighbor_idx][nei*real_N + ii];
				linearity_cwd[ii] = (codeword[ii] * cur_neighbor_cwd_bit); // punctured bit remains 0
				linearity_metric += -1*((float)linearity_cwd[ii] * receivedMessage[ii]);
			}
			// compare metrics and store current best message and metric
			if (linearity_metric < metric_to_beat){
				metric_to_beat = linearity_metric;
				// get linearity message
				for (int jj = 0; jj < K; jj++) {
					cur_neighbor_msg_bit = messageList[neighbor_idx][nei*K + jj];
					linearity_message[jj] = message[jj] ^ cur_neighbor_msg_bit;
				}
			}
		}
		// if we reach max list size and still no ELF-TB identified by SLVD, return SSD result
		if (numPathsSearched == this->listSize - 1){
			output.message = linearity_message;
		 	output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			output.TBListSize = TBPathsSearched + 1;
			output.SSD = true;
		 	return output;
		}

		numPathsSearched++;
		if(path[0] == path[lowrate_pathLength - 1])
			TBPathsSearched++;
	} // while(numPathsSearched < this->listSize)
	output.SSD = false;
	output.listSizeExceeded = true;
	return output;
}

// SSD-SLVD for zero-terminated convolutional codes without puncturing
MessageInformation LowRateListDecoder::SSD_SLVD_ZT(std::vector<float> receivedMessage){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis_ZT(receivedMessage);

	// start search
	MessageInformation output;
	//RBTree detourTree;
	MinHeap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	float metric_to_beat = INT_MAX;
	int cur_neighbor_msg_bit; 
    int cur_neighbor_cwd_bit; 
	std::vector<int> linearity_message(K, 0);
	int real_N = N;
	std::vector<int> linearity_cwd(real_N, 0);
  
	while(numPathsSearched < this->listSize){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		float forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			float suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			float prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			float suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			float currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			float prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		
		previousPaths.push_back(path);

		std::vector<int> message = pathToMessage(path);
		std::vector<int> codeword = pathToCodeword(path);
		int ED = crc::crc_remainder(message, crcDegree, crc);

		// return if SLVD finds ELF-TB codeword
		// ZT only needs to check ELF difference
		if(!ED){	
			output.message = message;
			output.path = path;
		 	output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			output.SSD = false;
		 	return output;
		}
		// apply offset table
		int neighbor_idx = ED;
        int numNeighbor = messageList[neighbor_idx].size()/K;
		for (int nei=0; nei<numNeighbor; nei ++){ 
			// get ELF neighbor codeword and calculate metric
            float linearity_metric = 0.0;
			for (int ii = 0; ii < real_N; ii++) {
				cur_neighbor_cwd_bit = codewordList[neighbor_idx][nei*real_N + ii];
				linearity_cwd[ii] = (codeword[ii] * cur_neighbor_cwd_bit); // punctured bit remains 0
				linearity_metric += -1* ((float)linearity_cwd[ii] * receivedMessage[ii]);
			}
			// compare metrics and store current best message and metric
			if (linearity_metric < metric_to_beat){
				metric_to_beat = linearity_metric;
				// get linearity message
				for (int jj = 0; jj < K; jj++) {
					cur_neighbor_msg_bit = messageList[neighbor_idx][nei*K + jj];
					linearity_message[jj] = message[jj] ^ cur_neighbor_msg_bit;
				}
			}
		}
		// if we reach max list size and still no ELF-TB identified by SLVD, return SSD result
		if (numPathsSearched == this->listSize - 1){
			output.message = linearity_message;
		 	output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			output.SSD = true;
		 	return output;
		}

		numPathsSearched++;

	} // while(numPathsSearched < this->listSize)
	output.SSD = false;
	output.listSizeExceeded = true;
	return output;
}

// construct ZT trellis
std::vector<std::vector<LowRateListDecoder::cell>> LowRateListDecoder::constructLowRateTrellis_ZT(std::vector<float> receivedMessage){
	std::vector<std::vector<cell>> trellisInfo;
	lowrate_pathLength = (receivedMessage.size() / lowrate_symbolLength) + 1;

	trellisInfo = std::vector<std::vector<cell>>(lowrate_numStates, std::vector<cell>(lowrate_pathLength));

	// initialize only 0 as the starting states
	trellisInfo[0][0].pathMetric = 0;
	trellisInfo[0][0].init = true;
	
	// building the trellis
	for(int stage = 0; stage < lowrate_pathLength - V - 1; stage++){
		for(int currentState = 0; currentState < lowrate_numStates; currentState++){
			// if the state / stage is invalid, we move on
			if(!trellisInfo[currentState][stage].init)
				continue;

			// otherwise, we compute the relevent information
			for(int forwardPathIndex = 0; forwardPathIndex < numForwardPaths; forwardPathIndex++){
				// since our transitions correspond to symbols, the forwardPathIndex has no correlation 
				// beyond indexing the forward path

				int nextState = lowrate_nextStates[currentState][forwardPathIndex];
				
				// if the nextState is invalid, we move on
				if(nextState < 0)
					continue;
				
				float branchMetric = 0;
				std::vector<int> output_point = crc::get_point(lowrate_outputs[currentState][forwardPathIndex], lowrate_symbolLength);
				
				for(int i = 0; i < lowrate_symbolLength; i++){
					branchMetric += -1 * (receivedMessage[lowrate_symbolLength * stage + i] * (float)output_point[i]);
				}
				float totalPathMetric = branchMetric + trellisInfo[currentState][stage].pathMetric;
				
				// dealing with cases of uninitialized states, when the transition becomes the optimal father state, and suboptimal father state, in order
				if(!trellisInfo[nextState][stage + 1].init){
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
					trellisInfo[nextState][stage + 1].init = true;
				}
				else if(trellisInfo[nextState][stage + 1].pathMetric > totalPathMetric){
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = trellisInfo[nextState][stage + 1].pathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = trellisInfo[nextState][stage + 1].optimalFatherState;
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
				}
				else{
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = currentState;
				}
			}

		}
	}
	// ZT stage
	for(int stage = lowrate_pathLength - V - 1; stage < lowrate_pathLength - 1; stage++){
		for(int currentState = 0; currentState < lowrate_numStates; currentState++){
			// if the state / stage is invalid, we move on
			if(!trellisInfo[currentState][stage].init)
				continue;

			// zero terminating
			int forwardPathIndex = 0;
			// since our transitions correspond to symbols, the forwardPathIndex has no correlation 
			// beyond indexing the forward path

			int nextState = lowrate_nextStates[currentState][forwardPathIndex];
			
			// if the nextState is invalid, we move on
			if(nextState < 0)
				continue;
			
			float branchMetric = 0;
			std::vector<int> output_point = crc::get_point(lowrate_outputs[currentState][forwardPathIndex], lowrate_symbolLength);
			
			for(int i = 0; i < lowrate_symbolLength; i++){
				branchMetric += -1 * (receivedMessage[lowrate_symbolLength * stage + i] * (float)output_point[i]);
			}
			float totalPathMetric = branchMetric + trellisInfo[currentState][stage].pathMetric;
			
			// dealing with cases of uninitialized states, when the transition becomes the optimal father state, and suboptimal father state, in order
			if(!trellisInfo[nextState][stage + 1].init){
				trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
				trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
				trellisInfo[nextState][stage + 1].init = true;
			}
			else if(trellisInfo[nextState][stage + 1].pathMetric > totalPathMetric){
				trellisInfo[nextState][stage + 1].suboptimalPathMetric = trellisInfo[nextState][stage + 1].pathMetric;
				trellisInfo[nextState][stage + 1].suboptimalFatherState = trellisInfo[nextState][stage + 1].optimalFatherState;
				trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
				trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
			}
			else{
				trellisInfo[nextState][stage + 1].suboptimalPathMetric = totalPathMetric;
				trellisInfo[nextState][stage + 1].suboptimalFatherState = currentState;
			}

		}
	}
	return trellisInfo;
}

// minimizing the metric of product of -1* (received bit * decoded bit), instead of squared euclidean distance
std::vector<std::vector<LowRateListDecoder::cell>> LowRateListDecoder::constructLowRateTrellis(std::vector<float> receivedMessage){
	/* ---- Code Begins ---- */
	std::vector<std::vector<cell>> trellisInfo;
	lowrate_pathLength = (receivedMessage.size() / lowrate_symbolLength) + 1;

	trellisInfo = std::vector<std::vector<cell>>(lowrate_numStates, std::vector<cell>(lowrate_pathLength));

	// initializes all the valid starting states
	for(int i = 0; i < lowrate_numStates; i++){
		trellisInfo[i][0].pathMetric = 0;
		trellisInfo[i][0].init = true;
	}
	
	// building the trellis
	for(int stage = 0; stage < lowrate_pathLength - 1; stage++){
		for(int currentState = 0; currentState < lowrate_numStates; currentState++){
			// if the state / stage is invalid, we move on
			if(!trellisInfo[currentState][stage].init)
				continue;

			// otherwise, we compute the relevent information
			for(int forwardPathIndex = 0; forwardPathIndex < numForwardPaths; forwardPathIndex++){
				// since our transitions correspond to symbols, the forwardPathIndex has no correlation 
				// beyond indexing the forward path

				int nextState = lowrate_nextStates[currentState][forwardPathIndex];
				
				// if the nextState is invalid, we move on
				if(nextState < 0)
					continue;
				
				float branchMetric = 0;
				std::vector<int> output_point = crc::get_point(lowrate_outputs[currentState][forwardPathIndex], lowrate_symbolLength);
				
				for(int i = 0; i < lowrate_symbolLength; i++){
					branchMetric += -1 * (receivedMessage[lowrate_symbolLength * stage + i] * (float)output_point[i]);
				}
				
				float totalPathMetric = branchMetric + trellisInfo[currentState][stage].pathMetric;
				
				// dealing with cases of uninitialized states, when the transition becomes the optimal father state, and suboptimal father state, in order
				if(!trellisInfo[nextState][stage + 1].init){
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
					trellisInfo[nextState][stage + 1].init = true;
				}
				else if(trellisInfo[nextState][stage + 1].pathMetric > totalPathMetric){
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = trellisInfo[nextState][stage + 1].pathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = trellisInfo[nextState][stage + 1].optimalFatherState;
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
				}
				else{
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = currentState;
				}
			}

		}
	}
	return trellisInfo;
}

std::vector<std::vector<LowRateListDecoder::cell>> LowRateListDecoder::constructLowRateTrellis_Punctured(std::vector<float> receivedMessage, std::vector<int> punctured_indices){
	/* Constructs a trellis for a low rate code, with puncturing
		Args:
			receivedMessage (std::vector<float>): the received message
			punctured_indices (std::vector<int>): the indices of the punctured bits

		Returns:
			std::vector<std::vector<cell>>: the trellis
	*/

	/* ---- Code Begins ---- */
	std::vector<std::vector<cell>> trellisInfo;
	lowrate_pathLength = (receivedMessage.size() / lowrate_symbolLength) + 1;

	trellisInfo = std::vector<std::vector<cell>>(lowrate_numStates, std::vector<cell>(lowrate_pathLength));

	// initializes all the valid starting states
	for(int i = 0; i < lowrate_numStates; i++){
		trellisInfo[i][0].pathMetric = 0;
		trellisInfo[i][0].init = true;
	}
	
	// building the trellis
	for(int stage = 0; stage < lowrate_pathLength - 1; stage++){
		for(int currentState = 0; currentState < lowrate_numStates; currentState++){
			// if the state / stage is invalid, we move on
			if(!trellisInfo[currentState][stage].init)
				continue;

			// otherwise, we compute the relevent information
			for(int forwardPathIndex = 0; forwardPathIndex < numForwardPaths; forwardPathIndex++){
				// since our transitions correspond to symbols, the forwardPathIndex has no correlation 
				// beyond indexing the forward path

				int nextState = lowrate_nextStates[currentState][forwardPathIndex];
				
				// if the nextState is invalid, we move on
				if(nextState < 0)
					continue;
				
				float branchMetric = 0;
				std::vector<int> output_point = crc::get_point(lowrate_outputs[currentState][forwardPathIndex], lowrate_symbolLength);
				
				for(int i = 0; i < lowrate_symbolLength; i++){
					if (std::find(punctured_indices.begin(), punctured_indices.end(), lowrate_symbolLength * stage + i) != punctured_indices.end()){
						branchMetric += 0;
					} else {
						branchMetric += -1 * (receivedMessage[lowrate_symbolLength * stage + i] * (float)output_point[i]);
					}
				}
				
				float totalPathMetric = branchMetric + trellisInfo[currentState][stage].pathMetric;
				
				// dealing with cases of uninitialized states, when the transition becomes the optimal father state, and suboptimal father state, in order
				if(!trellisInfo[nextState][stage + 1].init){
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
					trellisInfo[nextState][stage + 1].init = true;
				}
				else if(trellisInfo[nextState][stage + 1].pathMetric > totalPathMetric){
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = trellisInfo[nextState][stage + 1].pathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = trellisInfo[nextState][stage + 1].optimalFatherState;
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
				}
				else{
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = currentState;
				}
			}

		}
	}
	return trellisInfo;
}

// converts a path through the tb trellis to the binary message it corresponds with
std::vector<int> LowRateListDecoder::pathToMessage(std::vector<int> path){
	std::vector<int> message;
	for(int pathIndex = 0; pathIndex < path.size() - 1; pathIndex++){
		for(int forwardPath = 0; forwardPath < numForwardPaths; forwardPath++){
			if(lowrate_nextStates[path[pathIndex]][forwardPath] == path[pathIndex + 1])
				message.push_back(forwardPath);
		}
	}
	return message;
}

// converts a path through the tb trellis to the BPSK it corresponds with
// currently does NOT puncture the codeword
std::vector<int> LowRateListDecoder::pathToCodeword(std::vector<int> path){
	std::vector<int> nopunc_codeword;
	for(int pathIndex = 0; pathIndex < path.size() - 1; pathIndex++){
		for(int forwardPath = 0; forwardPath < numForwardPaths; forwardPath++){
			if(lowrate_nextStates[path[pathIndex]][forwardPath] == path[pathIndex + 1]){
				std::vector<int> output_bin;
				crc::dec_to_binary(lowrate_outputs[path[pathIndex]][forwardPath], output_bin, lowrate_symbolLength);
				for (int outbit = 0; outbit < lowrate_symbolLength; outbit ++){
					nopunc_codeword.push_back(-2 * output_bin[outbit] + 1);
				}
			}
		}
	}

	return nopunc_codeword;
}