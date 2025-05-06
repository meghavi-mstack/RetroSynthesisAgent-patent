#!/bin/bash

# This script demonstrates the three different retrieval modes for RetroSynthesisAgent

# Check if a material name was provided
if [ -z "$1" ]; then
  echo "Usage: $0 <material_name> [num_results]"
  echo "Example: $0 aspirin 10"
  exit 1
fi

MATERIAL=$1
NUM_RESULTS=${2:-10}  # Default to 10 if not provided

echo "===== Mode 1: Patent-Patent (Patents for both initial and expansion) ====="
python3 main.py --material "$MATERIAL" --num_results $NUM_RESULTS --alignment True --expansion True --filtration False --retrieval_mode patent-patent

echo "===== Mode 2: Paper-Paper (Academic papers for both initial and expansion) ====="
python3 main.py --material "$MATERIAL" --num_results $NUM_RESULTS --alignment True --expansion True --filtration False --retrieval_mode paper-paper

echo "===== Mode 3: Both-Both (Both patents and papers for initial and expansion) ====="
python3 main.py --material "$MATERIAL" --num_results $NUM_RESULTS --alignment True --expansion True --filtration False --retrieval_mode both-both

echo "All retrieval modes completed!"
