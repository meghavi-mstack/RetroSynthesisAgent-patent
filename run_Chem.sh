#!/bin/bash

# Array of chemicals to process
chemicals=(
#  "Cubane"
#  "Bicyclo[1.1.1]pentane"
#  "Dodecahedrane"
  "Azulene"
#  "Cubylamine"
#  "Cyclobutadiene"
#  "Bicyclo[2.2.2]octane"
#  "aspirin"
#  "Tropone"
#  "Cyclopropenone"
#  "Baking soda"
#  "Ethanol"
#  "Sodium bicarbonate"
#  "Calcium carbonate"
#  "Chloroform"
#  "Tetrafluoroborate"
#"1-Ethylpiperazine"
#"Ethyl acetate"
#"Bis(pinacolato)diboron"
#"Ammonium molybdate"
#"Aminoacetaldehyde dimethyl acetal"

)

# Settings for all runs
NUM_RESULTS=16
ALIGNMENT=True
EXPANSION=True
FILTRATION=False
RETRIEVAL_MODE="both-both"

# Process each chemical
for chemical in "${chemicals[@]}"; do
  echo "Processing: $chemical"
  echo "========================================"

  # Run the main.py script with the current chemical
  python3 main.py --material "$chemical" --num_results $NUM_RESULTS --alignment $ALIGNMENT --expansion $EXPANSION --filtration $FILTRATION --retrieval_mode $RETRIEVAL_MODE

  echo "Completed: $chemical"
  echo "========================================"
  echo ""
done

echo "All chemicals processed successfully!"