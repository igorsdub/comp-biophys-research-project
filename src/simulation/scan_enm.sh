#!/bin/bash
# Executes essential DDPT routines on PDB file
# DDPT must be installed and added to PATH

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <pdb-filepath> <results-filepath> <flags>"
    echo "Example: $0 pdb/processed/1c3b.1.pdb "-ca -het -c 8.0""
    exit 1
fi

echo "PDB filepath:"    $1
echo "Results dir:"     $2
echo "GENENMM flags:"   $3

# Shortcut to the files and binaries
ROOT=$PWD
PDB_FILEPATH=$ROOT/$1
RESULTS_DIR=$ROOT/$2
GENENMM_FLAGS=$3

# PDB filepath must be in a format PDB_ID.PDB_FORM.pdb
# Extract PDB filename and ID name
PDB_FORM=$(echo "$1" | sed 's/.*\.\([[:digit:]]*\)\.pdb$/\1/')
echo "PDB form:"        $PDB_FORM
PDB_ID=$(echo "$1" | sed 's/.*\(....\)\.[[:digit:]]*\.pdb$/\1/')
echo "PDB ID:"          $PDB_ID

# Create results and working dirs
# printf -v CUTOFF_PAD "%04.1f" $CUTOFF
WORK_DIR=$ROOT/tmp/working-$PDB_ID-$PDB_FORM
mkdir -p $RESULTS_DIR $WORK_DIR

echo -n "Moving to working directory:"
pushd $WORK_DIR

echo "Running: SPACING"
SPACING -pdb $PDB_FILEPATH -ca

echo "Running: GENENMM"
GENENMM -pdb $PDB_FILEPATH $GENENMM_FLAGS

echo "Running: DIAGSTD"
DIAGSTD

echo "Extract eigenvalues:"
grep "VECTOR" matrix.eigenfacs > matrix.eigenvals
awk '{ print $4 }' matrix.eigenvals > eigenvalues
NO_MODES=$(awk 'END{print NR}' matrix.eigenvals)

echo "Running: FREQEN"
FREQEN -s 1 -e $NO_MODES
# mv mode.energy mode.m$NO_MODES.energy

echo "Running: PROJECT"
PROJECT -pdb CAonly.pdb -s 1 -e 106 -scale 1

for MODE_NUM in 1 5 10 25 50 75 100
do
    echo "Running: CROSCOR and RMSCOL for $MODE_NUM non-trivial mode/s"
    CROSCOR -i matrix.eigenfacs -s 7 -e "$(($MODE_NUM+6))"
    RMSCOL -i matrix.eigenfacs -s 7 -e "$(($MODE_NUM+6))"
    printf -v MODE_NUM_PAD "%03d" $MODE_NUM
    mv crosscor.dat crosscor.m${MODE_NUM_PAD}.dat
    mv mode.bfactors mode.m${MODE_NUM_PAD}.bfactors 
    echo
done

# Write format flags to save files
CUTOFF=$(echo $GENENMM_FLAGS | sed 's/.*-c\([0-9]*\.[0-9]*\)/\1/')
PRINTF_STATEMENT=$(echo $GENENMM_FLAGS | sed "s/ //g" | sed "s/\(.*-c\).*/\1%05.2f/g")
printf -v FORMATTED_GENENMM_FLAGS $PRINTF_STATEMENT $CUTOFF

echo "Copy all simulation results to results dir:"
for file in *; do 
    cp -- "$file" "${RESULTS_DIR}/${PDB_ID}.${PDB_FORM}.${FORMATTED_GENENMM_FLAGS}.${file}"
done

echo "Moving out of working directory: "
popd
rm -rf $WORK_DIR
