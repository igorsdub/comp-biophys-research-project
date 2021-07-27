#!/bin/bash
# Executes essential DDPT routines on PDB file
# DDPT must be installed and added to PATH

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <pdb-filepath> <results-dir> <GENENMM-flags>"
    echo "Example: $0 pdb/processed/1.pdb "-ca -het -c 8.0""
    exit 1
fi

echo "PDB filepath:"    $1
echo "Results dir:"     $2
echo "GENENMM flags:"   $3

# Shortcut to the files and binaries
ROOT=$PWD
PDB_FILEPATH=${ROOT}/$1
RESULTS_DIR=${ROOT}/$2
GENENMM_FLAGS=$3

# PDB filepath must be in a format PDB_ID.FORM_IDX.pdb
# Extract PDB filename 
FILENAME=$(basename "${PDB_FILEPATH}")
FILENAME="${FILENAME%%.*}"
echo "Filename:"          ${FILENAME}

# Create results and working dirs
# printf -v CUTOFF_PAD "%04.1f" $CUTOFF
WORK_DIR=$ROOT/tmp/working-${FILENAME}
mkdir -p $RESULTS_DIR $WORK_DIR 

# Copy auxilary files, if any are present,
# for -mass -ca -res, -ccust, -spcust and -fcust flags. 
# mv -f $ROOT/misc/{resmass.dat,cutoff.radius,res.force,fix.springs} -t $WORK_DIR 2> /dev/null

echo -n "Moving to working directory:"
pushd $WORK_DIR

echo "Running: GENENMM"
GENENMM -pdb $PDB_FILEPATH $GENENMM_FLAGS
echo

echo "Running: DIAGSTD"
DIAGSTD -i matrix.sdijf
echo

echo "Extract eigenvalues:"
grep "VECTOR" matrix.eigenfacs | awk '{ print $4 }' > eigenvalues
echo

echo "Copy matrix.eigenfacs file to results dir:"
mv matrix.eigenfacs ${RESULTS_DIR}
echo

echo "Moving out of working directory: "
popd
rm -rf $WORK_DIR
echo
