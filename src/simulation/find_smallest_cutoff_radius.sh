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
echo
# Shortcut to the files and binaries
ROOT=$PWD
PDB_FILEPATH=$ROOT/$1
RESULTS_DIR=$ROOT/$2
GENENMM_FLAGS=$3

# PDB filepath must be in a format PDB_ID.PDB_FORM.pdb
# Extract PDB filename and ID name
PDB_ID=$(echo "$1" | sed 's/.*\(....\)\.[[:digit:]]*\.pdb$/\1/')
echo "PDB ID:"          $PDB_ID
PDB_FORM=$(echo "$1" | sed 's/.*\.\([[:digit:]]*\)\.pdb$/\1/')
echo "PDB form:"        $PDB_FORM


# Create results and working dirs
# printf -v CUTOFF_PAD "%04.1f" $CUTOFF
WORK_DIR=$ROOT/tmp/working-$PDB_ID-$PDB_FORM
mkdir -p $RESULTS_DIR $WORK_DIR

echo -n "Moving to working directory:"
pushd $WORK_DIR
echo

echo "Running: SPACING"
SPACING -pdb $PDB_FILEPATH -ca
echo

echo "Running: GENENMM"
GENENMM -pdb $PDB_FILEPATH $GENENMM_FLAGS
echo

echo "Running: DIAGSTD"
DIAGSTD
echo

echo "Extract eigenvalues:"
grep "VECTOR" matrix.eigenfacs > matrix.eigenvals
awk '{ print $4 }' matrix.eigenvals > eigenvalues


echo "Copy eigenvalues to results dir:"
cp eigenvalues dist.dat ${RESULTS_DIR}

echo "Moving out of working directory:"
echo
popd
rm -rf $WORK_DIR
