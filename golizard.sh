#!/bin/bash

# Exit on errors
set -e

# Define paths to Python scripts
SCRIPT_DIR=$(dirname "$0")
INPUT_CHECK="$SCRIPT_DIR/retrieval/check_input.py"
BLAST_API_SCRIPT="$SCRIPT_DIR/retrieval/get_blast_api.py"
BLAST_LOCAL_SCRIPT="$SCRIPT_DIR/retrieval/get_blast_local.py"
FS_API_SCRIPT="$SCRIPT_DIR/retrieval/get_fold_api.py"
ELM_SCRIPT="$SCRIPT_DIR/retrieval/get_elm.py"
MERGE_ID_SCRIPT="$SCRIPT_DIR/retrieval/merge_ids.py"
GO_SCRIPT="$SCRIPT_DIR/retrieval/get_go_terms.py"
MERGE_GO_SCRIPT="$SCRIPT_DIR/retrieval/merge_go.py"
REVIGO_SCRIPT="$SCRIPT_DIR/generate_results/get_revigo.py"

GO_OBO="$SCRIPT_DIR/meta/go-basic.obo"
ELM_TO_GO="$SCRIPT_DIR/retrieval/elm_goterms.tsv"

BLAST_DBS=("swissprot" "refseq_protein" "nr" 'genbank' 'gnomon' 'pdb')
BLAST_PROGS=("blastp" "tblastn")
FS_DBS=("afdb50" "afdb-swissprot" "afdb-proteome")
MAX_SEQ_LEN=400

# Function to show usage
usage() {
    echo "Usage: $0 --file <path_to_fasta_file> [--outputdir <output_directory>]"
    echo
    echo "Optional arguments:"
    echo
    echo "BLAST options:"
    echo "  --blast-dbs <database>          List of databases for remote BLAST (default: swissprot)"
    echo "                                  Valid options: ${BLAST_DBS[*]}"
    echo "  --blast-prog <program>          BLAST program (default: blastp)"
    echo "                                  Valid options: ${BLAST_PROGS[*]}"
    echo "  --blast-kMax <int>              Maximum rank of hits to consider (default: 10)"
    echo "  --blast-max_eval <float>        Maximum e-value (default: 1e-5)"
    echo "  --blast-min_ident <float>       Minimum identity percentage (default: 30.0)"
    echo "  --blast-cluster <float>         Clustering threshold (default: 0.9)"
    echo "  --blast-local <path>      Path to a local BLAST DB (no file extensions)."
    echo "                            If specified *alone*, only local BLAST runs."
    echo "                            If specified with --blast-db, both local and remote BLAST are run."
    echo "  --blast-res <path>        Directly input BLAST results (one protein ID per line) to skip BLAST"
    echo
    echo "FoldSeek options:"
    echo "  --fs-dbs <database>          List of FoldSeek databases (default: afdb50)"
    echo "                               Valid options: ${FS_DBS[*]}"
    echo "  --fs-max_eval <float>        Maximum e-value (default: 10)"
    echo "  --fs-min_ident <float>       Minimum identity percentage (default: 30.0)"
    echo "  --fs-kMax <int>              Maximum hits to retain (default: 10)"
    echo "  --fs-res <path>              Directly input FoldSeek results (one protein ID per line) to skip FoldSeek"
    echo
    echo "  --results <path>             Directly input BLAST and FoldSeek results (one protein ID per line) to skip querying"
    echo
    echo "GO Term options:"
    echo "  --GO <basic|plant>  Choose the GO ontology file to use (default: basic)"
    echo
    echo "Example usage:"
    echo "  $0 --file my_input.fasta --blast-local /path/to/blastdb/swissprot --subset-ids my_id_list.txt"
    echo "  $0 --file my_input.fasta --blast-local /path/to/blastdb/swissprot --subset-ids my_id_list.txt --blast-db swissprot"
    echo "  $0 --file my_input.fasta [no additional blast arguments => defaults to remote swissprot]"
    exit 1
}

# Parse arguments
if [[ $# -eq 0 ]]; then
    usage
fi

# Default parameters
BLAST_DB=()
BLAST_PROG="blastp"
BLAST_KMAX=10
BLAST_MAXEVAL=1e-5
BLAST_MINIDENTITY=30.0
BLAST_CLUSTER=0.9
BLAST_LOCAL_DB=()

FS_DB=()
FS_MAXEVAL=10
FS_MINIDENTITY=30.0
FS_KMAX=10
FS_CLUSTER=1.0

INPUT_FILE=""
OUTPUT="$SCRIPT_DIR/results"
BLAST_RESULTS=""
FS_API_RESULTS=""
COMBINED_RESULTS=""
ELM_RESULTS=""

ELM_RUN=false


# Read arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --file)
            INPUT_FILE="$2"
            shift 2
            ;;
        --outputdir)
            OUTPUT="$2"
            shift 2
            ;;
        --blast-dbs)
            shift
            while [[ $# -gt 0 && ! $1 == --* ]]; do
                BLAST_DB+=("$1")
                shift
            done
            ;;
        --blast-local)
            shift
            while [[ $# -gt 0 && ! $1 == --* ]]; do
                BLAST_LOCAL_DB+=("$1")
                shift
            done
            ;;
        --blast-prog)
            BLAST_PROG="$2"
            shift 2
            ;;
        --blast-kMax)
            BLAST_KMAX="$2"
            shift 2
            ;;
        --blast-max_eval)
            BLAST_MAXEVAL="$2"
            shift 2
            ;;
        --blast-min_ident)
            BLAST_MINIDENTITY="$2"
            shift 2
            ;;
        --blast-cluster)
            BLAST_CLUSTER="$2"
            shift 2
            ;;
        --blast-res)
            BLAST_RESULTS="$2"
            shift 2
            ;;
        --fs-dbs)
            shift
            while [[ $# -gt 0 && ! $1 == --* ]]; do
                FS_DB+=("$1")
                shift
            done
            ;;
        --fs-max_eval)
            FS_MAXEVAL="$2"
            shift 2
            ;;
        --fs-min_ident)
            FS_MINIDENTITY="$2"
            shift 2
            ;;
        --fs-kMax)
            FS_KMAX="$2"
            shift 2
            ;;
        --fs-res)
            FS_API_RESULTS="$2"
            shift 2
            ;;
        --results)
            COMBINED_RESULTS="$2"
            shift 2
            ;;
        --GO)
            case $2 in
                #full)
                #    GO_OBO="$SCRIPT_DIR/meta/go.obo"
                #    ;;
                basic)
                    GO_OBO="$SCRIPT_DIR/meta/go-basic.obo"
                    ;;
                plant)
                    GO_OBO="$SCRIPT_DIR/meta/goslim_plant.obo"
                    ;;
                *)
                    echo "[ERROR] Invalid GO option: $2"
                    usage
                    exit 1
                    ;;
            esac
            shift
            shift
            ;;
        *)
            usage
            ;;
    esac
done

if [[ ${#FS_DB[@]} -eq 0 ]]; then
    FS_DB=("afdb50")
fi

generate_subdir() {
    local base_dir="$1"
    local datetime
    local subdir

    datetime=$(date +"%Y-%m-%d_%H%M%S")
    echo "$base_dir/$datetime"
}

OUTPUT="${OUTPUT%/}"
OUTPUT_DIR=$(generate_subdir "$OUTPUT")

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"
echo

#### Step 1: Check input file ####
##########################################################################
if [[ -z "$COMBINED_RESULTS" ]]; then

    if [[ -z "$INPUT_FILE" ]]; then
        echo "[ERROR] Input file not provided"
        usage
    fi

    VALID_SEQ=$(python "$INPUT_CHECK" --file "$INPUT_FILE" --max_len $MAX_SEQ_LEN)
    EXIT_CODE=$?

    if [[ $EXIT_CODE -ne 0 ]]; then
        echo "[ERROR] FASTA validation failed"
        echo "$VALID_SEQ"
        exit 1
    fi

elif [[ -n "$COMBINED_RESULTS" && ! -f "$COMBINED_RESULTS" ]]; then
    echo "[ERROR] Provided BLAST and FoldSeek results file does not exist: $COMBINED_RESULTS"
    exit 1
elif [[ -n "$COMBINED_RESULTS" ]] && [[ ! -s "$COMBINED_RESULTS" ]]; then
    echo "[ERROR] Provided BLAST and FoldSeek results file is empty: $COMBINED_RESULTS"
    exit 1
fi

if [[ ! -z "$BLAST_RESULTS" ]]; then
    if [[ -n "$BLAST_RESULTS" && ! -s "$BLAST_RESULTS" ]]; then
        echo "[ERROR] Provided BLAST results file is empty: $BLAST_RESULTS"
        exit 1
    fi
fi

if [[ ! -z "$FS_API_RESULTS" ]]; then
    if [[ -n "$FS_API_RESULTS" && ! -s "$FS_API_RESULTS" ]]; then
        echo "[ERROR] Provided FoldSeek results file is empty: $FS_API_RESULTS"
        exit 1
    fi
fi
#### Step 2: Decide which BLAST steps to run ####
##########################################################################
# If no BLAST results are found, run ELM Domain Analysis
# Requirement:
# - If only blast-local is given => local only
# - If only blast-db is given => remote only
# - If both are given => run both
# - If none are given => default to remote SwissProt

if [[ -z "$BLAST_RESULTS" && -z "$COMBINED_RESULTS" ]]; then

    BLAST_API_RESULTS="$OUTPUT_DIR/blast_api_results.txt"
    BLAST_LOCAL_RESULTS="$OUTPUT_DIR/blast_local_results.txt"
    BLAST_RESULTS="$OUTPUT_DIR/blast_results.txt"

    LOCAL_ONLY=false
    API_ONLY=false
    BOTH_BLAST=false
    DEFAULT_REMOTE=false

    if [[ -n "$BLAST_LOCAL_DB" && -n "$BLAST_DB" ]]; then
        BOTH_BLAST=true
    elif [[ -n "$BLAST_LOCAL_DB" && -z "$BLAST_DB" ]]; then
        LOCAL_ONLY=true
    elif [[ -z "$BLAST_LOCAL_DB" && -n "$BLAST_DB" ]]; then
        API_ONLY=true
    else
        DEFAULT_REMOTE=true
    fi

    if [[ "$BOTH_BLAST" == true ]]; then

        echo "[INFO] Running BOTH local and remote BLAST..."

        # Local BLAST
        python "$BLAST_LOCAL_SCRIPT" \
            --fasta "$INPUT_FILE" \
            --db_paths "${BLAST_LOCAL_DB[@]}" \
            --prog "$BLAST_PROG" \
            --kMax "$BLAST_KMAX" \
            --max_eval "$BLAST_MAXEVAL" \
            --min_ident "$BLAST_MINIDENTITY" \
            --cluster "$BLAST_CLUSTER" \
            --output_file "$BLAST_LOCAL_RESULTS" || exit 1

        # Remote BLAST
        python "$BLAST_API_SCRIPT" \
            --fasta "$INPUT_FILE" \
            --dbs "${BLAST_DB[@]}" \
            --prog "$BLAST_PROG" \
            --kMax "$BLAST_KMAX" \
            --max_eval "$BLAST_MAXEVAL" \
            --min_ident "$BLAST_MINIDENTITY" \
            --cluster "$BLAST_CLUSTER" \
            --max_len $MAX_SEQ_LEN \
            --output_file "$BLAST_API_RESULTS" || exit 1


    elif [[ "$LOCAL_ONLY" == true ]]; then
        echo "[INFO] Running Local BLAST only..."
        python "$BLAST_LOCAL_SCRIPT" \
            --fasta "$INPUT_FILE" \
            --db_paths "${BLAST_LOCAL_DB[@]}" \
            --prog "$BLAST_PROG" \
            --kMax "$BLAST_KMAX" \
            --max_eval "$BLAST_MAXEVAL" \
            --min_ident "$BLAST_MINIDENTITY" \
            --cluster "$BLAST_CLUSTER" \
            --output_file "$BLAST_LOCAL_RESULTS" || exit 1

        rm -f "$BLAST_API_RESULTS"

    elif [[ "$API_ONLY" == true ]]; then
        echo "[INFO] Running Remote BLAST only..."

        python "$BLAST_API_SCRIPT" \
            --fasta "$INPUT_FILE" \
            --dbs "${BLAST_DB[@]}" \
            --prog "$BLAST_PROG" \
            --kMax "$BLAST_KMAX" \
            --max_eval "$BLAST_MAXEVAL" \
            --min_ident "$BLAST_MINIDENTITY" \
            --cluster "$BLAST_CLUSTER" \
            --max_len $MAX_SEQ_LEN \
            --output_file "$BLAST_API_RESULTS" || exit 1

        rm -f "$BLAST_LOCAL_RESULTS"

    else
        echo "[WARNING] No BLAST arguments provided. Running default remote BLAST with swissprot..."

        python "$BLAST_API_SCRIPT" \
            --fasta "$INPUT_FILE" \
            --dbs "swissprot" \
            --prog "$BLAST_PROG" \
            --kMax "$BLAST_KMAX" \
            --max_eval "$BLAST_MAXEVAL" \
            --min_ident "$BLAST_MINIDENTITY" \
            --cluster "$BLAST_CLUSTER" \
            --max_len $MAX_SEQ_LEN \
            --output_file "$BLAST_API_RESULTS" || exit 1

        rm -f "$BLAST_LOCAL_RESULTS"

    fi
#### Step 3: Merge local & remote BLAST results (if both exist) ####
# If no BLAST results are found, run ELM Domain Analysis #
##########################################################################

    python "$MERGE_ID_SCRIPT" \
        --local_file "$BLAST_LOCAL_RESULTS" \
        --api_file "$BLAST_API_RESULTS" \
        --output_file "$BLAST_RESULTS" || exit 1

    if [[ ! -s "$BLAST_RESULTS" ]]; then
        echo "[WARNING] No BLAST hits found. Running ELM domain analysis..."
        echo "--------------------------------------------------------------------"
        ELM_JSON="$OUTPUT_DIR/elm_go_terms.json"
        ELM_CSV="$OUTPUT_DIR/elm_go_terms.csv"

        python "$ELM_SCRIPT" \
            --fasta "$INPUT_FILE" \
            --go_map "$ELM_TO_GO" \
            --output_json "$ELM_JSON" \
            --output_csv "$ELM_CSV" || exit 1

        if [[ -f "$ELM_CSV" && -s "$ELM_CSV" ]]; then
            ELM_RUN=true
            echo "[INFO] ELM dict results saved => $ELM_JSON"
            echo "[INFO] ELM table results saved => $ELM_CSV"
            echo "--------------------------------------------------------------------"
        fi
    else
        echo "[INFO] BLAST hits found. ELM not required."
        echo "--------------------------------------------------------------------"
    fi

elif [[ -n "$BLAST_RESULTS" && ! -s "$BLAST_RESULTS" ]]; then
    echo "[ERROR] Provided BLAST results file is empty: $BLAST_RESULTS"
    exit 1

elif [[ -z "$COMBINED_RESULTS" ]]; then
    echo "Loading BLAST results from $BLAST_RESULTS..."
    echo "--------------------------------------------------------------------"
fi

#### Step 4: Run FoldSeek ####
##########################################################################
if [[ -z "$FS_API_RESULTS" && -z $COMBINED_RESULTS ]]; then

    echo "[INFO] Running FoldSeek..."
    FS_API_RESULTS="$OUTPUT_DIR/fold_api_results.txt"
    python "$FS_API_SCRIPT" \
            --fasta "$INPUT_FILE" \
            --dbs "${FS_DB[@]}" \
            --max_eval "$FS_MAXEVAL" \
            --min_ident "$FS_MINIDENTITY" \
            --kMax "$FS_KMAX" \
            --max_len $MAX_SEQ_LEN \
            --output_file "$FS_API_RESULTS" || exit 1

    if [[ -n "$FS_API_RESULTS" && -s "$FS_API_RESULTS" ]]; then
        echo "[INFO] FoldSeek results saved => $FS_API_RESULTS"
    fi

elif [[ -n "$FS_API_RESULTS" && ! -s "$FS_API_RESULTS" ]]; then
    echo "[ERROR] Provided FoldSeek results file is empty: $FS_API_RESULTS"
    exit 1

elif [[ -z "$COMBINED_RESULTS" ]]; then
    echo "Loading FoldSeek results from $FS_API_RESULTS..."
fi
echo "--------------------------------------------------------------------"

#### Step 5: Combine results ####
##########################################################################
if [[ -z "$COMBINED_RESULTS" ]]; then
    COMBINED_RESULTS="$OUTPUT_DIR/combined_results.txt"

    # 1️ BLAST present (BLAST + FoldSeek combined)
    if [[ "$ELM_RUN" == false ]]; then
        echo "[INFO] Combining BLAST and FoldSeek IDs..."

        # Not local but Blast | Not api but FoldSeek
        python "$MERGE_ID_SCRIPT" \
            --local_file "$BLAST_RESULTS" \
            --api_file "$FS_API_RESULTS" \
            --output_file "$COMBINED_RESULTS" || exit 1

    # 2 No BLAST (but ELM and FoldSeek present)
    elif [[ "$ELM_RUN" == true && -s "$FS_API_RESULTS"  ]]; then
        echo "[INFO] Combining ELM with FoldSeek IDs..."
        cp "$FS_API_RESULTS" "$COMBINED_RESULTS"

    # 3 Only ELM results
    elif [[ "$ELM_RUN" == true ]]; then
        echo "[INFO] No FoldSeek hits. Proceeding with ELM GO terms only."

    else
        echo "[INFO] No results from BLAST, FoldSeek and ELM found, exiting..."
        exit 1
    fi

elif [[ -n "$COMBINED_RESULTS" ]] && [[ ! -s "$COMBINED_RESULTS" ]]; then
    echo "[ERROR] Provided BLAST and FoldSeek results file is empty: $COMBINED_RESULTS"
    exit 1
else
    echo "Loading BLAST and FoldSeek results from $COMBINED_RESULTS"
fi
echo "--------------------------------------------------------------------"

#### Step 6: Get GO terms ####
##########################################################################
GO_DICT="$OUTPUT_DIR/go_terms.json"
GO_CSV="$OUTPUT_DIR/go_terms.csv"

# === Preliminary Check ===
if [[ ( ! -f "$ELM_JSON" || ! -s "$ELM_JSON" ) && ( ! -f "$COMBINED_RESULTS" || ! -s "$COMBINED_RESULTS" ) ]]; then
    echo "[ERROR] Neither ELM, BLAST nor FoldSeek results exist or are non-empty."
    echo "Exiting..."
    exit 1
fi

if [[ "$ELM_RUN" == false ]]; then
    echo "[INFO] Fetching GO terms for Protein IDs..."
    python "$GO_SCRIPT" \
        --file "$COMBINED_RESULTS" \
        --output_json "$GO_DICT" \
        --output_csv "$GO_CSV" || exit 1

elif [[ "$ELM_RUN" == true && -s "$COMBINED_RESULTS" ]]; then # FS results copied to COMBINED_RESULTS
        echo "[INFO] Fetching GO terms for FoldSeek IDs (ELM already mapped)..."
        FS_GO_DICT="$OUTPUT_DIR/fs_go_terms.json"
        FS_GO_CSV="$OUTPUT_DIRP/fs_go_terms.csv"

        python "$GO_SCRIPT" \
            --file "$COMBINED_RESULTS" \
            --output_json "$FS_GO_DICT" \
            --output_csv "$FS_GO_CSV" || exit 1
        python "$MERGE_GO_SCRIPT" \
            --elm "$ELM_RESULTS" \
            --fs "$FS_GO_TERMS" \
            --output "$GO_TERMS"
else
    cp "$ELM_RESULTS" "$GO_DICT"
fi

if [[ ! -f "$GO_DICT" || ! -s "$GO_DICT" ]]; then
    echo "[ERROR] GO terms not found or file is empty."
    exit 1
fi
echo "[INFO] GO terms saved => $GO_DICT"

#### Step 7: Graph GO terms ####
##########################################################################

echo "--------------------------------------------------------------------"
echo "[INFO] Creating GO graph using GOLizzard..."

python "$REVIGO_SCRIPT" \
    --uniprot_to_go "$GO_DICT" \
    --go_graph "$GO_OBO" \
    --output_dir "$OUTPUT_DIR" || exit 1


#----------------------------------
if [[ ! -d "$DOTFILES_DIR" || ! -s "$DOTFILES_DIR" ]]; then
    echo "[ERROR] GO graph not created or directory is empty."
    exit 1
fi

if [ "$(ls -A "$DOTFILES_DIR" 2>/dev/null)" ]; then
    echo "[INFO] GO graphs saved => $DOTFILES_DIR"
else
    echo "[ERROR] GO graph generation failed. Check the logs."
    exit 1
fi


# use pathlib import Path in path checks
# MAKE IT WORK ON WINDOWS -> installation guide
# LOCAL DATABASES (FOLD)
# ADD PDB MAPPING -> Use SIFTS or PDBe API to map PDB IDs to UniProt
# map_to_uniprot.py -> EMBL-GenBank-DDBJ not working
#./golizard.sh --file tst.fasta --blast-local /home/pospim/Desktop/work/GOLizard/blastdb/swissprot --blast-min_ident 99 --fs-min_ident 99
